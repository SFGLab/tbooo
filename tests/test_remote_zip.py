"""RemoteZip is exercised against a real archive served from a local HTTP server.

A local server is used rather than mocks because the whole point of the module is
correct HTTP range handling and ZIP offset arithmetic; mocking the transport would
test neither.
"""

from __future__ import annotations

import functools
import http.server
import random
import threading
import zipfile

import pytest

from tbooo.download.remote_zip import RemoteZip


class _RangeHandler(http.server.SimpleHTTPRequestHandler):
    """Minimal Range support; stdlib's handler has none, and Zenodo's does."""

    def do_GET(self):  # noqa: N802 - stdlib naming
        header = self.headers.get("Range")
        if not header:
            return super().do_GET()

        path = self.translate_path(self.path)
        with open(path, "rb") as handle:
            body = handle.read()
        start_text, _, end_text = header.partition("=")[2].partition("-")
        start = int(start_text) if start_text else len(body) - int(end_text)
        end = int(end_text) if end_text and start_text else len(body) - 1
        chunk = body[start:end + 1]

        self.send_response(206)
        self.send_header("Content-Type", "application/octet-stream")
        self.send_header("Content-Length", str(len(chunk)))
        self.send_header("Content-Range", f"bytes {start}-{end}/{len(body)}")
        self.end_headers()
        self.wfile.write(chunk)

    def log_message(self, *args):
        pass


@pytest.fixture(scope="module")
def served_zip(tmp_path_factory):
    """Serve a directory containing a small ZIP over HTTP; yield its URL."""
    root = tmp_path_factory.mktemp("zipserve")
    archive = root / "sample.zip"
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("dir/small.txt", "hello world\n")
        # Compressible payload large enough that deflate genuinely engages.
        zf.writestr("dir/big.tsv", "col_a\tcol_b\n" + "1\t2\n" * 50_000)
        zf.writestr("stored.bin", b"\x00\x01\x02\x03", compress_type=zipfile.ZIP_STORED)
        # Incompressible bulk, so the archive is much larger than the 64 KB tail
        # scan. Without this the selective-read test is vacuous.
        zf.writestr("bulk.bin", random.Random(0).randbytes(4_000_000),
                    compress_type=zipfile.ZIP_STORED)

    handler = functools.partial(_RangeHandler, directory=str(root))
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{server.server_address[1]}/sample.zip"
    server.shutdown()


@pytest.fixture(scope="module")
def served_zip_no_ranges(tmp_path_factory):
    """Same archive from a server that ignores Range headers."""
    root = tmp_path_factory.mktemp("norange")
    with zipfile.ZipFile(root / "sample.zip", "w") as zf:
        zf.writestr("a.txt", "x")
    handler = functools.partial(http.server.SimpleHTTPRequestHandler,
                                directory=str(root))
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    yield f"http://127.0.0.1:{server.server_address[1]}/sample.zip"
    server.shutdown()


def test_raises_when_the_server_ignores_range_requests(served_zip_no_ranges):
    """Slicing a full 200 response as a range would give wrong bytes silently."""
    with pytest.raises(OSError, match="ignored Range"):
        RemoteZip(served_zip_no_ranges).entries


def test_lists_every_member(served_zip):
    names = set(RemoteZip(served_zip).entries)
    assert names == {"dir/small.txt", "dir/big.tsv", "stored.bin", "bulk.bin"}


def test_reads_a_deflated_member(served_zip):
    assert RemoteZip(served_zip).read("dir/small.txt") == b"hello world\n"


def test_reads_a_large_deflated_member_intact(served_zip):
    data = RemoteZip(served_zip).read("dir/big.tsv")
    assert data.startswith(b"col_a\tcol_b\n")
    assert data.count(b"\n") == 50_001


def test_reads_a_stored_member(served_zip):
    """Stored entries have no deflate stream and must not be inflated."""
    assert RemoteZip(served_zip).read("stored.bin") == b"\x00\x01\x02\x03"


def test_entry_sizes_match_the_archive(served_zip):
    entry = RemoteZip(served_zip).entries["dir/big.tsv"]
    assert entry.uncompressed_size == len("col_a\tcol_b\n" + "1\t2\n" * 50_000)
    assert entry.compressed_size < entry.uncompressed_size


def test_missing_member_raises_with_a_useful_message(served_zip):
    with pytest.raises(KeyError, match="members available"):
        RemoteZip(served_zip).read("does/not/exist.txt")


def test_extract_writes_the_file(served_zip, tmp_path):
    dest = RemoteZip(served_zip).extract("dir/small.txt", tmp_path / "out" / "s.txt")
    assert dest.read_text() == "hello world\n"


def test_only_the_needed_bytes_are_fetched(served_zip, monkeypatch):
    """The entire point: reading one member must not download the whole archive."""
    zf = RemoteZip(served_zip)
    total = zf._total
    fetched = []
    original = zf._range

    def counting(start, length):
        fetched.append(length)
        return original(start, length)

    monkeypatch.setattr(zf, "_range", counting)
    zf.read("dir/small.txt")
    # Reading a 12-byte member from a 4 MB archive must not scale with the
    # archive: only the tail scan, central directory and member are fetched.
    assert sum(fetched) < total / 10, f"fetched {sum(fetched)} of {total} bytes"
