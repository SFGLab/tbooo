"""Read selected members from a remote ZIP archive over HTTP range requests.

MAGE distributes its processed data as a single 59.5 GB Zip64 archive on Zenodo,
but the files needed for population-structure work total ~165 MB. Downloading the
whole archive to use 0.3% of it is not acceptable, so this module parses the
central directory remotely and fetches only the requested members.

Requires a server that honours `Range` requests; Zenodo does.
"""

from __future__ import annotations

import io
import struct
import urllib.request
import zlib
from dataclasses import dataclass
from pathlib import Path

_EOCD64_LOCATOR = b"PK\x06\x07"
_EOCD64 = b"PK\x06\x06"
_EOCD = b"PK\x05\x06"
_CENTRAL_HEADER = b"PK\x01\x02"
_ZIP64_EXTRA_TAG = 0x0001
_SATURATED = 0xFFFFFFFF


@dataclass(frozen=True)
class ZipEntry:
    """One member of a remote archive."""

    name: str
    header_offset: int
    compressed_size: int
    uncompressed_size: int
    compress_type: int


def _parse_zip64_extra(extra: bytes, uncompressed: int, compressed: int,
                       header_offset: int) -> tuple[int, int, int]:
    """Pull real sizes/offset out of the Zip64 extra field.

    Values appear in a fixed order but only for the fields whose 32-bit slot is
    saturated, so they must be consumed conditionally rather than by position.
    """
    cursor = 0
    while cursor + 4 <= len(extra):
        tag, size = struct.unpack("<HH", extra[cursor:cursor + 4])
        if tag == _ZIP64_EXTRA_TAG:
            payload = extra[cursor + 4:cursor + 4 + size]
            values = list(struct.unpack("<" + "Q" * (len(payload) // 8),
                                        payload[:len(payload) // 8 * 8]))
            if uncompressed == _SATURATED and values:
                uncompressed = values.pop(0)
            if compressed == _SATURATED and values:
                compressed = values.pop(0)
            if header_offset == _SATURATED and values:
                header_offset = values.pop(0)
            break
        cursor += 4 + size
    return uncompressed, compressed, header_offset


class RemoteZip:
    """Random access to a ZIP archive served over HTTP."""

    def __init__(self, url: str, timeout: int = 300):
        self.url = url
        self.timeout = timeout
        self._total = self._content_length()
        self._entries: dict[str, ZipEntry] | None = None

    def _content_length(self) -> int:
        request = urllib.request.Request(self.url, method="HEAD")
        with urllib.request.urlopen(request, timeout=self.timeout) as response:
            return int(response.headers["Content-Length"])

    def _range(self, start: int, length: int) -> bytes:
        if length <= 0:
            return b""
        end = min(start + length - 1, self._total - 1)
        request = urllib.request.Request(
            self.url, headers={"Range": f"bytes={start}-{end}"}
        )
        with urllib.request.urlopen(request, timeout=self.timeout) as response:
            # A server that ignores Range answers 200 with the whole body. Slicing
            # that as though it were the requested window would silently yield the
            # wrong bytes, and for a 59 GB archive it would also download
            # everything. Fail loudly instead.
            if response.status != 206:
                raise OSError(
                    f"server ignored Range request (HTTP {response.status}); "
                    f"{self.url} cannot be read incrementally"
                )
            return response.read()

    def _central_directory_location(self) -> tuple[int, int]:
        """Return (offset, size) of the central directory."""
        tail_length = min(65536 + 22, self._total)
        tail = self._range(self._total - tail_length, tail_length)

        locator = tail.rfind(_EOCD64_LOCATOR)
        if locator != -1:
            eocd64_offset = struct.unpack("<Q", tail[locator + 8:locator + 16])[0]
            eocd64 = self._range(eocd64_offset, 56)
            if eocd64[:4] != _EOCD64:
                raise ValueError("Zip64 end-of-central-directory record not found")
            size = struct.unpack("<Q", eocd64[40:48])[0]
            offset = struct.unpack("<Q", eocd64[48:56])[0]
            return offset, size

        eocd = tail.rfind(_EOCD)
        if eocd == -1:
            raise ValueError("not a ZIP archive: no end-of-central-directory record")
        size = struct.unpack("<I", tail[eocd + 12:eocd + 16])[0]
        offset = struct.unpack("<I", tail[eocd + 16:eocd + 20])[0]
        return offset, size

    @property
    def entries(self) -> dict[str, ZipEntry]:
        """Member name -> entry, parsed from the central directory on first use."""
        if self._entries is not None:
            return self._entries

        offset, size = self._central_directory_location()
        buffer = io.BytesIO(self._range(offset, size))
        found: dict[str, ZipEntry] = {}
        while buffer.read(4) == _CENTRAL_HEADER:
            header = buffer.read(42)
            compress_type = struct.unpack("<H", header[6:8])[0]
            compressed, uncompressed = struct.unpack("<II", header[16:24])
            name_len, extra_len, comment_len = struct.unpack("<HHH", header[24:30])
            header_offset = struct.unpack("<I", header[38:42])[0]
            name = buffer.read(name_len).decode("utf-8", "replace")
            extra = buffer.read(extra_len)
            buffer.read(comment_len)

            if _SATURATED in (uncompressed, compressed, header_offset):
                uncompressed, compressed, header_offset = _parse_zip64_extra(
                    extra, uncompressed, compressed, header_offset
                )
            found[name] = ZipEntry(name, header_offset, compressed,
                                   uncompressed, compress_type)

        self._entries = found
        return found

    def read(self, name: str) -> bytes:
        """Return the decompressed bytes of one member."""
        try:
            entry = self.entries[name]
        except KeyError:
            raise KeyError(
                f"{name!r} not in archive; {len(self.entries)} members available"
            ) from None

        # The central directory's name/extra lengths may differ from the local
        # header's, so the local header must be read to locate the data.
        local = self._range(entry.header_offset, 30)
        name_len, extra_len = struct.unpack("<HH", local[26:30])
        data_start = entry.header_offset + 30 + name_len + extra_len
        payload = self._range(data_start, entry.compressed_size)

        if entry.compress_type == 0:
            return payload
        return zlib.decompressobj(-zlib.MAX_WBITS).decompress(payload)

    def extract(self, name: str, dest: Path) -> Path:
        """Write one member to `dest`, creating parent directories."""
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(self.read(name))
        return dest
