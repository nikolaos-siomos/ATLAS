"""Open static ATLAS PNG plots without invoking other image decoders."""
from contextlib import contextmanager
import os
import struct

from PIL import Image, UnidentifiedImageError


def _check_static_png(stream):
    """Check chunk framing and reject APNG before Pillow's plugin factory.

    Older Pillow can allocate an APNG disposal canvas during Image.open(),
    before its pixel-limit check. Checking image.is_animated afterwards is too
    late. Skip chunk payloads without decoding or allocating their contents;
    Pillow remains responsible for actual PNG validation and decoding.
    """
    size = os.fstat(stream.fileno()).st_size
    if stream.read(8) != b"\x89PNG\r\n\x1a\n":
        raise UnidentifiedImageError("ATLAS plots must be PNG images.")
    first = True
    while True:
        header = stream.read(8)
        if len(header) != 8:
            raise UnidentifiedImageError("Truncated PNG chunk header.")
        length, kind = struct.unpack(">I4s", header)
        end = stream.tell() + length + 4  # Payload and CRC.
        if length > 0x7fffffff or end > size:
            raise UnidentifiedImageError("Invalid PNG chunk length.")
        if first and (kind != b"IHDR" or length != 13):
            raise UnidentifiedImageError("PNG must start with an IHDR chunk.")
        first = False
        if kind in (b"acTL", b"fcTL", b"fdAT"):
            raise UnidentifiedImageError("ATLAS plots must be static PNG images; animated PNG is unsupported.")
        if kind == b"IEND":
            if length != 0:
                raise UnidentifiedImageError("Invalid PNG end chunk.")
            return
        stream.seek(end)


@contextmanager
def open_png(path):
    """Yield a static PNG, closing the image and file even if processing fails."""
    # Inspect and decode the same open file, rather than reopening a pathname
    # that could refer to a different file between the two operations.
    with open(path, "rb") as stream:
        _check_static_png(stream)
        stream.seek(0)
        with Image.open(stream, formats=["PNG"]) as image:
            yield image
