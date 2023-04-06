from pathlib import Path


def is_gzipped(path: str | Path) -> bool:
    with open(path, "rb") as f:
        magic = f.read(2)
    return magic[0] == 0x1F and magic[1] == 0x8B
