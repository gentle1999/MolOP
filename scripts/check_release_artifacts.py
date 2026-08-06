"""Validate MolOP wheel and source distribution release contracts."""

from __future__ import annotations

import argparse
import tarfile
import zipfile
from email.parser import BytesParser
from email.policy import default
from pathlib import Path, PurePosixPath


MAX_SDIST_BYTES = 5 * 1024 * 1024
ALLOWED_SDIST_ROOTS = {
    ".gitignore",
    "LICENSE",
    "PKG-INFO",
    "README.md",
    "README.zh.md",
    "pyproject.toml",
    "src",
}
REQUIRED_SDIST_FILES = {
    "LICENSE",
    "PKG-INFO",
    "README.md",
    "README.zh.md",
    "pyproject.toml",
    "src/molop/py.typed",
}


def _single_artifact(dist_dir: Path, pattern: str, label: str) -> Path:
    artifacts = sorted(dist_dir.glob(pattern))
    if len(artifacts) != 1:
        raise ValueError(f"Expected exactly one {label}, found {len(artifacts)}")
    return artifacts[0]


def _safe_archive_path(name: str) -> PurePosixPath:
    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"Unsafe archive member: {name}")
    return path


def _check_sdist(sdist: Path) -> str:
    size = sdist.stat().st_size
    if size > MAX_SDIST_BYTES:
        raise ValueError(f"sdist is {size} bytes; limit is {MAX_SDIST_BYTES}")

    with tarfile.open(sdist, "r:gz") as archive:
        files = {
            _safe_archive_path(member.name) for member in archive.getmembers() if member.isfile()
        }

    roots = {path.parts[0] for path in files}
    if len(roots) != 1:
        raise ValueError(f"sdist must have one archive root, found {sorted(roots)}")
    root = next(iter(roots))
    relative_files = {PurePosixPath(*path.parts[1:]) for path in files if len(path.parts) > 1}
    top_level = {path.parts[0] for path in relative_files}

    unexpected = top_level - ALLOWED_SDIST_ROOTS
    if unexpected:
        raise ValueError(f"Unexpected sdist paths: {sorted(unexpected)}")

    missing = {PurePosixPath(path) for path in REQUIRED_SDIST_FILES} - relative_files
    if missing:
        raise ValueError(f"Missing required sdist files: {sorted(map(str, missing))}")

    unexpected_sources = {
        str(path)
        for path in relative_files
        if path.parts[0] == "src" and (len(path.parts) < 2 or path.parts[1] != "molop")
    }
    if unexpected_sources:
        raise ValueError(f"Unexpected source trees: {sorted(unexpected_sources)}")

    return root.removeprefix("molop-")


def _check_wheel(wheel: Path) -> str:
    with zipfile.ZipFile(wheel) as archive:
        files = {_safe_archive_path(name) for name in archive.namelist() if not name.endswith("/")}
        dist_info_dirs = {path.parts[0] for path in files if path.parts[0].endswith(".dist-info")}
        if len(dist_info_dirs) != 1:
            raise ValueError(f"wheel must have one .dist-info directory, found {dist_info_dirs}")
        dist_info = next(iter(dist_info_dirs))

        required = {
            PurePosixPath("molop/py.typed"),
            PurePosixPath(f"{dist_info}/METADATA"),
            PurePosixPath(f"{dist_info}/WHEEL"),
            PurePosixPath(f"{dist_info}/entry_points.txt"),
        }
        missing = required - files
        if missing:
            raise ValueError(f"Missing required wheel files: {sorted(map(str, missing))}")

        top_level = {path.parts[0] for path in files}
        unexpected = top_level - {"molop", dist_info}
        if unexpected:
            raise ValueError(f"Unexpected wheel paths: {sorted(unexpected)}")

        metadata_bytes = archive.read(f"{dist_info}/METADATA")
        entry_points = archive.read(f"{dist_info}/entry_points.txt").decode("utf-8")

    metadata = BytesParser(policy=default).parsebytes(metadata_bytes)
    if metadata["Name"] != "molop":
        raise ValueError(f"Unexpected package name: {metadata['Name']}")
    if metadata["Requires-Python"] != ">=3.10":
        raise ValueError(f"Unexpected Requires-Python: {metadata['Requires-Python']}")
    if "Typing :: Typed" not in metadata.get_all("Classifier", []):
        raise ValueError("Wheel metadata does not declare Typing :: Typed")
    if "molop = molop.cli.app:app" not in entry_points:
        raise ValueError("Wheel does not expose the molop console entry point")

    version = metadata["Version"]
    if not version:
        raise ValueError("Wheel metadata does not contain a version")
    return version


def check_release_artifacts(dist_dir: Path) -> tuple[Path, Path]:
    if not dist_dir.is_dir():
        raise ValueError(f"Distribution directory does not exist: {dist_dir}")

    wheel = _single_artifact(dist_dir, "*.whl", "wheel")
    sdist = _single_artifact(dist_dir, "*.tar.gz", "source distribution")
    wheel_version = _check_wheel(wheel)
    sdist_version = _check_sdist(sdist)
    if wheel_version != sdist_version:
        raise ValueError(
            f"Wheel version {wheel_version!r} does not match sdist version {sdist_version!r}"
        )
    return wheel, sdist


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dist_dir", type=Path, help="Directory containing one wheel and one sdist")
    args = parser.parse_args()

    try:
        wheel, sdist = check_release_artifacts(args.dist_dir)
    except ValueError as exc:
        parser.error(str(exc))

    print(f"wheel: {wheel.name} ({wheel.stat().st_size} bytes)")
    print(f"sdist: {sdist.name} ({sdist.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
