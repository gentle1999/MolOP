from __future__ import annotations

from molop.io.base_models.SearchPattern import MolOPPattern


class ORCAInpPatterns:
    """ORCA input grammar fragments shared by file and frame parsers."""

    OUTPUT_PRINT_SETTING = MolOPPattern(
        content_pattern=r"(?i)^Print\s*\[\s*(?P<target>[^\]]+?)\s*\]\s*"
        r"(?:=|\s+)\s*(?P<value>.+)$",
    )
    COORD_HEADER = MolOPPattern(
        content_pattern=r"(?i)^\s*\*\s*"
        r"(?:xyz|cart|cartesian|int|internal|gzmt|xyzfile|gzmtfile|pdbfile)\b",
        content_repeat=0,
    )
    ATOM_TOKEN = MolOPPattern(
        content_pattern=r"^(?P<symbol>[A-Za-z]+)(?:\((?P<fragment>-?\d+)\))?$"
    )
    VARIABLE_OFFSET = MolOPPattern(
        content_pattern=r"^(?P<name>[A-Za-z_][A-Za-z0-9_]*)\s*"
        r"(?P<op>[+-])\s*(?P<offset>[+-]?\d+(?:\.\d*)?)$"
    )
    SOLVATION_TOKEN = MolOPPattern(
        content_pattern=r"^(?P<model>[A-Za-z][A-Za-z0-9-]*)\((?P<solvent>[^()]+)\)$"
    )
    MRCI_INLINE_EXCITATIONS = MolOPPattern(
        content_pattern=r"(?i)\bexcitations\s+(?P<excitations>\S+)"
    )
    MRCI_INLINE_REFS = MolOPPattern(content_pattern=r"(?i)\brefs\s+(?P<refs>.+?)\s+end\b")
    TRAILING_END = MolOPPattern(
        content_pattern=r"(?i)\s+end\s*$",
    )
    END_TOKEN_AT_LINE_END = MolOPPattern(
        content_pattern=r"(?i)\bend\s*$",
    )
    CAS_REFS = MolOPPattern(
        content_pattern=r"(?i)\bcas\s*\(\s*(?P<electrons>\d+)\s*,\s*"
        r"(?P<orbitals>\d+)\s*\)",
    )


orca_inp_patterns = ORCAInpPatterns()
