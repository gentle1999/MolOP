from __future__ import annotations

from molop.io.base_models.SearchPattern import MolOPPattern


ORCA_FLOAT_PATTERN = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"


class ORCALogPatterns:
    """ORCA output grammar fragments used by file and frame parsers."""

    FLOAT = MolOPPattern(content_pattern=rf"(?P<value>{ORCA_FLOAT_PATTERN})", content_repeat=0)
    BANNER = MolOPPattern(content_pattern=r"(?i)\*\s+O\s+R\s+C\s+A\s+\*")
    VERSION = MolOPPattern(content_pattern=r"Program Version\s+(?P<version>[0-9][^\s]*)")
    INPUT_BLOCK = MolOPPattern(
        content_pattern=r"(?s)=+\s*\n\s*INPUT FILE\s*\n=+\s*\n(?P<body>.*?)\n=+",
        content_repeat=0,
    )
    INPUT_NAME = MolOPPattern(content_pattern=r"^\s*NAME\s*=\s*(?P<name>.+?)\s*$")
    INPUT_LINE = MolOPPattern(content_pattern=r"^\|\s*\d+>\s?(?P<line>.*)$", content_repeat=0)

    COORD_HEADER = MolOPPattern(
        content_pattern=r"^-+\s*\nCARTESIAN COORDINATES \(ANGSTROEM\)\s*\n-+\s*$",
        content_repeat=0,
    )
    COORD_ROW = MolOPPattern(
        content_pattern=rf"^\s*(?P<symbol>[A-Z][a-z]?)\s+"
        rf"(?P<x>{ORCA_FLOAT_PATTERN})\s+"
        rf"(?P<y>{ORCA_FLOAT_PATTERN})\s+"
        rf"(?P<z>{ORCA_FLOAT_PATTERN})\s*$",
        content_repeat=0,
    )
    FINAL_ENERGY = MolOPPattern(
        content_pattern=rf"FINAL SINGLE POINT ENERGY\s+(?P<energy>{ORCA_FLOAT_PATTERN})",
        content_repeat=0,
    )
    SCF_ENERGY = MolOPPattern(
        content_pattern=rf"Total Energy\s*:\s*(?P<energy>{ORCA_FLOAT_PATTERN})\s+Eh",
        content_repeat=0,
    )
    MP2_ENERGY = MolOPPattern(
        content_pattern=rf"(?i)(?:MP2 TOTAL ENERGY:\s*(?P<mp2_total>{ORCA_FLOAT_PATTERN})\s*Eh|"
        rf"E\(MP2\)\s*=\s*(?P<mp2_corr>{ORCA_FLOAT_PATTERN}))",
        content_repeat=0,
    )
    MP3_ENERGY = MolOPPattern(
        content_pattern=rf"(?i)E\(MP3\)\s*=\s*(?P<energy>{ORCA_FLOAT_PATTERN})",
        content_repeat=0,
    )
    CCSD_ENERGY = MolOPPattern(
        content_pattern=rf"(?i)E\(CCSD(?:\(T\))?\)\s*=\s*(?P<energy>{ORCA_FLOAT_PATTERN})",
        content_repeat=0,
    )
    GRADIENT_HEADER = MolOPPattern(
        content_pattern=r"^-+\s*\nCARTESIAN GRADIENT(?: \(NUMERICAL\))?\s*\n-+\s*$",
        content_repeat=0,
    )
    GRADIENT_ROW = MolOPPattern(
        content_pattern=rf"^\s*\d+\s+[A-Z][a-z]?\s*:\s+"
        rf"(?P<x>{ORCA_FLOAT_PATTERN})\s+"
        rf"(?P<y>{ORCA_FLOAT_PATTERN})\s+"
        rf"(?P<z>{ORCA_FLOAT_PATTERN})\s*$",
        content_repeat=0,
    )
    MULLIKEN_CHARGE_ROW = MolOPPattern(
        content_pattern=rf"^\s*\d+\s+[A-Z][a-z]?\s*:\s*(?P<value>{ORCA_FLOAT_PATTERN})\s*$",
        content_repeat=0,
    )
    HIRSHFELD_ROW = MolOPPattern(
        content_pattern=rf"^\s*\d+\s+[A-Z][a-z]?\s+"
        rf"(?P<charge>{ORCA_FLOAT_PATTERN})\s+"
        rf"(?P<spin>{ORCA_FLOAT_PATTERN})\s*$",
        content_repeat=0,
    )
    DIPOLE = MolOPPattern(
        content_pattern=rf"Total Dipole Moment\s*:\s*(?P<x>{ORCA_FLOAT_PATTERN})\s+"
        rf"(?P<y>{ORCA_FLOAT_PATTERN})\s+(?P<z>{ORCA_FLOAT_PATTERN})",
        content_repeat=0,
    )
    POLAR_ISOTROPIC = MolOPPattern(
        content_pattern=rf"(?i)Isotropic polarizability\s*:\s*(?P<value>{ORCA_FLOAT_PATTERN})",
        content_repeat=0,
    )
    RUN_TIME = MolOPPattern(
        content_pattern=r"(?i)TOTAL RUN TIME:\s*"
        r"(?P<days>\d+)\s+days\s+"
        r"(?P<hours>\d+)\s+hours\s+"
        r"(?P<minutes>\d+)\s+minutes?\s+"
        r"(?P<seconds>\d+)\s+seconds?\s+"
        r"(?P<msec>\d+)\s+msec",
        content_repeat=0,
    )
    FREQUENCY = MolOPPattern(
        content_pattern=rf"^\s*(?P<idx>\d+):\s*(?P<frequency>{ORCA_FLOAT_PATTERN})\s+"
        r"cm\*\*-1(?:\s+(?P<label>\S+))?\s*$",
        content_repeat=0,
    )
    STATE = MolOPPattern(
        content_pattern=rf"(?i)STATE\s+(?P<root>\d+):\s+E=\s*"
        rf"(?P<au>{ORCA_FLOAT_PATTERN})\s+au\s+"
        rf"(?P<ev>{ORCA_FLOAT_PATTERN})\s+eV(?P<tail>.*)$",
        content_repeat=0,
    )
    STATE_MULTIPLICITY = MolOPPattern(content_pattern=r"\bMult\s+(?P<multiplicity>\d+)")
    STATE_IRREP = MolOPPattern(content_pattern=r"\bSym:\s*(?P<irrep>\S+)")
    ABS_TRANSITION = MolOPPattern(
        content_pattern=rf"^\s*\S+\s+->\s+(?P<root>\d+)-(?P<label>\S+)\s+"
        rf"(?P<ev>{ORCA_FLOAT_PATTERN})\s+(?P<cm>{ORCA_FLOAT_PATTERN})\s+"
        rf"(?P<nm>{ORCA_FLOAT_PATTERN})\s+(?P<fosc>{ORCA_FLOAT_PATTERN})",
        content_repeat=0,
    )
    ROOT_TRANSITION = MolOPPattern(
        content_pattern=rf"^\s*(?P<root>\d+)\s+(?P<ev>{ORCA_FLOAT_PATTERN})\s+"
        rf"(?P<cm>{ORCA_FLOAT_PATTERN})\s+(?P<nm>{ORCA_FLOAT_PATTERN})\s+"
        rf"(?P<fosc>{ORCA_FLOAT_PATTERN})",
        content_repeat=0,
    )
    SOLVENT_NAME = MolOPPattern(content_pattern=r"Solvent:\s*\.\.\.\s*(?P<solvent>\S+)")
    SOLVENT_EPSILON = MolOPPattern(
        content_pattern=rf"Epsilon\s*\.\.\.\s*(?P<epsilon>{ORCA_FLOAT_PATTERN})"
    )

    @staticmethod
    def optimization_metric(label: str) -> MolOPPattern:
        return MolOPPattern(
            content_pattern=rf"(?i){MolOPPattern.escape_literal(label)}\s+(?P<value>{ORCA_FLOAT_PATTERN})",
            content_repeat=0,
        )


orca_log_patterns = ORCALogPatterns()
