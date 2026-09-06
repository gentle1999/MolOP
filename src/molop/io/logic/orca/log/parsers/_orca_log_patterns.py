from __future__ import annotations

from molop.io.base_models.SearchPattern import MolOPPattern


ORCA_FLOAT_PATTERN = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"


class ORCALogPatterns:
    """ORCA output grammar fragments used by file and frame parsers."""

    FLOAT = MolOPPattern(content_pattern=rf"(?P<value>{ORCA_FLOAT_PATTERN})", content_repeat=0)
    BANNER = MolOPPattern(content_pattern=r"(?i)\*\s+O\s+R\s+C\s+A\s+\*")
    JOB_NUMBER_LINE = MolOPPattern(
        content_pattern=(
            r"^[^\S\r\n]*\*+[^\S\r\n]*JOB\s+NUMBER\s+(?P<number>\d+)"
            r"[^\S\r\n]*\*+[^\S\r\n]*\r?$"
        ),
        content_repeat=0,
    )
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
    COORD_AU_HEADER = MolOPPattern(
        content_pattern=r"^[ \t]*-+[ \t]*\r?\nCARTESIAN COORDINATES \(A\.U\.\)[ \t]*\r?\n"
        r"[ \t]*-+[ \t]*$",
        content_repeat=0,
    )
    COORD_ROW = MolOPPattern(
        content_pattern=rf"^\s*(?P<symbol>[A-Z][a-z]?)\s+"
        rf"(?P<x>{ORCA_FLOAT_PATTERN})\s+"
        rf"(?P<y>{ORCA_FLOAT_PATTERN})\s+"
        rf"(?P<z>{ORCA_FLOAT_PATTERN})\s*$",
        content_repeat=0,
    )
    COORD_AU_MASS_ROW = MolOPPattern(
        content_pattern=rf"^[ \t]*(?P<index>\d+)[ \t]+[A-Z][a-z]?[ \t]+"
        rf"{ORCA_FLOAT_PATTERN}[ \t]+\d+[ \t]+(?P<mass>{ORCA_FLOAT_PATTERN})[ \t]+"
        rf"{ORCA_FLOAT_PATTERN}[ \t]+{ORCA_FLOAT_PATTERN}[ \t]+{ORCA_FLOAT_PATTERN}[ \t]*$",
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
        rf"E\(MP2\)\s*=\s*(?P<mp2_equation_total>{ORCA_FLOAT_PATTERN}))",
        content_repeat=0,
    )
    MP2_CORRELATION_ENERGY = MolOPPattern(
        content_pattern=(
            rf"(?i)(?:MP2 CORRELATION ENERGY\s*:\s*"
            rf"(?P<labeled>{ORCA_FLOAT_PATTERN})\s*Eh|"
            rf"EC\(MP2\)\s*=\s*(?P<equation>{ORCA_FLOAT_PATTERN})|"
            rf"E\(MP2\)\s*\.{{3}}\s*(?P<component>{ORCA_FLOAT_PATTERN}))"
        ),
        content_repeat=0,
    )
    MP3_ENERGY = MolOPPattern(
        content_pattern=rf"(?i)E\(MP3\)\s*=\s*(?P<energy>{ORCA_FLOAT_PATTERN})",
        content_repeat=0,
    )
    MP3_CORRELATION_ENERGY = MolOPPattern(
        content_pattern=rf"(?i)EC\(MP3\)\s*=\s*(?P<energy>{ORCA_FLOAT_PATTERN})",
        content_repeat=0,
    )
    MP3_COMPONENT_ENERGY = MolOPPattern(
        content_pattern=rf"(?i)\bE3\s*=\s*(?P<energy>{ORCA_FLOAT_PATTERN})",
        content_repeat=0,
    )
    CCSD_ENERGY = MolOPPattern(
        content_pattern=rf"(?i)E\(CCSD\)\s*(?:=|\.{{3}})\s*(?P<energy>{ORCA_FLOAT_PATTERN})",
        content_repeat=0,
    )
    CCSD_TOTAL_ENERGY = MolOPPattern(
        content_pattern=rf"(?i)^E\(TOT\)\s*(?:=|\.{{3}})\s*"
        rf"(?P<energy>{ORCA_FLOAT_PATTERN})\s*$",
        content_repeat=0,
    )
    CCSD_CORRELATION_ENERGY = MolOPPattern(
        content_pattern=rf"(?i)^E\(CORR\)\s*(?:=|\.{{3}})\s*"
        rf"(?P<energy>{ORCA_FLOAT_PATTERN})\s*$",
        content_repeat=0,
    )
    CCSD_T_ENERGY = MolOPPattern(
        content_pattern=rf"(?i)E\(CCSD\(T\)\)\s*(?:=|\.{{3}})\s*(?P<energy>{ORCA_FLOAT_PATTERN})",
        content_repeat=0,
    )
    TRIPLES_CORRECTION = MolOPPattern(
        content_pattern=rf"(?i)^\s*triples correction \(T\)\s*\.{{3}}\s*"
        rf"(?P<energy>{ORCA_FLOAT_PATTERN})\s*$",
        content_repeat=0,
    )
    SCALED_TRIPLES_CORRECTION = MolOPPattern(
        content_pattern=rf"(?i)^\s*scaled triples correction \(T\)\s*\.{{3}}\s*"
        rf"(?P<energy>{ORCA_FLOAT_PATTERN})\s*$",
        content_repeat=0,
    )
    FINAL_CORRELATION_ENERGY = MolOPPattern(
        content_pattern=rf"(?i)Final correlation energy\s*\.{{3}}\s*"
        rf"(?P<energy>{ORCA_FLOAT_PATTERN})",
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
    MP2_GRADIENT_HEADER = MolOPPattern(
        content_pattern=r"^\s*The final MP2 gradient\s*$",
        content_repeat=0,
    )
    MP2_GRADIENT_ROW = MolOPPattern(
        content_pattern=rf"^\s*\d+\s*:\s*"
        rf"(?P<x>{ORCA_FLOAT_PATTERN})\s+"
        rf"(?P<y>{ORCA_FLOAT_PATTERN})\s+"
        rf"(?P<z>{ORCA_FLOAT_PATTERN})\s*$",
        content_repeat=0,
    )
    SCF_CONVERGED = MolOPPattern(
        content_pattern=r"(?i)\bSCF\s+CONVERGED(?:\s+AFTER\s+\d+\s+CYCLES?)?",
        content_repeat=0,
    )
    SCF_FAILED = MolOPPattern(
        content_pattern=(
            r"(?i)(?:\bSCF\s+NOT\s+CONVERGED(?:\s+AFTER\s+\d+\s+CYCLES?)?|"
            r"\bSCF\s+DID\s+NOT\s+CONVERGE|"
            r"\bSCF\s+CONVERGENCE(?:\s+HAS)?\s+FAIL(?:ED|URE)|"
            r"\bSCF\s+HAS\s+NOT\s+CONVERGED)"
        ),
        content_repeat=0,
    )
    OPTIMIZATION_CONVERGENCE_METRIC = MolOPPattern(
        content_pattern=(
            rf"^\s*(?P<label>Energy change|RMS gradient|MAX gradient|RMS step|MAX step)\s+"
            rf"(?P<value>{ORCA_FLOAT_PATTERN})\s+"
            rf"(?P<threshold>{ORCA_FLOAT_PATTERN})\s+"
            r"(?P<converged>YES|NO)\s*$"
        ),
        content_repeat=0,
    )
    MULLIKEN_CHARGE_ROW = MolOPPattern(
        content_pattern=(
            rf"^\s*\d+\s+[A-Z][a-z]?\s*:\s*(?P<value>{ORCA_FLOAT_PATTERN})"
            rf"(?:\s+(?P<spin>{ORCA_FLOAT_PATTERN}))?\s*$"
        ),
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
