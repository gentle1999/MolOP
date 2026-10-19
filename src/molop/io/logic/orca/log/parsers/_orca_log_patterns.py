"""ORCA output grammar.

Basis: ORCA-OUTPUT and the output fixture manifest. See
``docs/en/reference/formats/orcaout.md`` for verified
official manuals and tutorials. For field mappings and MolOP conventions, see
``docs/en/reference/data_model.md``.
"""

from __future__ import annotations

from molop.io.base_models.SearchPattern import MolOPPattern


ORCA_FLOAT_PATTERN = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"


class ORCALogPatterns:
    """ORCA output grammar fragments used by file and frame parsers."""

    SECTION_HEADER = MolOPPattern(
        content_pattern=r"^[ \t]*(?P<title>[A-Z][A-Z0-9 ()*/=.,:+-]*)[ \t]*\r?\n"
        r"[ \t]*-{3,}[ \t]*$",
        content_repeat=0,
    )

    # ORCA-OUTPUT / ORCA-METHOD: ORBITAL ENERGIES, including optional Irrep.
    HF_TYPE = MolOPPattern(
        content_pattern=r"(?i)Hartree-Fock type\s+HFTyp\s+\.{3,}\s+(?P<kind>\S+)",
        content_repeat=0,
    )
    ROHF_CASE = MolOPPattern(content_pattern=r"(?i)\bROHF_case\s+(?P<case>\w+)", content_repeat=0)
    REFERENCE_MULTIPLICITY = MolOPPattern(
        content_pattern=r"Multiplicity\s+Mult\s+\.{3,}\s+(?P<value>\d+)",
        content_repeat=0,
    )
    ORBITAL_ROW = MolOPPattern(
        content_pattern=rf"^[ \t]*(?P<index>\d+)[ \t]+(?P<occupancy>{ORCA_FLOAT_PATTERN})"
        rf"[ \t]+(?P<energy>{ORCA_FLOAT_PATTERN})[ \t]+{ORCA_FLOAT_PATTERN}"
        r"(?:[ \t]+(?P<symmetry>\S+))?[ \t]*$",
        content_repeat=0,
    )
    # ORCA-POP: thresholded B(atom, atom) rows may contain multiple pairs.
    MAYER_THRESHOLD = MolOPPattern(
        content_pattern=rf"Mayer bond orders larger than\s+(?P<value>{ORCA_FLOAT_PATTERN})",
        content_repeat=0,
    )
    MAYER_BOND = MolOPPattern(
        content_pattern=rf"B\(\s*(?P<first>\d+)-[A-Za-z]+\s*,\s*"
        rf"(?P<second>\d+)-[A-Za-z]+\s*\)\s*:\s*(?P<value>{ORCA_FLOAT_PATTERN})",
        content_repeat=0,
    )
    # ORCA-SPIN: observed expectation, not the adjacent ideal S*(S+1) value.
    SPIN_SQUARE = MolOPPattern(
        content_pattern=rf"Expectation value of <S\*\*2>\s*:\s*"
        rf"(?P<value>{ORCA_FLOAT_PATTERN})",
        content_repeat=0,
    )
    # ORCA-THERMO: printed energy terms (including T*S) retain their labels.
    FREQUENCY_SCALE = MolOPPattern(
        content_pattern=rf"Scaling factor for frequencies\s*=\s*(?P<value>{ORCA_FLOAT_PATTERN})",
        content_repeat=0,
    )
    THERMAL_CUTOFF = MolOPPattern(
        content_pattern=rf"Cut-Off Frequency\s*\.{{3}}\s*(?P<value>{ORCA_FLOAT_PATTERN})\s+cm\^-1",
        content_repeat=0,
    )
    THERMAL_ROW = MolOPPattern(
        content_pattern=rf"^[ \t]*(?P<label>[A-Za-z][A-Za-z ()-]*?)"
        rf"[ \t]+(?:\.{{3}}[ \t]*)?(?P<value>{ORCA_FLOAT_PATTERN})[ \t]+"
        r"(?P<unit>Eh|K|atm|AMU)\b",
        content_repeat=0,
    )
    GIBBS_CORRECTION = MolOPPattern(
        content_pattern=rf"G-E\(el\)\s*\.{{3}}\s*(?P<value>{ORCA_FLOAT_PATTERN})\s+Eh",
        content_repeat=0,
    )
    ROTATIONAL_CONSTANTS = MolOPPattern(
        content_pattern=rf"Rotational constants in cm-1:\s*(?P<x>{ORCA_FLOAT_PATTERN})"
        rf"\s+(?P<y>{ORCA_FLOAT_PATTERN})\s+(?P<z>{ORCA_FLOAT_PATTERN})",
        content_repeat=0,
    )
    SYMMETRY_NUMBER = MolOPPattern(
        content_pattern=r"Symmetry Number:\s*(?P<value>\d+)",
        content_repeat=0,
    )
    # ORCA-VIB: blocked normal modes and source-indexed IR/Raman rows.
    MODE_COLUMNS = MolOPPattern(
        content_pattern=r"^[ \t]*(?P<indices>\d+(?:[ \t]+\d+)*)[ \t]*$",
        content_repeat=0,
    )
    MODE_ROW = MolOPPattern(
        content_pattern=rf"^[ \t]*(?P<index>\d+)[ \t]+"
        rf"(?P<values>{ORCA_FLOAT_PATTERN}(?:[ \t]+{ORCA_FLOAT_PATTERN})*)[ \t]*$",
        content_repeat=0,
    )
    SPECTRUM_ROW = MolOPPattern(
        content_pattern=rf"^[ \t]*(?P<index>\d+):[ \t]+(?P<frequency>{ORCA_FLOAT_PATTERN})"
        rf"[ \t]+(?P<values>{ORCA_FLOAT_PATTERN}(?:[ \t]+{ORCA_FLOAT_PATTERN})*)",
        content_repeat=0,
    )
    # ORCA-NMR: zero-based output indices; fixtures define exact row spelling.
    NMR_NUCLEUS = MolOPPattern(
        content_pattern=r"^[ \t]*Nucleus\s+(?P<index>\d+)\s*(?P<symbol>[A-Z][a-z]?)\s*:",
        content_repeat=0,
    )
    TENSOR_ROW = MolOPPattern(
        content_pattern=rf"^[ \t]*(?P<x>{ORCA_FLOAT_PATTERN})[ \t]+"
        rf"(?P<y>{ORCA_FLOAT_PATTERN})[ \t]+(?P<z>{ORCA_FLOAT_PATTERN})[ \t]*$",
        content_repeat=0,
    )
    SHIELDING_PRINCIPAL_VALUES = MolOPPattern(
        content_pattern=rf"^[ \t]*Total[ \t]+(?P<x>{ORCA_FLOAT_PATTERN})[ \t]+"
        rf"(?P<y>{ORCA_FLOAT_PATTERN})[ \t]+(?P<z>{ORCA_FLOAT_PATTERN})"
        rf"[ \t]+iso=\s*(?P<isotropic>{ORCA_FLOAT_PATTERN})",
        content_repeat=0,
    )
    SHIELDING_SUMMARY_ROW = MolOPPattern(
        content_pattern=rf"^[ \t]*(?P<index>\d+)[ \t]+[A-Z][a-z]?[ \t]+"
        rf"(?P<isotropic>{ORCA_FLOAT_PATTERN})[ \t]+(?P<anisotropy>{ORCA_FLOAT_PATTERN})[ \t]*$",
        content_repeat=0,
    )
    NMR_GAUGE = MolOPPattern(
        content_pattern=r"Choice of magnetic origin\s*\.{3}\s*(?P<gauge>\S+)",
        content_repeat=0,
    )
    COUPLING_ISOTROPIC = MolOPPattern(
        content_pattern=rf"J\[\s*(?P<first>\d+)\s*,\s*(?P<second>\d+)\s*\]"
        rf"\((?P<component>DSO|PSO|FC|SD|Total)\)[^\r\n]*?"
        rf"iso=\s*(?P<value>{ORCA_FLOAT_PATTERN})",
        content_repeat=0,
    )
    COUPLING_COLUMNS = MolOPPattern(
        content_pattern=r"(?P<index>\d+)[ \t]+[A-Z][a-z]?",
        content_repeat=0,
    )
    COUPLING_ROW = MolOPPattern(
        content_pattern=rf"^[ \t]*(?P<index>\d+)[ \t]+[A-Z][a-z]?[ \t]+"
        rf"(?P<values>{ORCA_FLOAT_PATTERN}(?:[ \t]+{ORCA_FLOAT_PATTERN})*)[ \t]*$",
        content_repeat=0,
    )

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

    # ORCA-COORD: header units apply to coordinates; MASS in the A.U. table is amu.
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
        content_pattern=rf"^[ \t]*(?P<index>\d+)[ \t]+(?P<symbol>[A-Z][a-z]?)[ \t]+"
        rf"{ORCA_FLOAT_PATTERN}[ \t]+\d+[ \t]+(?P<mass>{ORCA_FLOAT_PATTERN})[ \t]+"
        rf"{ORCA_FLOAT_PATTERN}[ \t]+{ORCA_FLOAT_PATTERN}[ \t]+{ORCA_FLOAT_PATTERN}[ \t]*$",
        content_repeat=0,
    )
    # ORCA-ENERGY: totals, references and correlation corrections are distinct.
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
    # ORCA-ENERGY / ORCA-OPT: analytic, numerical and MP2 gradient variants.
    GRADIENT_HEADER = MolOPPattern(
        content_pattern=r"^-+\s*\nCARTESIAN GRADIENT(?: \(NUMERICAL\))?\s*\n-+\s*$",
        content_repeat=0,
    )
    GRADIENT_ROW = MolOPPattern(
        content_pattern=rf"^\s*(?P<index>\d+)\s+(?P<symbol>[A-Z][a-z]?)\s*:\s+"
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
        content_pattern=rf"^\s*(?P<index>\d+)\s*:\s*"
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
    # ORCA-OPT: preserve printed thresholds and source YES/NO decisions.
    OPTIMIZATION_CONVERGENCE_METRIC = MolOPPattern(
        content_pattern=(
            rf"^\s*(?P<label>Energy change|RMS gradient|MAX gradient|RMS step|MAX step)\s+"
            rf"(?P<value>{ORCA_FLOAT_PATTERN})\s+"
            rf"(?P<threshold>{ORCA_FLOAT_PATTERN})\s+"
            r"(?P<converged>YES|NO)\s*$"
        ),
        content_repeat=0,
    )
    # ORCA-POP: charge-only and open-shell charge/spin tables.
    MULLIKEN_CHARGE_ROW = MolOPPattern(
        content_pattern=(
            rf"^\s*(?P<index>\d+)\s+(?P<symbol>[A-Z][a-z]?)\s*:\s*(?P<value>{ORCA_FLOAT_PATTERN})"
            rf"(?:\s+(?P<spin>{ORCA_FLOAT_PATTERN}))?\s*$"
        ),
        content_repeat=0,
    )
    HIRSHFELD_ROW = MolOPPattern(
        content_pattern=rf"^\s*(?P<index>\d+)\s+(?P<symbol>[A-Z][a-z]?)\s+"
        rf"(?P<charge>{ORCA_FLOAT_PATTERN})\s+"
        rf"(?P<spin>{ORCA_FLOAT_PATTERN})\s*$",
        content_repeat=0,
    )
    # ORCA-ELECTRIC: vector in e*bohr; the adjacent magnitude has both au and Debye.
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
    # ORCA-EXCITED: observed STATE/ROOT/absorption layouts in corpus fixtures.
    STATE = MolOPPattern(
        content_pattern=rf"(?i)STATE\s+(?P<root>\d+):\s+E=\s*"
        rf"(?P<au>{ORCA_FLOAT_PATTERN})\s+au\s+"
        rf"(?P<ev>{ORCA_FLOAT_PATTERN})\s+eV(?P<tail>.*)$",
        content_repeat=0,
    )
    STATE_MULTIPLICITY = MolOPPattern(content_pattern=r"\bMult\s+(?P<multiplicity>\d+)")
    STATE_IRREP = MolOPPattern(content_pattern=r"\bSym:\s*(?P<irrep>\S+)")
    ABS_TRANSITION = MolOPPattern(
        content_pattern=rf"^[ \t]*0-(?P<initial_label>\S+)[ \t]+->[ \t]+(?P<root>\d+)-"
        rf"(?P<label>(?P<multiplicity>\d+)(?P<irrep>[^\s\d]\S*))[ \t]+"
        rf"(?P<ev>{ORCA_FLOAT_PATTERN})\s+(?P<cm>{ORCA_FLOAT_PATTERN})\s+"
        rf"(?P<nm>{ORCA_FLOAT_PATTERN})\s+(?P<fosc>{ORCA_FLOAT_PATTERN})"
        rf"(?P<dipole>(?:[ \t]+{ORCA_FLOAT_PATTERN}){{4}})?[ \t]*$",
        content_repeat=0,
    )
    ROOT_TRANSITION = MolOPPattern(
        content_pattern=rf"^\s*(?P<root>\d+)\s+(?P<ev>{ORCA_FLOAT_PATTERN})\s+"
        rf"(?P<cm>{ORCA_FLOAT_PATTERN})\s+(?P<nm>{ORCA_FLOAT_PATTERN})\s+"
        rf"(?P<fosc>{ORCA_FLOAT_PATTERN})"
        rf"(?P<dipole>(?:[ \t]+{ORCA_FLOAT_PATTERN}){{4}})?[ \t]*$",
        content_repeat=0,
    )
    LEGACY_ROOT_TRANSITION = MolOPPattern(
        content_pattern=rf"^[ \t]*(?P<root>\d+)[ \t]+(?P<cm>{ORCA_FLOAT_PATTERN})"
        rf"[ \t]+(?P<nm>{ORCA_FLOAT_PATTERN})[ \t]+(?P<fosc>{ORCA_FLOAT_PATTERN})"
        rf"(?P<dipole>(?:[ \t]+{ORCA_FLOAT_PATTERN}){{4}})?[ \t]*$",
        content_repeat=0,
    )
    # ORCA-SOLVENT: CPCM/SMD output labels, not inferred solvent defaults.
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
