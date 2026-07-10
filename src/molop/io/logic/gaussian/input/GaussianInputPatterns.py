from __future__ import annotations

from molop.io.base_models.SearchPattern import MolOPPattern


class G16InputPatterns:
    """Gaussian input file grammar fragments."""

    OPTIONS = MolOPPattern(
        content_pattern=r"^\s*%(?P<key>[a-zA-Z0-9]+)(?:=(?P<value>[^\s]+))?",
        content_repeat=0,
        end_pattern=r"^\s*\n",
    )
    ROUTE = MolOPPattern(
        content_pattern=r"^\s*(?P<route>#[a-zA-Z0-9\(\)=,\/\\\s\n-]+)\n\n",
        end_pattern=r"^\s*\n",
    )
    TITLE = MolOPPattern(end_pattern=r"^\s*\n", content_pattern=r"^(?P<title>[^\n]+)")
    CHARGE_MULTIPLICITY = MolOPPattern(
        content_pattern=r"^\s*(?P<charge>-?\d+)\s+(?P<multiplicity>\d+)\s*\n",
        description="The charge and multiplicity of the Gaussian calculation. link 1",
    )
    CHARGE_MULTIPLICITY_LINE = MolOPPattern(
        content_pattern=r"^[+-]?\d+(?:[\s,/]+[+-]?\d+)+$",
    )
    INTEGER_TOKEN = MolOPPattern(
        content_pattern=r"^[+-]?\d+$",
    )
    CARTESIAN_COORD_START = MolOPPattern(
        content_pattern=r"^[-+]?\d+(?:\.\d+)?$",
    )
    ATOM_SPECIFICATION = MolOPPattern(
        content_pattern=r"^(?P<element_label>[A-Z][A-Za-z0-9]*)"
        r"(?:-(?P<atom_type>[A-Za-z][A-Za-z0-9]*)"
        r"(?:-(?P<charge>-?\d+(?:\.\d+)?))?)?"
        r"(?:\((?P<params>[A-Za-z][A-Za-z0-9]*=\w+"
        r"(?:,[A-Za-z][A-Za-z0-9]*=\w+)*)\))?$",
    )
    ELEMENT_LABEL_SYMBOL = MolOPPattern(
        content_pattern=r"(?P<symbol>[A-Z][a-z]?)\d+",
    )
    ZMAT_VARIABLE_ASSIGNMENT = MolOPPattern(
        content_pattern=r"^(?P<name>[A-Za-z][A-Za-z0-9_]*)\s+"
        r"(?P<value>[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)$",
    )
    ATOMS = MolOPPattern(
        content_pattern=r"^\s*(?P<symbol>[A-Z][a-z]*)"
        r"(?P<x>\s*-?\d+(?:\.\d*)?)"
        r"(?P<y>\s*-?\d+(?:\.\d*)?)"
        r"(?P<z>\s*-?\d+(?:\.\d*)?)",
        content_repeat=0,
    )
    LINK1_MARKER = MolOPPattern(
        content_pattern=r"(?i)--link1--",
        content_repeat=0,
    )
    GIC_FUNCTION_EXPRESSION = MolOPPattern(
        content_pattern=r"^(?P<name>[A-Za-z][A-Za-z0-9_]*)\((?P<args>.*)\)$",
    )
    GIC_STANDALONE_ATOM = MolOPPattern(
        content_pattern=r"^atom\s+\S+",
    )
    GIC_LEFT_ASSIGNMENT = MolOPPattern(
        content_pattern=r"^(?P<label>[A-Za-z][A-Za-z0-9]*)"
        r"(?:\s*[\(\[\{](?P<opts>[^\)\]\}]*)[\)\]\}])?$",
    )
    MODREDUNDANT_ATOM_REF = MolOPPattern(
        content_pattern=r"^(?P<atom_ref>\*|-?\d+)$",
    )
    MODREDUNDANT_COORDINATE_TYPE = MolOPPattern(
        content_pattern=r"^[A-Za-z]$",
    )
    ADDITIONAL_SECTION_SEPARATOR = MolOPPattern(
        content_pattern=r"\n\s*\n",
        content_repeat=0,
    )


g16_input_patterns = G16InputPatterns()
