"""ORCA lexical conversion and comments.

Basis: ORCA-INPUT. See
``docs/en/reference/formats/orcainp.md`` for verified
official manuals and tutorials. For field mappings and MolOP conventions, see
``docs/en/reference/data_model.md``.
"""

from __future__ import annotations


def parse_orca_float(token: str) -> float | None:
    cleaned = token.strip().strip('"').strip("'").rstrip(",")
    if "," in cleaned:
        cleaned = cleaned.split(",", 1)[0]
    if not cleaned:
        return None
    try:
        return float(cleaned.replace("D", "E").replace("d", "e"))
    except ValueError:
        return None


def parse_orca_int(token: str) -> int | None:
    cleaned = token.strip().strip('"').strip("'").rstrip(",")
    if not cleaned:
        return None
    try:
        return int(cleaned)
    except ValueError:
        try:
            as_float = float(cleaned)
        except ValueError:
            return None
        if not as_float.is_integer():
            return None
        return int(as_float)


def strip_orca_inline_comment(line: str) -> str:
    """ORCA-INPUT: # comments end at a second # or the end of the line.

    Quoted filenames retain literal # characters. A removed comment separates
    tokens, so text on its two sides cannot accidentally become one keyword.
    """
    result: list[str] = []
    quote: str | None = None
    in_comment = False
    for character in line:
        if in_comment:
            if character == "#":
                in_comment = False
            elif character in "\r\n":
                in_comment = False
                result.append(character)
        elif quote is not None:
            result.append(character)
            if character == quote:
                quote = None
        elif character == "#":
            in_comment = True
            result.append(" ")
        else:
            result.append(character)
            if character in {"'", '"'}:
                quote = character
    return "".join(result)


def split_orca_key_value(line: str) -> tuple[str, str]:
    content = strip_orca_inline_comment(line).strip()
    if "=" in content:
        key, value = content.split("=", 1)
        return key.strip(), value.strip()
    tokens = content.split(maxsplit=1)
    if not tokens:
        return "", ""
    if len(tokens) == 1:
        return tokens[0], ""
    return tokens[0], tokens[1]
