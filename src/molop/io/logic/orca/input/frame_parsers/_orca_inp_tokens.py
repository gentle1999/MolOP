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
    stripped = line.lstrip()
    if stripped.startswith("#"):
        return ""
    if "#" in line:
        return line.split("#", 1)[0]
    return line


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
