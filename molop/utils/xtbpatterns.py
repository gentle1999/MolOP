'''
Author: TMJ
Date: 2024-03-23 21:07:37
LastEditors: TMJ
LastEditTime: 2024-03-25 23:14:49
Description: 请填写简介
'''
import re
from typing import Dict, Tuple

xtboutpatterns: Dict[str, re.Pattern] = {
    "method": re.compile(r"Hamiltonian\s+([GFN\d\-xTB]+)"),
    "solvent model": re.compile(r"Solvation model:\s+([a-zA-Z\d]+)"),
    "solvent": re.compile(r"Solvent\s+([a-zA-Z\d\-]+)\n"),
    "temperature": re.compile(r"electronic temp.\s+(\d+.\d+)\s+K"),
    "version": re.compile(
        r"xtb (version \d+\.\d+\.\d+\s+\([0-9a-z]+\) compiled by ['0-9a-zA-Z@_-]+ on \d{4}-\d{2}-\d{2})"
    ),
    "old version": re.compile(
        r"(Version\s+[0-9\.]+\s+[0-9a-zA-Z]+\s+[\(\)a-zA-Z0-9]+)"
    ),
    "charge": re.compile(r"charge\s+:\s+([-+0-9]+)"),
    "multiplicity": re.compile(r"(--uhf|-u) (\d+)"),
    "total charge": re.compile(r"total charge\s+([\-\+0-9]+)"),
    "parameter": re.compile(r"program call\s+:\s+([0-9a-zA-Z\\/_\-\s\.]+)\n"),
    "n atoms": re.compile(r"number of atoms\s+:\s+(\d+)"),
    "state": re.compile(r"GEOMETRY OPTIMIZATION CONVERGED"),
    "homo": re.compile(r"([\+\-0-9.]+)\s+\(HOMO\)"),
    "lumo": re.compile(r"([\+\-0-9.]+)\s+\(LUMO\)"),
    "wiberg_start": re.compile(r"#   Z sym  total"),
    "wiberg_end": re.compile(r"molecular dipole"),
    "index": re.compile(r" \d+\s+\d+\s+[a-zA-Z]+"),
    "bond": re.compile(r"(\d+)\s+[a-zA-Z]+"),
    "value": re.compile(r"[\s-]\d+\.\d+"),
    "dipole_start": re.compile(r"molecular dipole"),
    "dipole_end": re.compile(r"molecular quadrupole"),
    "coords_start": re.compile(r"final structure"),
    "coords_new": re.compile(
        r"([a-zA-Z]+)\s+([\s\-]\d+\.\d+)\s+([\s\-]\d+\.\d+)\s+([\s\-]\d+\.\d+)"
    ),
    "coords_old": re.compile(
        r"([\s\-]\d+\.\d+)\s+([\s\-]\d+\.\d+)\s+([\s\-]\d+\.\d+)\s+([a-zA-Z]+)"
    ),
    "coords_end": re.compile(r"Bond Distances"),
    "coords_type": re.compile(r"coordinate file\s*:\s*[^\s]+.(xyz|sdf|mol)"),
    "zp correct": re.compile(r"zero point energy\s+([\-\+0-9.]+)"),
    "H correct": re.compile(
        r"\d+.\d+\s+\d+.\d+E[+-]\d+\s+(\d+.\d+E[+-]\d+)\s+\d+.\d+E[+-]\d+\s+\d+.\d+E[+-]\d+"
    ),
    "G correct": re.compile(
        r"\d+.\d+\s+\d+.\d+E[+-]\d+\s+\d+.\d+E[+-]\d+\s+\d+.\d+E[+-]\d+\s+(\d+.\d+E[+-]\d+)"
    ),
    "sum H": re.compile(r"TOTAL ENTHALPY\s+([\-\+0-9.]+)"),
    "sum G": re.compile(r"TOTAL FREE ENERGY\s+([\-\+0-9.]+)"),
    "sum E": re.compile(r"TOTAL ENERGY\s+([\-\+0-9.]+)"),
}
