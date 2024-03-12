"""
Author: TMJ
Date: 2023-11-02 15:36:39
LastEditors: TMJ
LastEditTime: 2024-02-21 21:07:05
Description: 请填写简介
"""


import re
from typing import Dict, Tuple


def parameter_comment_parser(
    parameter_comment: str,
) -> Tuple[Dict, Dict, str, str, str]:
    """
    Parse parameter comment from G16 log file.

    Args:
        parameter_comment (str): The parameter comment from the G16 log file.

    Returns:
        - link0 (Dict): A dictionary containing parameter-value pairs extracted from the first line of the parameter comment.
        - route_params (Dict): A dictionary containing parameter-value pairs extracted from the route section of the parameter comment.
        - dieze_tag (str): The dieze tag extracted from the route section of the parameter comment.
        - functional (str): The functional extracted from the route section of the parameter comment.
        - basis_set (str): The basis set extracted from the route section of the parameter comment.
    """
    link0 = {
        f"{para.split('=')[0]}": f"{para.split('=')[1]}"
        for para in g16logpatterns["link0"].findall(parameter_comment.split("\n")[0])
    }
    route = parameter_comment.split("\n")[1]
    scrf_patt = g16logpatterns["scrf_patt"]
    multi_params_patt = g16logpatterns["multi_params_patt"]
    functional = basis_set = None
    route_params = {}
    dieze_tag = None

    if route:
        if "/" in route:
            tok = route.split("/")
            functional = tok[0].split()[-1]
            basis_set = tok[1].split()[0]
            for tok in [functional, basis_set, "/"]:
                route = route.replace(tok, "")

        for tok in route.split():
            if scrf_patt.match(tok):
                m = scrf_patt.match(tok)
                route_params[m.group(1)] = m.group(2)
            elif tok.upper() in ["#", "#N", "#P", "#T"]:
                # does not store # in route to avoid error in input
                dieze_tag = "#N" if tok == "#" else tok
                continue
            else:
                m = re.match(multi_params_patt, tok.strip("#"))
                if m:
                    pars = {}
                    for par in m.group(2).split(","):
                        p = par.split("=")
                        pars[p[0]] = None if len(p) == 1 else p[1]
                    route_params[m.group(1)] = pars
                else:
                    d = tok.strip("#").split("=")
                    route_params[d[0]] = None if len(d) == 1 else d[1]

    return link0, route_params, dieze_tag, functional, basis_set


g16logpatterns: Dict[str, re.Pattern] = {
    "n atoms": re.compile(r"NAtoms=\s*(\d+)"),
    "coords start": re.compile(r"Standard orientation:"),
    "coords": re.compile(
        r"\s+\d+\s+(\d+)\s+\d+\s+([\s\-]\d+\.\d+)\s+([\s\-]\d+\.\d+)\s+([\s\-]\d+\.\d+)"
    ),
    "link0": re.compile(r"\%[a-z]+=[^\s]+"),
    "scrf_patt": re.compile(r"^([sS][cC][rR][fF])\s*=\s*(.+)"),
    "multi_params_patt": re.compile(r"^([A-z]+[0-9]*)[\s=]+\((.*)\)$"),
    "mulliken start": re.compile(
        r"(Mulliken charges:|Mulliken atomic charges|Mulliken charges and spin densities:)"
    ),
    "mulliken match": re.compile(r"\d+\s+[A-Z][a-z]?\s+([\s-]\d+.\d+)"),
    "mulliken end": re.compile(r"(Sum of Mulliken )(.*)(charges)\s*=\s*(\D)"),
    "spin density start": re.compile(r"Mulliken charges and spin densities:"),
    "spin density match": re.compile(
        r"\d+\s+[A-Z][a-z]?\s+[\s-]\d+.\d+\s+([\s-]\d+.\d+)"
    ),
    "spin density end": re.compile(r"(Sum of Mulliken )(.*)(charges)\s*=\s*(\D)"),
    "gradients start": re.compile(r"Center\s+Atomic\s+Forces\s+\(Hartrees/Bohr\)"),
    "gradients match": re.compile(
        r"\s+\d+\s+\d+\s+([\s\-]\d+\.\d+)\s+([\s\-]\d+\.\d+)\s+([\s\-]\d+\.\d+)"
    ),
    "gradients end": re.compile(r"Cartesian\s+Forces:\s+Max.*RMS.*"),
    "opt stat": re.compile(
        r"(Maximum\s+Force|RMS\s+Force|Maximum\s+Displacement|RMS\s+Displacement)\s+[\d.]+\s+[\d.]+\s+(NO|YES)"
    ),
    "termination": re.compile(r"(Normal|Error) termination"),
    "failure reason": re.compile(r"(! Non-Optimized Parameters !|Convergence failure)"),
    "energy": re.compile(r"E\(.*\)\s*=\s*([-\.\d]+)\s+"),
    "energy alter": re.compile(r"HF=([\-0-9\s+\.]+)\\"),
    "spins": re.compile(r"<S\*\*2>=([\s\-\d.]+)\s+S=([\s\-\d.]+)"),
    "orbital": re.compile(r"(Alpha|Beta)\s*(occ.|virt.)\s*eigenvalues -- (.*)"),
    "freq start": re.compile(
        r"Harmonic\sfrequencies\s+\(cm\*\*-1\),\sIR\sintensities.*Raman.*"
    ),
    "freq": re.compile(
        r"Frequencies -- \s+([\s-]\d+.\d+)\s+([\s-]\d+.\d+)\s+([\s-]\d+.\d+)"
    ),
    "freq red. masses": re.compile(
        r"Red. masses -- \s+([\s-]\d+.\d+)\s+([\s-]\d+.\d+)\s+([\s-]\d+.\d+)"
    ),
    "freq frc consts": re.compile(
        r"Frc consts  -- \s+([\s-]\d+.\d+)\s+([\s-]\d+.\d+)\s+([\s-]\d+.\d+)"
    ),
    "freq IR Inten": re.compile(
        r"IR Inten    -- \s+([\s-]\d+.\d+)\s+([\s-]\d+.\d+)\s+([\s-]\d+.\d+)"
    ),
    "freq mode": re.compile(
        r"\s\s\s\s([\s-]\d.\d\d)\s\s([\s-]\d.\d\d)\s\s([\s-]\d.\d\d)"
    ),
    "sum energies": re.compile(
        r"Sum of electronic and (thermal Free Energies|thermal Enthalpies|thermal Energies|zero-point Energies)=\s+([\-0-9.]+)"
    ),
    "corrections": re.compile(r"(Zero-point|Thermal) correction(.*)=\s+([\d\.-]+)"),
    "wiberg_start": re.compile(r"Wiberg bond index matrix in the NAO basis"),
    "wiberg_end": re.compile(r"Wiberg bond index, Totals by atom"),
    "atom_atom_overlap_start": re.compile(r"Atom-atom overlap-weighted NAO bond order"),
    "atom_atom_overlap_end": re.compile(r"Atom-atom overlap-weighted NAO bond order, Totals by atom"),
    "mo_start": re.compile(r"MO bond order"),
    "mo_end": re.compile(r"MO atomic valencies"),
    "digit": re.compile(r"[\s-]\d+.\d+"),
    "dipole": re.compile(
        r"X=\s+([\s-]\d+\.\d+)\s+Y=\s+([\s-]\d+\.\d+)\s+Z=\s+([\s-]\d+\.\d+)\s+"
    ),
    "nbo charge start": re.compile(r"Summary of Natural Population Analysis:"),
    "nbo charge end": re.compile(r"=\n\s+\* Total \*"),
    "nbo charge match": re.compile(r"[A-Za-z]+\s+\d+\s+([-\d.]+)"),
    "nbo summary start": re.compile(r"Natural Bond Orbitals \(Summary\):"),
    "nbo summary end": re.compile(r"-\n\s+Total Lewis"),
    "nbo summary match": re.compile(r"BD\s+\(\s+\d\)\s+[A-Za-z]+\s+(\d+)\s+-[\sA-Za-z]+\s+(\d+)\s+(\d+.\d+)\s+([-\d.]+)")
}

g16fchkpatterns: Dict[str, re.Pattern] = {
    "charge": re.compile(r"Charge\s+[ICRU]+\s+([\-0-9]+)"),
    "multi": re.compile(r"Multiplicity\s+[ICRU]+\s+([\-0-9]+)"),
    "n_atoms": re.compile(r"Number of atoms\s+[ICRU]+\s+([0-9]+)"),
    "version": re.compile(
        r"Gaussian Version\s+[ICRU]\s+[A-Z=]+\s+[\-0-9]+\s+([a-zA-Z0-9\-\.]+)"
    ),
    "route": re.compile(
        r"Route\s+[ICRU]\s+[A-Z=]+\s+[0-9]+\n([a-zA-Z%0-9\*#\/=\+\- ]+)\n"
    ),
    "total energy": re.compile(r"Total Energy\s+[A-Z]+\s+([\-\+0-9\.E]+)"),
    "scf energy": re.compile(r"SCF Energy\s+[A-Z]+\s+([\-\+0-9\.E]+)"),
    "thermal energy": re.compile(r"Thermal Energy\s+[A-Z]+\s+([\-\+0-9\.E]+)"),
    "thermal enthalpy": re.compile(r"Thermal Enthalpy\s+[A-Z]+\s+([\-\+0-9\.E]+)"),
    "thermal free energy": re.compile(
        r"Thermal Free Energy\s+[A-Z]+\s+([\-\+0-9\.E]+)"
    ),
    "orbital": re.compile(r"(Alpha|Beta) Orbital Energies"),
    "job status": re.compile(r"Job Status\s+[A-Z]+\s+([\-\+0-9\.E]+)"),
    "spin": re.compile(r"S\*\*2\s+R\s+([\+\-\.0-9E]+)"),
    "freq num": re.compile(r"Number of Normal Modes\s+[A-Z]+\s+([\-\+0-9\.E]+)"),
    "dipole": re.compile(r"Dipole Moment\s+[ICRU]\s+N=\s+\d+\n\s+([-\d.E]+)\s+([-\d.E]+)\s+([-\d.E]+)"),
}

xtboutpatterns: Dict[str, re.Pattern] = {
    "version": re.compile(
        r"xtb (version \d+\.\d+\.\d+\s+\([0-9a-z]+\) compiled by ['0-9a-zA-Z@_-]+ on \d{4}-\d{2}-\d{2})"
    ),
    "old version": re.compile(
        r"(Version\s+[0-9\.]+\s+[0-9a-zA-Z]+\s+[\(\)a-zA-Z0-9]+)"
    ),
    "charge": re.compile(r"charge\s+:\s+([-+0-9]+)"),
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
    "coords_new": re.compile(r"([a-zA-Z]+)\s+([\s\-]\d+\.\d+)\s+([\s\-]\d+\.\d+)\s+([\s\-]\d+\.\d+)"),
    "coords_old": re.compile(r"([\s\-]\d+\.\d+)\s+([\s\-]\d+\.\d+)\s+([\s\-]\d+\.\d+)\s+([a-zA-Z]+)"),
    "coords_end": re.compile(r"Bond Distances"),
    "zp correct": re.compile(r"zero point energy\s+([\-\+0-9.]+)"),
    "H correct": re.compile(r"\d+.\d+\s+\d+.\d+E[+-]\d+\s+(\d+.\d+E[+-]\d+)\s+\d+.\d+E[+-]\d+\s+\d+.\d+E[+-]\d+"),
    "G correct": re.compile(r"\d+.\d+\s+\d+.\d+E[+-]\d+\s+\d+.\d+E[+-]\d+\s+\d+.\d+E[+-]\d+\s+(\d+.\d+E[+-]\d+)"),
    "sum H": re.compile(r"TOTAL ENTHALPY\s+([\-\+0-9.]+)"),
    "sum G": re.compile(r"TOTAL FREE ENERGY\s+([\-\+0-9.]+)"),
    "sum E": re.compile(r"TOTAL ENERGY\s+([\-\+0-9.]+)"),
}


def is_metal(number: int):
    if number in (
        1,
        2,
        5,
        6,
        7,
        8,
        9,
        10,
        14,
        15,
        16,
        17,
        18,
        33,
        34,
        35,
        36,
        52,
        53,
        54,
        85,
        86,
        118,
    ):
        return False
    else:
        return True
