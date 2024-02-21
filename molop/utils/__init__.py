'''
Author: TMJ
Date: 2023-11-02 15:36:39
LastEditors: TMJ
LastEditTime: 2024-02-21 21:07:05
Description: 请填写简介
'''


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
        for para in patterns["link0"].findall(parameter_comment.split("\n")[0])
    }
    route = parameter_comment.split("\n")[1]
    scrf_patt = patterns["scrf_patt"]
    multi_params_patt = patterns["multi_params_patt"]
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


patterns = {
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
}