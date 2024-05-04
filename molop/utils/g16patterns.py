"""
Author: TMJ
Date: 2024-03-23 21:02:58
LastEditors: TMJ
LastEditTime: 2024-03-23 21:03:37
Description: 请填写简介
"""

import re
from typing import Dict, Tuple


def link0_parser(
    link0: str,
) -> Dict:
    """
    Parse link0 from G16 log file.

    Args:
        link0 (str): The link0 section from the G16 log file.

    Returns:
        - link0 (Dict): A dictionary containing parameter-value pairs extracted from the link0 section.
    """
    return {
        f"{para.split('=')[0]}": f"{para.split('=')[1]}"
        for para in g16logpatterns["link0"].findall(link0)
    }


def parameter_comment_parser(
    route: str,
) -> Tuple[Dict, str]:
    """
    Parse parameter comment from G16 log file.

    Args:
        parameter_comment (str): The parameter comment from the G16 log file.

    Returns:
        - route_params (Dict): A dictionary containing parameter-value pairs extracted from the route section of the parameter comment.
        - dieze_tag (str): The dieze tag extracted from the route section of the parameter comment.
    """
    _route = (
        route.replace(" =", "=")
        .replace("= ", "=")
        .replace(" ,", ",")
        .replace(", ", ",")
        .replace("/", " ")
        .lower()
    )

    multi_params_patt = g16logpatterns["multi_params_patt"]
    route_params = {}
    dieze_tag = None

    if _route:
        for tok in _route.split():
            if tok.upper() in ["#", "#N", "#P", "#T"]:
                # does not store # in route to avoid error in input
                dieze_tag = "#N" if tok == "#" else tok
                continue
            else:
                m = re.search(multi_params_patt, tok.strip("#"))
                if m:
                    pars = {}
                    for par in m.group(3).split(","):
                        p = par.split("=")
                        pars[p[0]] = None if len(p) == 1 else p[1].lower()
                    route_params[m.group(1).lower()] = pars
                else:
                    d = tok.strip("#").split("=")
                    route_params[d[0]] = None if len(d) == 1 else d[1].lower()

    return route_params, dieze_tag


def get_solvent_model(route_params: dict) -> str:
    if "scrf" in route_params:
        if isinstance(route_params["scrf"], str):
            return route_params["scrf"]
        elif isinstance(route_params["scrf"], dict):
            model = [
                f"{key}={value}" if value else f"{key}"
                for key, value in route_params["scrf"].items()
                if key != "solvent"
            ]

            if len(model) == 0:
                return "smd"
            else:
                return ",".join(model)
    return ""


def get_solvent(route_params: dict) -> str:
    if "scrf" in route_params:
        if isinstance(route_params["scrf"], str):
            return "water"
        elif isinstance(route_params["scrf"], dict):
            for key, value in route_params["scrf"].items():
                if key == "solvent":
                    return value
    return ""


SHELL_ORBITALS = {
    0: ["S"],
    1: ["PX", "PY", "PZ"],
    -1: ["S", "PX", "PY", "PZ"],
    2: ["D1", "D2", "D3", "D4", "D5", "D6"],
    -2: ["D1", "D2", "D3", "D4", "D5"],
    3: ["F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10"],
    -3: ["F1", "F2", "F3", "F4", "F5", "F6", "F7"],
    4: [
        "G1",
        "G2",
        "G3",
        "G4",
        "G5",
        "G6",
        "G7",
        "G8",
        "G9",
        "G10",
        "G11",
        "G12",
        "G13",
    ],
    -4: ["G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9"],
}

SHELL_START = {0: 1, 1: 2, -1: 2, 2: 3, -2: 3, 3: 4, -3: 4}


def shell_to_orbitals(shell_type, offset):
    return [
        f"{SHELL_START[shell_type] + offset}{x}" for x in SHELL_ORBITALS[shell_type]
    ]


g16logpatterns: Dict[str, re.Pattern] = {
    "n atoms": re.compile(r"NAtoms=\s*(\d+)"),
    "input_coords_start": re.compile(r"Input orientation:"),
    "standard_coords_start": re.compile(r"Standard orientation:"),
    "coords_end": re.compile(
        r"(Basis read|Standard basis|Rotational constants \(GHZ\)|Symmetry turned off)"
    ),
    "coords": re.compile(
        r"\s+\d+\s+(\d+)\s+\d+\s+([\s\-]\d+\.\d+)\s+([\s\-]\d+\.\d+)\s+([\s\-]\d+\.\d+)"
    ),
    "link0": re.compile(r"\%[a-zA-Z]+=[^\s]+"),
    "scrf_patt": re.compile(r"^([sS][cC][rR][fF])\s*=\s*(.+)"),
    "multi_params_patt": re.compile(r"([A-z\d]+)([\s=]+)*\(([a-zA-Z\d,=\-\(\)]+)\)"),
    "rotation_consts": re.compile(
        r"Rotational constants \(GHZ\):\s*(\d+.\d+)\s*(\d+.\d+)\s*(\d+.\d+)"
    ),
    "basis": re.compile(r"Standard basis:\s+([^\s]+)"),
    "solvent_start": re.compile(r"Polarizable Continuum Model \(PCM\)"),
    "solvent_model": re.compile(r"Model\s+:\s+([^\s]+)."),
    "solvent_atomic_radii": re.compile(r"Atomic radii\s+:\s+([^\s]+)."),
    "solvent_type": re.compile(r"Solvent\s*:\s*([^\s]+),"),
    "solvent_eps": re.compile(r"Eps\s*=\s*(\d+.\d+)"),
    "solvent_eps_inf": re.compile(r"Eps\((inf|infinity)\)\s*=\s*(\d+.\d+)"),
    "solvent_end": re.compile(r"Spheres list:"),
    "functional": re.compile(r"SCF Done:  E\(([^\s]+)\)"),
    "scf_energy": re.compile(r"SCF Done:\s*E\(.*\)\s*=\s*([-\s]\d+.\d+)"),
    "mp2-4": re.compile(
        r"E\d[\s\(\)SDTQ]*=\s*-*\d+.\d+D[+-]\d+\s*E*U(MP\d)[\(\)SDTQ]*\s*=\s*(-*\d+.\d+D[+-]\d+)"
    ),
    "mp5": re.compile(r"MP5\s*=\s*(-*\d+.\d+D[+-]\d+)\s*MP5\s*=\s*(-*\d+.\d+D[+-]\d+)"),
    "ccsd": re.compile(r"DE\(Corr\)=\s*[^\s]+\s*E\(Corr\)=\s*([^\s]+)"),
    "ccsd(t)": re.compile(r"CCSD\(T\)\s*=\s*(-*\d+.\d+D[+-]\d+)"),
    "Pseudopotential": re.compile(r"Pseudopotential Parameters"),
    "spins": re.compile(r"<S\*\*2>=([\s\-\d.]+)\s+S=([\s\-\d.]+)"),
    "orbital_start": re.compile(r"Population analysis using the SCF Density"),
    "orbital_end": re.compile(r"Condensed to atoms"),
    "orbital": re.compile(r"(Alpha|Beta)\s*(occ.|virt.)\s*eigenvalues -- (.*)"),
    "mulliken start": re.compile(
        r"(Mulliken charges:|Mulliken atomic charges|Mulliken charges and spin densities:)"
    ),
    "mulliken match": re.compile(r"\d+\s+[A-Z][a-z]?\s+([\s-]\d+.\d+)"),
    "mulliken end": re.compile(r"(Sum of Mulliken )(.*)(charges)\s*=\s*(\D)"),
    "mulliken spin density start": re.compile(r"Mulliken charges and spin densities:"),
    "mulliken spin density match": re.compile(
        r"\d+\s+[A-Z][a-z]?\s+[\s-]\d+.\d+\s+([\s-]\d+.\d+)"
    ),
    "mulliken spin density end": re.compile(
        r"(Sum of Mulliken )(.*)(charges)\s*=\s*(\D)"
    ),
    "apt_start": re.compile(r"APT\s*charges:"),
    "apt_match": re.compile(r"\d+\s+[A-Z][a-z]?\s+([\s-]\d+.\d+)"),
    "apt_end": re.compile(r"Sum of APT charges"),
    "lowdin_start": re.compile(r"Lowdin\s*charges"),
    "lowdin_match": re.compile(r"\d+\s+[A-Z][a-z]?\s+([\s-]\d+.\d+)"),
    "lowdin_end": re.compile(r"Sum of Lowdin charges"),
    "dipole_start": re.compile(r"Dipole moment \(field-independent basis"),
    "dipole": re.compile(r"[XYZ]=\s+([\s-]\d+\.\d+)"),
    "quadrupole_start": re.compile(r"Quadrupole moment \(field-independent basis"),
    "quadrupole": re.compile(r"[XYZ][XYZ]=\s+([\s-]\d+\.\d+)"),
    "octapole_start": re.compile(
        r"Traceless Quadrupole moment \(field-independent basis"
    ),
    "octapole": re.compile(r"[XYZ][XYZ][XYZ]=\s+([\s-]\d+\.\d+)"),
    "hexadecapole_start": re.compile(r"Hexadecapole moment \(field-independent basis"),
    "hexadecapole": re.compile(r"[XYZ][XYZ][XYZ][XYZ]=\s+([\s-]\d+\.\d+)"),
    "hexadecapole_end": re.compile(
        r"N-N=[\s-]\d+.\d+D[+-]\d+\s*E-N=[\s-]\d+.\d+D[+-]\d+\s*KE=[\s-]\d+.\d+D[+-]\d+\s*"
    ),
    "hirshfeld_start": re.compile(
        r"Hirshfeld charges, spin densities, dipoles, and CM5 charges"
    ),
    "hirshfeld": re.compile(
        r"\d\s*[A-Za-z]\s*([-\s]\d+.\d+)\s*([-\s]\d+.\d+)\s*[-\s]\d+.\d+\s*[-\s]\d+.\d+\s*[-\s]\d+.\d+\s*([-\s]\d+.\d+)"
    ),
    "hirshfeld_end": re.compile(
        r"Hirshfeld charges with hydrogens summed into heavy atoms:"
    ),
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
    "freq end": re.compile(r"- Thermochemistry -"),
    "Temperature": re.compile(r"Temperature\s+(\d+.\d+)"),
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
    "scf_energy_alter": re.compile(r"HF=([\-0-9\s+\.]+)\\"),
    "sum energies": re.compile(
        r"Sum of electronic and (thermal Free Energies|thermal Enthalpies|thermal Energies|zero-point Energies)=\s+([\-0-9.]+)"
    ),
    "corrections": re.compile(r"(Zero-point|Thermal) correction(.*)=\s+([\d\.-]+)"),
    "wiberg_start": re.compile(r"Wiberg bond index matrix in the NAO basis"),
    "wiberg_end": re.compile(r"Wiberg bond index, Totals by atom"),
    "atom_atom_overlap_start": re.compile(r"Atom-atom overlap-weighted NAO bond order"),
    "atom_atom_overlap_end": re.compile(
        r"Atom-atom overlap-weighted NAO bond order, Totals by atom"
    ),
    "mo_start": re.compile(r"MO bond order"),
    "mo_end": re.compile(r"MO atomic valencies"),
    "digit": re.compile(r"[\s-]\d+.\d+"),
    "nbo charge start": re.compile(r"Summary of Natural Population Analysis:"),
    "nbo charge end": re.compile(r"=\n\s+\* Total \*"),
    "nbo charge match": re.compile(r"[A-Za-z]+\s+\d+\s+([-\d.]+)"),
    "nbo summary start": re.compile(r"Natural Bond Orbitals \(Summary\):"),
    "nbo summary end": re.compile(r"-\n\s+Total Lewis"),
    "nbo summary match": re.compile(
        r"BD\s+\(\s+\d\)\s+[A-Za-z]+\s+(\d+)\s+-[\sA-Za-z]+\s+(\d+)\s+(\d+.\d+)\s+([-\d.]+)"
    ),
    "tail_start": re.compile(
        r"(Unable to Open any file for archive entry.|Test job not archived.)"
    ),
    "tail_end": re.compile(r"\\@"),
    "tail_match": re.compile(r"(HF|MP2|MP3|MP4[SDTQ]*|CCSD[\(\)T]*)=([\s-]\d+.\d+)"),
    "tail_thermal_match": re.compile(
        r"(ZeroPoint|Thermal|ETot|HTot|GTot)=([\s-]\d+.\d+)"
    ),
    "ratation_consts": re.compile(
        r"Rotational constants \(GHZ\):\s*(\d+.\d+|\*+)\s*(\d+.\d+|\*+)\s*(\d+.\d+|\*+)"
    ),
    "thermal_Cv_S_start": re.compile(r"E \(Thermal\)\s*CV\s*S"),
    "electronic_spatial_extent": re.compile(
        r"Electronic spatial extent \(au\):  <R\*\*2>=\s*(\d+.\d+)"
    ),
    "isotropic_polarizability": re.compile(r"(\d+.\d+)\s*Bohr\*\*3."),
}


g16fchkpatterns: Dict[str, re.Pattern] = {
    "charge": re.compile(r"Charge\s+[ICRU]+\s+([\-0-9]+)"),
    "multi": re.compile(r"Multiplicity\s+[ICRU]+\s+([\-0-9]+)"),
    "n_atoms": re.compile(r"Number of atoms\s+[ICRU]+\s+([0-9]+)"),
    "version": re.compile(r"Gaussian Version\s*[ICRU]\s*[A-Z=]+\s*[\-0-9]+\n([^\s]+)"),
    "route": re.compile(
        r"Route\s+[ICRU]\s+[A-Z=]+\s+[0-9]+\n([a-zA-Z%0-9\*#\/=\+\- \n\(\),]+)\nCharge"
    ),
    "atomic_number_start": re.compile(r"Atomic numbers\s*[ICRU]\s*N=\s*(\d+)"),
    "atomic_number_end": re.compile(r"Nuclear charges\s*[ICRU]\s*N=\s*(\d+)"),
    "coords_start": re.compile(r"Current cartesian coordinates\s*[ICRU]\s*N=\s*(\d+)"),
    "coords_end": re.compile(r"Number of symbols in /Mol/"),
    "total energy": re.compile(r"Total Energy\s*[ICRU]\s*([\s-]\d+.\d+E[+-]\d+)"),
    "scf energy": re.compile(r"SCF Energy\s*[ICRU]\s*([\s-]\d+.\d+E[+-]\d+)"),
    "cluster energy": re.compile(r"Cluster Energy\s*[ICRU]\s*([\s-]\d+.\d+E[+-]\d+)"),
    "mp2-4": re.compile(r"(MP[234])[SDTQ]* Energy\s*[ICRU]\s*([\s-]\d+.\d+E[+-]\d+)"),
    "spin": re.compile(r"S\*\*2\s+R\s+([\+\-\.0-9E]+)"),
    "job status": re.compile(r"Job Status\s+[A-Z]+\s+([\-\+0-9\.E]+)"),
    "thermal energy": re.compile(r"Thermal Energy\s+[A-Z]+\s+([\-\+0-9\.E]+)"),
    "thermal enthalpy": re.compile(r"Thermal Enthalpy\s+[A-Z]+\s+([\-\+0-9\.E]+)"),
    "thermal free energy": re.compile(
        r"Thermal Free Energy\s+[A-Z]+\s+([\-\+0-9\.E]+)"
    ),
    "alpha_elec": re.compile(r"Number of alpha electrons\s*[ICRU]\s*(\d+)"),
    "alpha_start": re.compile(r"Alpha Orbital Energies\s*[ICRU]\s*N=\s*(\d+)"),
    "alpha_end": re.compile(r"Alpha MO coefficients\s*[ICRU]\s*N=\s*(\d+)"),
    "beta_elec": re.compile(r"Number of beta electrons\s*[ICRU]\s*(\d+)"),
    "beta_start": re.compile(r"Beta Orbital Energies\s*[ICRU]\s*N=\s*(\d+)"),
    "beta_end": re.compile(r"Beta MO coefficients\s*[ICRU]\s*N=\s*(\d+)"),
    "mulliken_start": re.compile(r"Mulliken Charges\s*[ICRU]\s*N=\s*(\d+)"),
    "gradient_start": re.compile(r"Cartesian Gradient\s*[ICRU]\s*N=\s*(\d+)"),
    "gradient_end": re.compile(r"Nonadiabatic coupling\s*[ICRU]\s*N=\s*(\d+)"),
    "freq num": re.compile(r"Number of Normal Modes\s+[A-Z]+\s+([\-\+0-9\.E]+)"),
    "vib_e2_start": re.compile(r"Vib-E2\s*[ICRU]\s*N=\s*(\d+)"),
    "vib_mode_start": re.compile(r"Vib-Modes\s*[ICRU]\s*N=\s*(\d+)"),
    "vib_mode_end": re.compile(r"ClPar MaxAn"),
    "dipole": re.compile(
        r"Dipole Moment\s+[ICRU]\s+N=\s+\d+\n\s+([-\d.E]+)\s+([-\d.E]+)\s+([-\d.E]+)"
    ),
    "int_digits": re.compile(r"\d+"),
    "float_digits": re.compile(r"[\s-]\d+.\d+E[+-]\d+"),
}
