"""
Author: TMJ
Date: 2024-10-19 09:57:26
LastEditors: TMJ
LastEditTime: 2024-10-19 21:37:34
Description: 请填写简介
"""

import re
from typing import Dict, Tuple

xtboutpatterns: Dict[str, re.Pattern] = {
    "time": re.compile(
        r"(wall-time:|cpu-time:)\s*(\d+)\s*d,\s*(\d+)\s*h,\s*(\d+)\s*min,\s*(\d+.\d+)\s*sec"
    ),
    "method": re.compile(r"Hamiltonian\s+([GFN\d\-xTBF]+)"),
    "solvent model": re.compile(r"Solvation model:\s+([a-zA-Z\d]+)"),
    "solvent": re.compile(r"Solvent\s+([a-zA-Z\d\-]+)\n"),
    "temperature": re.compile(r"(\d+.\d+)\s*VIB"),
    "version": re.compile(
        r"xtb (version \d+\.\d+\.\d+\s+\([0-9a-z]+\) compiled by ['0-9a-zA-Z@_-]+ on \d{4}-\d{2}-\d{2})"
    ),
    "old version": re.compile(
        r"(Version\s+[0-9\.]+\s+[0-9a-zA-Z]+\s+[\(\)a-zA-Z0-9]+)"
    ),
    "charge": re.compile(r"charge\s+:\s+([-+0-9]+)"),
    "multiplicity": re.compile(r"(--uhf|-u)\s+(\d+)"),
    "total charge": re.compile(r"total charge\s+([\-\+\d]+.\d+)"),
    "parameter": re.compile(r"program call\s+:\s+([0-9a-zA-Z\\/_\-\s\.]+)\n"),
    "n atoms": re.compile(r"number of atoms\s+:\s+(\d+)"),
    "geometric_optimization_state": re.compile(r"GEOMETRY OPTIMIZATION CONVERGED"),
    "orbitals_start": re.compile(r"Orbital Energies and Occupations"),
    "orbitals_match": re.compile(
        r"(\d+|\.\.\.)\s+(\d.\d{4}|\s|\.\.\.)\s+(-*\d+.\d+[^\n]|\.\.\.)\s*(-*\d+.\d+|\.\.\.)\s*(\(HOMO\)|\(LUMO\)|)"
    ),
    "orbitals_end": re.compile(r"HL-Gap\s*(\d+.\d+)\s*Eh\s*(\d+.\d+)\s*eV"),
    "wiberg_start": re.compile(r"#   Z sym  total"),
    "wiberg_end": re.compile(r"Writing (corrected )*topology"),
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
    "coords_type": re.compile(r"coordinate file\s*:\s*[^\s^.]+.([^\s]+)"),
    "optimized_coords": re.compile(r"optimized geometry written to: ([a-zA-Z].xyz)"),
    "attached_coords": re.compile(r"coordinate file\s*:\s*([^\s]+)"),
    "zp correct": re.compile(r"zero point energy\s+([\-\+0-9.]+)"),
    "H correct": re.compile(
        r"\d+.\d+\s+\d+.\d+E[+-]\d+\s+(\d+.\d+E[+-]\d+)\s+\d+.\d+E[+-]\d+\s+\d+.\d+E[+-]\d+"
    ),
    "G correct": re.compile(
        r"\d+.\d+\s+\d+.\d+E[+-]\d+\s+\d+.\d+E[+-]\d+\s+\d+.\d+E[+-]\d+\s+(\d+.\d+E[+-]\d+)"
    ),
    "sum H": re.compile(r"TOTAL ENTHALPY\s*(-*\d+.\d+)\s*Eh"),
    "sum G": re.compile(r"TOTAL FREE ENERGY\s*(-*\d+.\d+)\s*Eh"),
    "sum E": re.compile(r"TOTAL ENERGY\s*(-*\d+.\d+)\s*Eh"),
    "thermal_Cv_S_start": re.compile(r"heat capacity  entropy"),
    "rotation_consts": re.compile(
        r"rotational constants/cm⁻¹ :\s*(\d+.\d+E[+-]\d+)\s*(\d+.\d+E[+-]\d+)\s*(\d+.\d+E[+-]\d+)"
    ),
    "freq_start": re.compile(r"projected vibrational frequencies \(cm⁻¹\)\n"),
    "reduced_masses_start": re.compile(r"reduced masses \(amu\)\n"),
    "ir_start": re.compile(r"IR intensities \(km·mol⁻¹\)\n"),
    "raman_start": re.compile(r"Raman intensities \(amu\)\n"),
    "total_energy": re.compile(r"TOTAL ENERGY\s*(-*\d+.\d+)\s*Eh"),
    "total_E": re.compile(r"(total E|\* total energy)\s*:\s*(-*\d+.\d+)"),
    "gfn1_charges_start": re.compile(
        r"Mulliken/CM5 charges\s*n\(s\)\s*n\(p\)\s*n\(d\)\n"
    ),
    "gfn2_charges_start": re.compile(r"#\s*Z\s*covCN\s*q\s*C6AA\s*α\(0\)\n"),
    "success_tag": re.compile(r"\* finished run on"),
    "VIP": re.compile(r"delta SCC IP \(eV\):\s*(-*\d+.\d+)"),
    "VEA": re.compile(r"delta SCC EA \(eV\):\s*(-*\d+.\d+)"),
    "GEI": re.compile(r"Global electrophilicity index \(eV\):\s*(-*\d+.\d+)"),
    "fukui_start": re.compile(r"#\s*f\(\+\)\s*f\(-\)\s*f\(0\)\n"),
    "fukui_match": re.compile(r"\d+[a-zA-Z]\s*(-*\d+.\d+)\s*(-*\d+.\d+)\s*(-*\d+.\d+)"),
    "energy_convergence": re.compile(r"energy convergence\s*(\d+.\d+[E+\-\d]*)"),
    "grad_convergence": re.compile(r"grad. convergence\s*(\d+.\d+[E+\-\d]*)"),
    "energy_change": re.compile(r"change\s*(-*\d+.\d+[E+\-\d]*)\s*Eh"),
    "grad_norm": re.compile(r"gradient norm\s*:\s*(\d+.\d+[E+\-\d]*)\s*Eh"),
}
