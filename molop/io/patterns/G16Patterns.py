"""
Author: TMJ
Date: 2025-02-16 19:20:04
LastEditors: TMJ
LastEditTime: 2025-02-17 12:41:34
Description: 请填写简介
"""

import re
from typing import Any, Dict, Tuple

from molop.io.base_models.SearchPattern import MolOPPattern

SEMI_EMPIRICAL_METHODS = ("am1", "pm3", "pm6", "pm7", "pddg", "indo", "cndo", "pm3mm")


class G16InputPatterns:
    """
    G16 input file patterns
    """

    OPTIONS = MolOPPattern(
        content_pattern=r"^\s*(%[a-zA-Z0-9]+)=([^\s]+)",
        content_repeat=0,
        end_pattern=r"^\s*\n",
    )
    ROUTE = MolOPPattern(
        content_pattern=r"^\s*(#[a-zA-Z0-9\(\)=,\/\\\s\n-]+)\n\n", end_pattern=r"^\s*\n"
    )
    TITLE = MolOPPattern(end_pattern=r"^\s*\n", content_pattern=r"^([^\n]+)")
    CHARGE_MULTIPLICITY = MolOPPattern(
        content_pattern=r"^\s*(-?\d+)\s+(\d+)\s*\n",
        description="The charge and multiplicity of the Gaussian calculation. link 1",
    )
    ATOMS = MolOPPattern(
        content_pattern=r"^\s*([A-Z][a-z]*)(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)",
        content_repeat=0,
    )


class G16LogPatterns:
    """
    G16 log file patterns
    """

    VERSION = MolOPPattern(
        start_pattern=r"^\s\*+\n",
        end_pattern=r"^\s\*+\n",
        end_offset=1,
        content_pattern=r"(Gaussian\s+\d+:\s*[^\s]*\s*\d+-[A-Za-z]{3}-\d{4}\n\s+\d+-[A-Za-z]{3}-\d{4})",
        description="The exact version of Gaussian used. link 1",
    )
    OPTIONS = MolOPPattern(
        start_pattern=r"^\s\*+\n",
        end_pattern=r"^\s*-+\n",
        content_pattern=r"^\s*(%[a-zA-Z0-9]+)=([^\s]+)",
        content_repeat=0,
        description="The options used in the Gaussian calculation. Such as: MEMORY, CPU, etc. link 1",
    )
    KEYWORDS = MolOPPattern(
        start_pattern=r"^\s*-+\n",
        end_pattern=r"^\s*-+\n",
        end_offset=1,
        content_pattern=r"^\s*(#[a-zA-Z0-9\(\)=,\/\\\s\n\-]+)\n\s*\-+",
        description="The keywords used in the Gaussian calculation. link 1",
    )
    TITLE = MolOPPattern(
        start_pattern=r"^\s*-+\n",
        start_offset=1,
        end_pattern=r"^\s*-+\n",
        end_offset=1,
        content_pattern=r"----\n\s*((.\n*)+)\n\s*----",
        description="The title card of the Gaussian calculation task. link 101",
    )
    CHARGE_MULTIPLICITY = MolOPPattern(
        content_pattern=r"Charge\s*=\s*([-\+\d]+)\s+Multiplicity\s*=\s*(\d+)",
        description="The charge and multiplicity of the Gaussian calculation. link 101",
    )
    INITIAL_INPUT_COORDS = MolOPPattern(
        start_pattern="Symbolic Z-matrix:",
        strat_regex=False,
        end_pattern=r"^\s*\n",
        content_pattern=r"^\s*([A-Z][a-z]*)\s+(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)",
        content_repeat=0,
        description="The initial input coordinates of the Gaussian calculation. link 101",
    )
    INPUT_COORDS = MolOPPattern(
        start_pattern="Input orientation:",
        strat_regex=False,
        end_pattern=r"^\s*-+\n",
        end_offset=2,
        content_pattern=r"\s+\d+\s+(\d+)\s+\d+\s+(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)",
        content_repeat=0,
        description="The input coordinates of the Gaussian calculation. link 202",
    )
    STANDARD_COORDS = MolOPPattern(
        start_pattern="Standard orientation:",
        strat_regex=False,
        end_pattern=r"^\s*-+\n",
        end_offset=2,
        content_pattern=r"\s+\d+\s+(\d+)\s+\d+\s+(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)",
        content_repeat=0,
        description="The standard coordinates of the Gaussian calculation. link 202",
    )
    COORDS = MolOPPattern(
        end_pattern=r"(Basis read|Standard basis|Rotational constants \(GHZ\)|Symmetry turned off|The archive entry for this job was punched.)",
        content_pattern=r"^\s*\d+\s+(\d+)\s+\d+\s+(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)",
        content_repeat=0,
        description="The coordinates of the atoms in specific frame. link 202",
    )
    ROTATIONAL_CONST = MolOPPattern(
        content_pattern=r"^\s*Rotational constants \(GHZ\):\s*(\d+.\d+|\*+)\s*(\d+.\d+|\*+)\s*(\d+.\d+|\*+)",
        description="The rotational constants of the coordinates. link 202",
    )
    BASIS_SET = MolOPPattern(
        content_pattern=r"^\s*Standard basis:\s+([^\s]+)",
        description="The basis set used in the Gaussian calculation. link 301",
    )
    PSEUDOPOTENTIAL_PARAMETERS = MolOPPattern(
        start_pattern=r"^\s*={2,}\n",
        end_pattern=r"^\s*={2,}\n",
        end_offset=2,
        description="The parameters of the pseudopotential used in the Gaussian calculation. link 301",
    )
    DISPERSION_CORRECTION = MolOPPattern(
        content_pattern=r"^\s*R6Disp:\s*(.+)\s*Dispersion energy=(\s*-?\d+\.\d*)\s*Hartrees.",
        description="The dispersion correction used in the Gaussian calculation. link 301",
    )
    SOLVENT_PARAMETERS = MolOPPattern(
        start_pattern="Polarizable Continuum Model (PCM)",
        strat_regex=False,
        end_pattern=r"^\s*-{2,}\n",
        description="The solvent parameters used in the Gaussian calculation. link 301",
    )
    SOLVENT_MODEL = MolOPPattern(
        content_pattern=r"Model\s*:\s*(.+)\.",
        description="The solvent model used in the Gaussian calculation. link 301",
    )
    SOLVENT_ATOM_RADII = MolOPPattern(
        content_pattern=r"Atomic radii\s*:\s*(.+)\.",
        description="The radii of the atoms in the solvent model. link 301",
    )
    SOLVENT_TYPE = MolOPPattern(
        content_pattern=r"Solvent\s*:\s*(.+),\s",
        description="The type of the solvent used in the Gaussian calculation. link 301",
    )
    SOLVENT_EPS = MolOPPattern(
        content_pattern=r"Eps\s*=\s*(\d+.\d*)",
        description="The solvent dielectric constant used in the Gaussian calculation. link 301",
    )
    SOLVENT_EPS_INF = MolOPPattern(
        content_pattern=r"Eps\((inf|infinity)\)\s*=\s*(\d+.\d*)",
        description="The solvent dielectric constant at infinity used in the Gaussian calculation. link 301",
    )
    SCF_ENERGIES = MolOPPattern(
        start_pattern=r"^\s*SCF Done:",
        description="The SCF energies of the Gaussian calculation. link 502",
    )
    SCF_ENERGY_AND_FUNCTIONAL = MolOPPattern(
        content_pattern=r"^\s*SCF Done:\s*E\((.*)\)\s*=\s*(\s*-?\d+\.\d*)",
        description="The SCF energy and functional used in the Gaussian calculation. link 502",
    )
    SPIN_SPIN_SQUERE = MolOPPattern(
        content_pattern=r"<S\*\*2>=(\s*-?\d+\.\d*)\s+S=(\s*-?\d+\.\d*)",
        description="The total spin and spin squere exactly after the SCF calculation. link 502",
    )
    ENERGY_MP2_4 = MolOPPattern(
        content_pattern=r"E\d[\s\(\)SDTQ]*=\s*-*\d+.\d+D[+-]\d+\s*E*U(MP\d)[\(\)SDTQ]*\s*=\s*(-?\d+.\d+D[+-]\d+)",
        content_repeat=0,
        description="The MP2-4 energy of the Gaussian calculation. link 804 for mp2 and link 913 for mp3-4",
    )
    ENERGY_MP5 = MolOPPattern(
        content_pattern=r"MP5\s*=\s*(-*\d+.\d+D[+-]\d+)\s*MP5\s*=\s*(-?\d+.\d+D[+-]\d+)",
        description="The MP5 energy of the Gaussian calculation. link 913",
    )
    ENERGY_CCSD = MolOPPattern(
        content_pattern=r"Wavefunction amplitudes converged. E\(Corr\)=\s*(-*\d+.\d*)",
        description="The CCSD energy of the Gaussian calculation. link 913",
    )
    ENERGY_CCSD_T = MolOPPattern(
        content_pattern=r"CCSD\(T\)\s*=\s*(-?\d+.\d+D[+-]\d+)",
        description="The CCSD(T) energy of the Gaussian calculation. link 913",
    )
    ISOTROPIC_POLARIZABILITY = MolOPPattern(
        start_pattern="Isotropic polarizability for W=",
        strat_regex=False,
        end_pattern=r"^\s*Isotropic polarizability for W=\s*\d+.\d+\s*(\d+.\d+)\s*Bohr\*\*3",
        content_pattern=r"\d+.\d+\s*(\d+.\d+)",
        description="The isotropic polarizability of the Gaussian calculation. link 1002",
    )
    POPULATION_ANALYSIS = MolOPPattern(
        start_pattern=r"^\s*Population analysis using the (SCF|CC) [Dd]ensity.",
        end_pattern=r"^\s*N-N=.*",
        description="The population analysis of the Gaussian calculation. link 601",
    )
    MOLECULAR_ORBITALS_SYMMETRY_ALPHA = MolOPPattern(
        start_pattern="Alpha Orbitals:",
        strat_regex=False,
        end_pattern="Beta  Orbitals:",
        end_regex=False,
        content_pattern=r"^\s*(Occupied|Virtual|)\s*((\([A-Z0-9]+\)\s*){1,12})\n",
        content_repeat=0,
        description="The alpha molecular orbitals symmetry of the Gaussian calculation. link 601",
    )
    MOLECULAR_ORBITALS_SYMMETRY_BETA = MolOPPattern(
        start_pattern="Beta  Orbitals:",
        strat_regex=False,
        end_pattern=r"^\s*The electronic state is (.*)\.",
        content_pattern=r"^\s*(Occupied|Virtual|)\s*((\([A-Z0-9]+\)\s*){1,12})\n",
        content_repeat=0,
        description="The beta molecular orbitals symmetry of the Gaussian calculation. link 601",
    )
    MOLECULAR_ORBITALS_SYMMETRY = MolOPPattern(
        start_pattern="Orbital symmetries:",
        strat_regex=False,
        end_pattern=r"^\s*The electronic state is (.*)\.",
        content_pattern=r"^\s*(Occupied|Virtual|)\s*((\([A-Z0-9]+\)\s*){1,12})\n",
        content_repeat=0,
        description="The molecular orbitals symmetry of the Gaussian calculation. link 601",
    )
    ELECTRONIC_STATE = MolOPPattern(
        content_pattern=r"^\s*The electronic state is (.*)\.",
        description="The electronic state of the Gaussian calculation. link 601",
    )
    MOLECULAR_ORBITALS = MolOPPattern(
        content_pattern=r"^\s*(Alpha|Beta)\s*(occ.|virt.)\s*eigenvalues -- (.*)",
        content_repeat=0,
        description="The molecular orbitals of the Gaussian calculation. link 601",
    )
    MULLIKEN_POPULATION = MolOPPattern(
        start_pattern=r"(Mulliken charges:|Mulliken atomic charges|Mulliken charges and spin densities:)",
        end_pattern=r"(Sum of Mulliken )(.*)(charges)\s*=\s*(\D)",
        content_pattern=r"\d+\s+[A-Z][a-z]?\s+(\s*-?\d+\.\d*)",
        content_repeat=0,
        description="The Mulliken population analysis of the Gaussian calculation. link 601",
    )
    MULLIKEN_SPIN_DENSITY = MolOPPattern(
        start_pattern="Mulliken charges and spin densities:",
        strat_regex=False,
        end_pattern=r"(Sum of Mulliken )(.*)(charges)\s*=\s*(\D)",
        content_pattern=r"\d+\s+[A-Z][a-z]?(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)",
        content_repeat=0,
        description="The Mulliken spin density analysis of the Gaussian calculation. link 601",
    )
    APT_POPULATION = MolOPPattern(
        start_pattern="APT charges:",
        strat_regex=False,
        end_pattern="Sum of APT charges",
        end_regex=False,
        content_pattern=r"\d+\s+[A-Z][a-z]?\s+(\s*-?\d+\.\d*)",
        content_repeat=0,
        description="The APT population analysis of the Gaussian calculation. link 601",
    )
    LOWDIN_POPULATION = MolOPPattern(
        start_pattern="Lowdin charges",
        strat_regex=False,
        end_pattern="Sum of Lowdin charges",
        end_regex=False,
        content_pattern=r"\d+\s+[A-Z][a-z]?\s+(\s*-?\d+\.\d*)",
        content_repeat=0,
        description="The Lowdin population analysis of the Gaussian calculation. link 601",
    )
    ELECTRONIC_SPATIAL_EXTENT = MolOPPattern(
        content_pattern=r"Electronic spatial extent \(au\):  <R\*\*2>=\s*(\d+\.\d*)",
        description="The electronic spatial extent (<R**2>) of the Gaussian calculation. link 601",
    )
    DIPOLE_MOMENT = MolOPPattern(
        start_pattern="Dipole moment (field-independent basis",
        strat_regex=False,
        content_pattern=r"[XYZ]=\s+(\s*-?\d+\.\d*)",
        content_repeat=3,
        description="The dipole moment of the Gaussian calculation. link 601",
    )
    QUADRUPOLE_MOMENT = MolOPPattern(
        start_pattern="Quadrupole moment (field-independent basis",
        strat_regex=False,
        content_pattern=r"\s[XYZ]{2}=\s+(\s*-?\d+\.\d*)",
        content_repeat=6,
        description="The quadrupole moment of the Gaussian calculation. link 601",
    )
    TRACELESS_QUADRUPOLE_MOMENT = MolOPPattern(
        start_pattern="Traceless Quadrupole moment (field-independent basis",
        strat_regex=False,
        content_pattern=r"\s[XYZ]{2}=\s+(\s*-?\d+\.\d*)",
        content_repeat=6,
        description="The traceless quadrupole moment of the Gaussian calculation. link 601",
    )
    OCTAPOLE_MOMENT = MolOPPattern(
        start_pattern="Octapole moment (field-independent basis",
        strat_regex=False,
        content_pattern=r"\s[XYZ]{3}=\s+(\s*-?\d+\.\d*)",
        content_repeat=10,
        description="The octapole moment of the Gaussian calculation. link 601",
    )
    HEXADECAPOLE_MOMENT = MolOPPattern(
        start_pattern="Hexadecapole moment (field-independent basis",
        strat_regex=False,
        content_pattern=r"\s[XYZ]{4}=\s+(\s*-?\d+\.\d*)",
        content_repeat=15,
        description="The hexadecapole moment of the Gaussian calculation. link 601",
    )
    HIRSHFELD_POPULATION = MolOPPattern(
        start_pattern="Hirshfeld charges, spin densities, dipoles, and CM5 charges",
        strat_regex=False,
        end_pattern=r"^\s*Hirshfeld charges( and spin densities|) with hydrogens",
        content_pattern=r"^\s*\d+\s+[A-Z][a-z]?\s+(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)",
        content_repeat=0,
        description="The Hirshfeld population analysis of the Gaussian calculation. link 601",
    )
    EXACT_POLARIZABILITY = MolOPPattern(
        content_pattern=r"Exact\s*polarizability:(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)",
        description="The exact polarizability of the Gaussian calculation. link 601",
    )
    APPROX_POLARIZABILITY = MolOPPattern(
        content_pattern=r"Approx\s*polarizability:(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)",
        description="The approximate polarizability of the Gaussian calculation. link 601",
    )
    ISOTROPIC_FERMI_CONTACT_COUPLING = MolOPPattern(
        start_pattern="Isotropic Fermi Contact Couplings",
        strat_regex=False,
        content_pattern=r"\s*\d+\s*[A-Z][a-z]?\((\d+)\)(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)",
        description="The isotropic Fermi contact coupling of the Gaussian calculation. link 601",
    )
    ESP_POPULATION = MolOPPattern(
        start_pattern="ESP charges:",
        strat_regex=False,
        end_pattern="Sum of ESP charges",
        end_regex=False,
        content_pattern=r"\d+\s+[A-Z][a-z]?\s+(\s*-?\d+\.\d*)",
        content_repeat=0,
        description="The ESP population analysis of the Gaussian calculation. link 601",
    )
    DIPOLE_BEFORE_FORCE = MolOPPattern(
        start_pattern="Dipole        =",
        strat_regex=False,
        content_pattern=r"(\s*-*\d+.\d+D[+-]\d+)",
        content_repeat=3,
        description="The dipole moment before the force calculation. link 716",
    )
    POLARIZIABILITIES_BEFORE_FORCE = MolOPPattern(
        start_pattern="Polarizability=",
        strat_regex=False,
        content_pattern=r"(\s*-*\d+.\d+D[+-]\d+)",
        content_repeat=6,
        description="The polarizabilities before the force calculation. link 716",
    )
    # TODO: add Bond order analysis patterns, link 607
    # TODO: add TDDFT patterns, link 914
    TDDFT_ORBITALS = MolOPPattern(description="TODO, link 914")
    FREQUENCY_ANALYSIS = MolOPPattern(
        start_pattern="Harmonic frequencies (cm**-1)",
        strat_regex=False,
        end_pattern=r" \n",
        end_regex=False,
        description="The frequency analysis of the Gaussian calculation. link 716",
    )
    DIAGONAL_VIBRATIONAL_POLARIZABILITY = MolOPPattern(
        content_pattern=r"^\s*Diagonal vibrational polarizability:\n(\s*-?\d*\.\d*)(\s*-?\d*\.\d*)(\s*-?\d*\.\d*)",
        description="The diagonal vibrational polarizability of the Gaussian calculation. link 716",
    )
    FREQUENCIES = MolOPPattern(
        content_pattern=r"Frequencies -- ((\s+-?\d+\.\d*){1,3})",
        content_repeat=0,
        description="The frequencies of each vibration in the frequency analysis of the Gaussian calculation. link 716",
    )
    FREQUENCIES_REDUCED_MASS = MolOPPattern(
        content_pattern=r"Red. masses -- ((\s+-?\d+\.\d*){1,3})",
        content_repeat=0,
        description="The reduced masses of each vibration in the frequency analysis of the Gaussian calculation. link 716",
    )
    FREQUENCIES_FORCE_CONSTANTS = MolOPPattern(
        content_pattern=r"Frc consts  -- ((\s+-?\d+\.\d*){1,3})",
        content_repeat=0,
        description="The force constants of each vibration in the frequency analysis of the Gaussian calculation. link 716",
    )
    FREQUENCIES_IR_INTENSITIES = MolOPPattern(
        content_pattern=r"IR Inten    -- ((\s+-?\d+\.\d*){1,3})",
        content_repeat=0,
        description="The IR intensities of each vibration in the frequency analysis of the Gaussian calculation. link 716",
    )
    FREQUENCIES_MODE = MolOPPattern(
        end_pattern="Thermochemistry",
        end_regex=False,
        content_pattern=r"^\s+\d+\s+\d+\s*((\s*-?\d+\.\d*){3,9})",
        content_repeat=0,
        description="The mode of each vibration in the frequency analysis of the Gaussian calculation. link 716",
    )
    TEMPEREATURE_PRESSURE = MolOPPattern(
        start_pattern="- Thermochemistry",
        strat_regex=False,
        end_pattern=r"^\s*Temperature\s*(\d+\.\d*)\s*Kelvin\.\s*Pressure\s*(\d+\.\d*)\s*Atm\.",
        content_pattern=r"^\s*Temperature\s*(\d+\.\d*)\s*Kelvin\.\s*Pressure\s*(\d+\.\d*)\s*Atm\.",
        description="The temperature and pressure of the Gaussian calculation. link 716",
    )
    ROTATIONAL_TEMPERATURE = MolOPPattern(
        content_pattern=r"^\s*Rotational temperatures \(Kelvin\)(\s*\d+\.\d*)(\s*\d+\.\d*)(\s*\d+\.\d*)",
        description="The rotational temperatures of the Gaussian calculation. link 716",
    )
    ROTATIONAL_CONST_IN_FREQUENCY_ANALYSIS = MolOPPattern(
        content_pattern=r"^\s*Rotational constants \(GHZ\):(\s*\d+\.\d*)(\s*\d+\.\d*)(\s*\d+\.\d*)",
        description="The rotational constants of the Gaussian calculation. link 716",
    )
    VIBRATIONAL_TEMPERATURE = MolOPPattern(
        content_pattern=r"^\s*(Vibrational temperatures:|\(Kelvin\)|)(\s*\d+\.\d*){1,5}\n",
        content_repeat=0,
        description="The vibrational temperatures of the Gaussian calculation. link 716",
    )
    THERMOCHEMISTRY_PART = MolOPPattern(
        start_pattern="Zero-point correction",
        strat_regex=False,
        end_pattern=r"^\s*Rotational(\s*-*\d+.\d+D[+-]\d+)\s*(-?\d+\.\d*)\s*(-?\d+\.\d*)",
        description="The thermochemistry part of the Gaussian calculation. link 716",
    )
    THERMOCHEMISTRY_CORRECTION = MolOPPattern(
        content_pattern=r"^\s*(Zero-point|Thermal) correction(.*)=\s*(-?\d+\.\d*)",
        content_repeat=0,
        description="The thermochemistry corrections of the Gaussian calculation. link 716",
    )
    THERMOCHEMISTRY_SUM = MolOPPattern(
        content_pattern=r"^\s*Sum of electronic and (thermal Free Energies|thermal Enthalpies|thermal Energies|zero-point Energies)=\s*(-?\d+.\d*)",
        content_repeat=0,
        description="The thermochemistry sums of the Gaussian calculation. link 716",
    )
    THERMOCHEMISTRY_CV_S = MolOPPattern(
        start_pattern="E (Thermal)             CV                S",
        strat_regex=False,
        content_pattern=r"^\s*Total(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)",
        description="The thermochemistry CV and S of the Gaussian calculation. link 716",
    )
    FORCES_IN_CARTESIAN = MolOPPattern(
        start_pattern="Center     Atomic                   Forces (Hartrees/Bohr)",
        strat_regex=False,
        end_pattern=r"^\s*(-){3,}",
        end_offset=1,
        content_pattern=r"^\s*\d+\s+\d+(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)",
        content_repeat=0,
        description="The forces in Cartesian coordinates of the Gaussian calculation. link 716",
    )
    HESSIAN_IN_CARTESIAN = MolOPPattern(
        start_pattern="Force constants in Cartesian coordinates",
        strat_regex=False,
        end_pattern=r"^\s*FormGI is forming",
        content_pattern=r"^\s*(\d+)((\s*-?\d+.\d+D[+-]\d+){1,5})",
        content_repeat=0,
        description="The Hessian in Cartesian coordinates of the Gaussian calculation. link 716",
    )
    BERNY_STATE_MAJOR_PART = MolOPPattern(
        start_pattern="Item               Value     Threshold  Converged?",
        strat_regex=False,
        end_pattern=r"(?r)^\s*^\s*Leave Link  103.*MaxMem=\s*(\d+)\s*cpu:\s*(\d+\.\d*)",
        description="The Berny optimization state. Patterns showed when #p used. link 103",
    )
    BERNY_STATE_BACKUP_PART = MolOPPattern(
        start_pattern="Item               Value     Threshold  Converged?",
        strat_regex=False,
        end_pattern=r"^\s*(Grad){3,}",
        description="The Berny optimization state. link 103",
    )
    BERNY_STATE = MolOPPattern(
        content_pattern=r"(?r)^\s*(Maximum\s+Force|RMS\s+Force|Maximum\s+Displacement|RMS\s+Displacement)\s+(\d+\.\d*)\s+(\d+\.\d*)\s+(NO|YES)",
        content_repeat=4,
        description="The Berny optimization state. link 103",
    )
    ENERGY_CHANGE = MolOPPattern(
        content_pattern=r"(?r)^\s*^\s*Predicted change in Energy=(\s*-?\d+.\d+D[+-]\d+)",
        description="The energy change of the Berny optimization state. link 103",
    )
    BERNY_CONCLUSION = MolOPPattern(
        content_pattern=r"(?r)^\s*^\s*Optimization completed\.",
        description="The Berny optimization conclusion. link 103",
    )
    ELECTRIC_DIPOLE_PART = MolOPPattern(
        start_pattern="Electric dipole moment (input orientation):",
        strat_regex=False,
        end_pattern=r"^\s*(-){3,}",
        description="The electric dipole moment part of the Gaussian calculation. link 9999",
    )
    ELECTRIC_DIPOLE_MOMENT = MolOPPattern(
        start_pattern="Electric dipole moment (input orientation):",
        strat_regex=False,
        content_pattern=r"^\s*(Tot|x|y|z)(\s*-?\d+\.\d*D[+-]\d+)(\s*-?\d+\.\d*D[+-]\d+)(\s*-?\d+\.\d*D[+-]\d+)",
        content_repeat=4,
        description="The electric dipole moment of the Gaussian calculation. link 9999",
    )
    DIPOLE_POLARIZABILITY = MolOPPattern(
        start_pattern="Dipole polarizability, Alpha (input orientation)",
        strat_regex=False,
        end_pattern=r"^\s*(-){3,}",
        content_pattern=r"^\s*(iso|aniso|xx|yx|yy|zx|zy|zz)(\s*-?\d+\.\d*D[+-]\d+)(\s*-?\d+\.\d*D[+-]\d+)(\s*-?\d+\.\d*D[+-]\d+)",
        content_repeat=8,
        description="The dipole polarizability of the Gaussian calculation. link 9999",
    )
    ARCHIVE_TAIL = MolOPPattern(
        start_pattern=r"^\s*1[\\|]1[\\|]GINC",
        end_pattern=r"(\\@\n|@\n)",
        description="The archive tail of the Gaussian calculation. link 9999",
    )
    JOB_TYPE_IN_ARCHIVE_TAIL = MolOPPattern(
        start_pattern="\\",
        strat_regex=False,
        start_offset=2,
        end_pattern="\\",
        end_regex=False,
        end_offset=1,
        content_pattern=r"\\(.*)\\",
        description="The job type in the 4th block of the archive tail of the Gaussian calculation. link 9999",
    )
    FUNCTIONAL_IN_ARCHIVE_TAIL = MolOPPattern(
        start_pattern="\\",
        strat_regex=False,
        end_pattern="\\",
        end_regex=False,
        end_offset=1,
        content_pattern=r"\\(.*)\\",
        description="The functional in the 5th block of the archive tail of the Gaussian calculation. link 9999",
    )
    BASIS_SET_IN_ARCHIVE_TAIL = MolOPPattern(
        start_pattern="\\",
        strat_regex=False,
        end_pattern="\\",
        end_regex=False,
        end_offset=1,
        content_pattern=r"\\(.*)\\",
        description="The basis set in the 6th block of the archive tail of the Gaussian calculation. link 9999",
    )
    KEYWORDS_IN_ARCHIVE_TAIL = MolOPPattern(
        start_pattern="\\\\",
        strat_regex=False,
        end_pattern="\\\\",
        end_regex=False,
        end_offset=1,
        content_pattern=r"\\\\(.*)\\\\",
        description="The keywords in the 12th block of the archive tail of the Gaussian calculation. link 9999",
    )
    TITLE_IN_ARCHIVE_TAIL = MolOPPattern(
        start_pattern="\\\\",
        strat_regex=False,
        end_pattern="\\\\",
        end_regex=False,
        end_offset=1,
        content_pattern=r"\\\\(.*)\\\\",
        description="The title in the 14th block of the archive tail of the Gaussian calculation. link 9999",
    )
    CHARGE_SPIN_MULTIPLICITY_IN_ARCHIVE_TAIL = MolOPPattern(
        start_pattern="\\\\",
        strat_regex=False,
        end_pattern="\\",
        end_regex=False,
        end_offset=2,
        content_pattern=r"\\(-?\d+),(-?\d+)\\",
        description="The charge and spin multiplicity in the 16th block of the archive tail of the Gaussian calculation. link 9999",
    )
    COORS_IN_ARCHIVE_TAIL = MolOPPattern(
        start_pattern="\\",
        strat_regex=False,
        end_pattern="\\\\",
        end_regex=False,
        content_pattern=r"\\([A-Z][a-z]?),(-?\d+\.\d*),(-?\d+\.\d*),(-?\d+\.\d*)",
        content_repeat=0,
        description="The coordinates in the 17th block of the archive tail of the Gaussian calculation. link 9999",
    )
    VERSION_IN_ARCHIVE_TAIL = MolOPPattern(
        start_pattern="\\\\",
        strat_regex=False,
        end_pattern="\\",
        end_regex=False,
        end_offset=2,
        content_pattern=r"\\\\Version=(.*)\\",
        description="The version in the archive tail of the Gaussian calculation. link 9999",
    )
    # TODO: electronic state
    ENERGIES_IN_ARCHIVE_TAIL = MolOPPattern(
        content_pattern=r"(HF|MP2|MP3|MP4[SDTQ]*|CCSD[\(\)T]*)=(\s*-?\d+\.\d*)",
        content_repeat=0,
        description="The energies in the archive tail of the Gaussian calculation. link 9999",
    )
    THERMOCHEMISTRY_IN_ARCHIVE_TAIL = MolOPPattern(
        content_pattern=r"(ZeroPoint|Thermal|ETot|HTot|GTot)=(\s*-?\d+\.\d*)",
        description="The thermochemistry in the archive tail of the Gaussian calculation. link 9999",
    )
    DIPOLE_IN_ARCHIVE_TAIL = MolOPPattern(
        content_pattern=r"Dipole\s*=(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)(\s*-?\d+\.\d*)",
        description="The dipole in the archive tail of the Gaussian calculation. link 9999",
    )
    POLAR_IN_ARCHIVE_TAIL = MolOPPattern(
        content_pattern=r"Polar\s*=\s*(-?\d+\.\d*)(-?\d+\.\d*)(-?\d+\.\d*)(-?\d+\.\d*)(-?\d+\.\d*)(-?\d+\.\d*)",
        description="The polarizability in the archive tail of the Gaussian calculation. link 9999",
    )
    QUADRUPOLE_IN_ARCHIVE_TAIL = MolOPPattern(
        content_pattern=r"Quadrupole\s*=\s*(-?\d+\.\d*)(-?\d+\.\d*)(-?\d+\.\d*)(-?\d+\.\d*)(-?\d+\.\d*)(-?\d+\.\d*)",
        description="The quadrupole in the archive tail of the Gaussian calculation. link 9999",
    )
    HESSIAN_IN_ARCHIVE_TAIL = MolOPPattern(
        start_pattern=r"NImag=\d+\\\\",
        content_pattern=r"(-?\d+\.\d*)",
        end_pattern=r"\\\\",
        content_repeat=0,
        description="The Hessian in the archive tail of the Gaussian calculation. link 9999",
    )
    JOB_TIME = MolOPPattern(
        content_pattern=r"(?r)^\s*(Job cpu time|Elapsed time):\s*(\d+)\s*days\s*(\d+)\s*hours\s*(\d+)\s*minutes\s*(\d+\.\d*)\s*seconds",
        content_repeat=0,
        description="The job cpu or elapsed time of the Gaussian calculation. link 9999",
    )
    TERMINATION_STATUS = MolOPPattern(
        content_pattern=r"(?r)^\s*(Normal|Error) termination (.+)\s*(Mon|Tue|Wed|Thu|Fri|Sat|Sun)\s*(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\s*(\d+)\s*(\d+)\s*:\s*(\d+)\s*:\s*(\d+)\s*(\d+)",
        description="The termination status of the Gaussian calculation. link 9999",
    )
    PROCEDURE_TIME = MolOPPattern(
        content_pattern=r"MaxMem=\s*(\d+)\s*cpu:\s*(\d+.\d+)\s*elap:\s*(\d+.\d+)",
        content_repeat=0,
        description="The procedure time of the end of link. link any",
    )


class G16FchkPatterns:
    """
    G16 fchk file patterns
    """


def options_parser(
    link0: str,
) -> Dict[str, str]:
    """
    Parse link0 from G16 log file.

    Args:
        link0 (str): The link0 section from the G16 log file.

    Returns:
        - link0 (Dict): A dictionary containing parameter-value pairs extracted from the link0 section.
    """
    return {
        f"{para.split('=')[0]}": f"{para.split('=')[1]}"
        for para in re.compile(r"(%[a-zA-Z0-9]+)=([^\s]+)").findall(link0)
    }


def parameter_comment_parser(
    route: str,
) -> Tuple[Dict[str, Any], str | None]:
    """
    Parse parameter comment from G16 log file.

    Args:
        parameter_comment (str): The parameter comment from the G16 log file.

    Returns:
        route_params (Dict):
            A dictionary containing parameter-value pairs extracted from the route section of the parameter comment.

        dieze_tag (str):
            The dieze tag extracted from the route section of the parameter comment.
    """
    _route = (
        route.replace(" =", "=")
        .replace("= ", "=")
        .replace(" ,", ",")
        .replace(", ", ",")
        .lower()
    )
    while split_match := re.search(r"/[a-z]", _route):
        _route = _route[: split_match.start()] + " " + _route[split_match.start() + 1 :]
    multi_params_patt = re.compile(
        r"([a-zA-Z\d\+\-]+)([\s=]+)*\(([a-zA-Z\d,=\-\(\)/]+)\)"
    )
    route_params: Dict[str, Any] = {}
    dieze_tag = None

    if _route:
        for tok in _route.split():
            if tok.upper() in ["#", "#N", "#P", "#T"]:
                # does not store # in route to avoid error in input
                dieze_tag = "#N" if tok == "#" else tok
                continue
            else:
                if m := re.search(multi_params_patt, tok.strip("#")):
                    pars = {}
                    for par in m.group(3).split(","):
                        p = par.split("=")
                        pars[p[0]] = None if len(p) == 1 else p[1].lower()
                    route_params[m.group(1).lower()] = pars
                else:
                    d = tok.strip("#").split("=")
                    route_params[d[0]] = None if len(d) == 1 else d[1].lower()

    return route_params, dieze_tag
