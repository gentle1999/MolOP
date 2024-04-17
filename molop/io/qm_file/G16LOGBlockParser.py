import math
from itertools import chain
from typing import Literal

import numpy as np

from molop.io.bases.molblock_base import QMBaseBlockParser
from molop.logger.logger import logger
from molop.unit import atom_ureg
from molop.utils.g16patterns import (
    g16logpatterns,
    get_solvent,
    get_solvent_model,
    link0_parser,
    parameter_comment_parser,
)


class G16LOGBlockParser(QMBaseBlockParser):
    """
    Parser for G16 log Blocks.
    """

    _block_type = "G16 LOG"

    def __init__(
        self,
        block: str,
        charge=0,
        multiplicity=1,
        n_atom=1,
        file_path="",
        version: str = None,
        parameter_comment: str = None,
        functional: str = None,
        basis: str = None,
        only_extract_structure=False,
    ):
        super().__init__(block, only_extract_structure)
        self.qm_software = "Gaussian"
        self.__block = block
        self._file_path = file_path
        self._charge = charge
        self._multiplicity = multiplicity
        self.__n_atom = n_atom
        self.version = version
        self.parameter_comment = parameter_comment
        link, route = self.parameter_comment.split("#")
        link = link.replace("\n", " ")
        route = "#" + route.replace("\n", " ")
        self._link0 = link0_parser(link)

        (
            self._route_params,
            self._dieze_tag,
        ) = parameter_comment_parser(route)
        self.functional = functional
        self.basis = basis

        self._parse_coords()

        if not self._only_extract_structure:
            self._parse()

    @property
    def link0(self) -> dict:
        return self._link0

    @property
    def route_params(self) -> dict:
        return self._route_params

    @property
    def dieze_tag(self) -> Literal["#N", "#P", "#T"]:
        return self._dieze_tag

    def _parse_coords(self):
        corrds_end = g16logpatterns["coords_end"].search(self.__block)
        if corrds_end:
            block = self.__block[: corrds_end.end()]
            coords = g16logpatterns["coords"].findall(
                block,
            )[-self.__n_atom :]
            temp_coords = []
            for atom_num, x, y, z in coords:
                self._atoms.append(int(atom_num))
                temp_coords.append((float(x), float(y), float(z)))
            self._coords = np.array(temp_coords, np.float32) * atom_ureg.angstrom

    def _parse(self):
        self._parse_rotation_consts()
        self._parse_functional_basis()
        self._parse_solvent()
        self._parse_energy()
        self._parse_spins()
        self._parse_orbitals()
        self._parse_mulliken_charges()
        self._parse_mulliken_spin_density()
        self._parse_apt_charges()
        self._parse_lowdin_charges()
        self._parse_dipole()
        self._parse_quadrupole()
        self._parse_octapole()
        self._parse_hexadecapole()
        self._parse_hirshfeld_spin_charges()
        self.wiberg_bond_order = self._parse_bond_order("wiberg")
        self.atom_atom_overlap_bond_order = self._parse_bond_order("atom_atom_overlap")
        self.mo_bond_order = self._parse_bond_order("mo")
        self._parse_nbo_charges()
        self._parse_nbo_bond_order()
        self._parse_frequencies()
        self._parse_sum_energy()
        self._parse_gradient()
        self._parse_state()
        self._parse_tail()

    def _parse_rotation_consts(self):
        rotation_consts = g16logpatterns["rotation_consts"].search(
            self.__block,
        )
        if rotation_consts:
            self._rotation_consts = (
                np.array([float(x) for x in rotation_consts.groups()]) * atom_ureg.GHz
            )
        else:
            self._rotation_consts = np.array([]) * atom_ureg.GHz

    def _parse_functional_basis(self):
        functional_match = g16logpatterns["functional"].search(self.__block)
        if functional_match:
            self.functional = functional_match.group(1).lower()
        basis_match = g16logpatterns["basis"].search(self.__block)
        if basis_match:
            self.basis = basis_match.group(1)
        if g16logpatterns["Pseudopotential"].search(self.__block):
            self.basis = "pseudopotential"

    def _parse_solvent(self):
        solvent_start = g16logpatterns["solvent_start"].search(self.__block)
        solvent_end = g16logpatterns["solvent_end"].search(self.__block)
        if solvent_start and solvent_end:
            solvent_block = self.__block[solvent_start.start() : solvent_end.end()]
            self.__block = self.__block[solvent_end.start() :]
            solvent_model = (
                g16logpatterns["solvent_model"].search(solvent_block).group(1)
            )
            solvent_atom_radii = (
                g16logpatterns["solvent_atomic_radii"].search(solvent_block).group(1)
            )
            if solvent_model == "PCM":
                self.solvent_model = "IEFPCM"
            elif solvent_model == "C-PCM":
                self.solvent_model = "C-PCM"
            elif solvent_atom_radii == "SMD-Coulomb":
                self.solvent_model = "SMD-IEFPCM"
            else:
                self.solvent_model = solvent_model
            self.solvent = g16logpatterns["solvent_type"].search(solvent_block).group(1)
            self._eps = float(
                g16logpatterns["solvent_eps"].search(solvent_block).group(1)
            )
            self._eps_inf = float(
                g16logpatterns["solvent_eps_inf"].search(solvent_block).group(2)
            )
        else:
            self.solvent_model = ""
            self.solvent = ""
            self.solvent_eps = None
            self.solvent_eps_inf = None

    def _parse_energy(self):
        scf_energy = g16logpatterns["scf_energy"].search(self.__block)
        if scf_energy:
            self.__block = self.__block[scf_energy.start() :]
            self.scf_energy = (
                round(float(scf_energy.group(1)), 6)
                * atom_ureg.hartree
                / atom_ureg.particle
            )
            self.status["SCF Done"] = True
        mp2_4_energy = g16logpatterns["mp2-4"].findall(self.__block)
        for mp_energy in mp2_4_energy:
            setattr(
                self,
                f"{mp_energy[0].lower()}_energy",
                round(float(mp_energy[1].replace("D", "E")), 6)
                * atom_ureg.hartree
                / atom_ureg.particle,
            )
            self.functional = mp_energy[0].lower()
        ccsd_energy = g16logpatterns["ccsd"].search(self.__block)
        if ccsd_energy:
            self.ccsd_energy = (
                round(float(ccsd_energy.group(1)), 6)
                * atom_ureg.hartree
                / atom_ureg.particle
            )
            self.functional = "ccsd"
        ccsd_energy = g16logpatterns["ccsd(t)"].search(self.__block)
        if ccsd_energy:
            self.ccsd_energy = (
                round(float(ccsd_energy.group(1).replace("D", "E")), 6)
                * atom_ureg.hartree
                / atom_ureg.particle
            )
            self.functional = "ccsd(t)"

    def _parse_spins(self):
        s = g16logpatterns["spins"].search(self.__block)
        if s:
            self.spin_multiplicity = float(s.group(1))
            self.spin_eigenvalue = float(s.group(2))
            self._multiplicity = int(round(2 * self.spin_eigenvalue + 1, 0))

    def _parse_orbitals(self):
        orbital_start = g16logpatterns["orbital_start"].search(self.__block)
        orbital_end = g16logpatterns["orbital_end"].search(self.__block)
        if orbital_start and orbital_end:
            block = self.__block[orbital_start.start() : orbital_end.end()]
            self.__block = self.__block[orbital_end.start() :]
            matches = g16logpatterns["orbital"].findall(block)
            temp_alpha_orbitals = []
            temp_beta_orbitals = []
            homo_idx = None
            for orbital_type, occ_stat, energies in matches:
                if orbital_type == "Alpha":
                    temp_alpha_orbitals.extend(
                        (
                            round(e, 6)
                            for e in map(
                                float,
                                [
                                    energies[j : j + 10]
                                    for j in range(0, len(energies), 10)
                                ],
                            )
                        )
                    )
                    if occ_stat == "occ.":
                        homo_idx = len(temp_alpha_orbitals) - 1
                elif orbital_type == "Beta":
                    temp_beta_orbitals.extend(
                        (
                            round(e, 6)
                            for e in map(
                                float,
                                [
                                    energies[j : j + 10]
                                    for j in range(0, len(energies), 10)
                                ],
                            )
                        )
                    )
                    if occ_stat == "occ.":
                        homo_idx = len(temp_beta_orbitals) - 1
            if temp_alpha_orbitals and homo_idx:
                self.alpha_FMO_orbits = (
                    np.array(temp_alpha_orbitals, dtype=np.float32)
                    * atom_ureg.hartree
                    / atom_ureg.particle
                )
                self.alpha_energy["homo"] = self.alpha_FMO_orbits[homo_idx]
                self.alpha_energy["lumo"] = self.alpha_FMO_orbits[homo_idx + 1]
                self.alpha_energy["gap"] = (
                    round((self.alpha_energy["lumo"] - self.alpha_energy["homo"]).m, 6)
                    * atom_ureg.hartree
                    / atom_ureg.particle
                )
            if temp_beta_orbitals and homo_idx:
                self.beta_FMO_orbits = (
                    np.array(temp_beta_orbitals, dtype=np.float32)
                    * atom_ureg.hartree
                    / atom_ureg.particle
                )
                self.beta_energy["homo"] = self.beta_FMO_orbits[homo_idx]
                self.beta_energy["lumo"] = self.beta_FMO_orbits[homo_idx + 1]
                self.beta_energy["gap"] = (
                    round((self.beta_energy["lumo"] - self.beta_energy["homo"]).m, 6)
                    * atom_ureg.hartree
                    / atom_ureg.particle
                )

    def _parse_mulliken_charges(self):
        mulliken_start = g16logpatterns["mulliken start"].search(self.__block)
        mulliken_end = g16logpatterns["mulliken end"].search(self.__block)
        if mulliken_start and mulliken_end:
            charges = g16logpatterns["mulliken match"].findall(
                self.__block[mulliken_start.start() : mulliken_end.end()]
            )
            self.mulliken_charges = [float(charge) for charge in charges]

    def _parse_mulliken_spin_density(self):
        spin_densities_start = g16logpatterns["mulliken spin density start"].search(
            self.__block
        )
        spin_densities_end = g16logpatterns["mulliken spin density end"].search(
            self.__block
        )
        if spin_densities_start and spin_densities_end:
            spin_densities = g16logpatterns["mulliken spin density match"].findall(
                self.__block[spin_densities_start.start() : spin_densities_end.end()]
            )
            self.__block = self.__block[spin_densities_end.end() :]
            self.mulliken_spin_densities = [
                float(spin_density) for spin_density in spin_densities
            ]

    def _parse_apt_charges(self):
        apt_charges_start = g16logpatterns["apt_start"].search(self.__block)
        apt_charges_end = g16logpatterns["apt_end"].search(self.__block)
        if apt_charges_start and apt_charges_end:
            apt_charges = g16logpatterns["apt_match"].findall(
                self.__block[apt_charges_start.start() : apt_charges_end.end()]
            )
            self.apt_charges = [float(charge) for charge in apt_charges]

    def _parse_lowdin_charges(self):
        lowdin_charges_start = g16logpatterns["lowdin_start"].search(self.__block)
        lowdin_charges_end = g16logpatterns["lowdin_end"].search(self.__block)
        if lowdin_charges_start and lowdin_charges_end:
            lowdin_charges = g16logpatterns["lowdin_match"].findall(
                self.__block[lowdin_charges_start.start() : lowdin_charges_end.end()]
            )
            self.lowdin_charges = [float(charge) for charge in lowdin_charges]

    def _parse_dipole(self):
        start = g16logpatterns["dipole_start"].search(self.__block)
        end = g16logpatterns["quadrupole_start"].search(self.__block)
        if start and end:
            dipole = np.array(
                [
                    float(l)
                    for l in g16logpatterns["dipole"].findall(
                        self.__block[start.start() : end.end()]
                    )
                ]
            )
            self.dipole = dipole * atom_ureg.debye

    def _parse_quadrupole(self):
        start = g16logpatterns["quadrupole_start"].search(self.__block)
        end = g16logpatterns["octapole_start"].search(self.__block)
        if start and end:
            quadrupole = np.array(
                [
                    float(l)
                    for l in g16logpatterns["quadrupole"].findall(
                        self.__block[start.start() : end.end()]
                    )
                ]
            )
            self.quadrupole = quadrupole * atom_ureg.debye * atom_ureg.angstrom

    def _parse_octapole(self):
        start = g16logpatterns["octapole_start"].search(self.__block)
        end = g16logpatterns["hexadecapole_start"].search(self.__block)
        if start and end:
            octapole = np.array(
                [
                    float(l)
                    for l in g16logpatterns["octapole"].findall(
                        self.__block[start.start() : end.end()]
                    )
                ]
            )
            self.octapole = (
                octapole * atom_ureg.debye * atom_ureg.angstrom * atom_ureg.angstrom
            )

    def _parse_hexadecapole(self):
        start = g16logpatterns["hexadecapole_start"].search(self.__block)
        end = g16logpatterns["hexadecapole_end"].search(self.__block)
        if start and end:
            hexadecapole = np.array(
                [
                    float(l)
                    for l in g16logpatterns["hexadecapole"].findall(
                        self.__block[start.start() : end.end()]
                    )
                ]
            )
            self.hexadecapole = (
                hexadecapole
                * atom_ureg.debye
                * atom_ureg.angstrom
                * atom_ureg.angstrom
                * atom_ureg.angstrom
            )

    def _parse_hirshfeld_spin_charges(self):
        start = g16logpatterns["hirshfeld_start"].search(self.__block)
        end = g16logpatterns["hirshfeld_end"].search(self.__block)
        if start and end:
            hirshfeld = g16logpatterns["hirshfeld"].findall(
                self.__block[start.start() : end.end()]
            )
            self.hirshfeld_charges = [float(charge) for charge, _, _ in hirshfeld]
            self.hirshfeld_spins = [float(spin) for _, spin, _ in hirshfeld]
            self.hirshfeld_q_cm5 = [float(q_cm5) for _, _, q_cm5 in hirshfeld]

    def _parse_bond_order(self, _type: Literal["wiberg", "atom_atom_overlap", "mo"]):
        start_matches = g16logpatterns[f"{_type}_start"].search(self.__block)
        end_matches = g16logpatterns[f"{_type}_end"].search(self.__block)
        if start_matches and end_matches:
            temp_block = self.__block[start_matches.start() : end_matches.end()]
            digits = list(map(float, g16logpatterns["digit"].findall(temp_block)))
            blks = []
            for idx in range(math.ceil(len(digits) / (self.__n_atom * 9))):
                block_length = (
                    min(self.__n_atom * 9, len(digits) - idx * self.__n_atom * 9)
                    // self.__n_atom
                )
                blks.append(
                    np.array(
                        digits[
                            idx * self.__n_atom * 9 : idx * self.__n_atom * 9
                            + min(
                                self.__n_atom * 9, len(digits) - idx * self.__n_atom * 9
                            )
                        ]
                    ).reshape(-1, block_length)
                )
            return np.concatenate(blks, axis=1)
        return np.array([[]])

    def _parse_nbo_charges(self):
        start_matches = g16logpatterns[f"nbo charge start"].search(self.__block)
        end_matches = g16logpatterns[f"nbo charge end"].search(self.__block)
        if start_matches and end_matches:
            temp_block = self.__block[start_matches.start() : end_matches.end()]
            charges = list(
                map(float, g16logpatterns["nbo charge match"].findall(temp_block))
            )
            self.nbo_charges = charges

    def _parse_nbo_bond_order(self):
        start_matches = g16logpatterns[f"nbo summary start"].search(self.__block)
        end_matches = g16logpatterns[f"nbo summary end"].search(self.__block)
        if start_matches and end_matches:
            temp_block = self.__block[start_matches.start() : end_matches.end()]
            temp_dict = {}
            for idx1, idx2, occ, energy in g16logpatterns["nbo summary match"].findall(
                temp_block
            ):
                _idx1, _idx2, _occ, _energy = (
                    int(idx1),
                    int(idx2),
                    float(occ),
                    float(energy),
                )
                if (_idx1, _idx2) in temp_dict:
                    temp_dict[(_idx1, _idx2)]["energy"] += _energy * _occ
                    temp_dict[(_idx1, _idx2)]["bond_order"] += 1
                else:
                    temp_dict[(_idx1, _idx2)] = {
                        "energy": _energy * _occ,
                        "bond_order": 1,
                    }
            self.nbo_bond_order = [
                (
                    idx1 - 1,
                    idx2 - 1,
                    temp_dict[(idx1, idx2)]["bond_order"],
                    temp_dict[(idx1, idx2)]["energy"]
                    * atom_ureg.hartree
                    / atom_ureg.particle,
                )
                for idx1, idx2 in temp_dict
            ]

    def _parse_frequencies(self):
        start = g16logpatterns["freq start"].search(self.__block)
        end = g16logpatterns["freq end"].search(self.__block)
        if start and end:
            block = self.__block[start.end() :]
            self.__block = self.__block[end.start() :]
            freqs = g16logpatterns["freq"].findall(block)
            red_masses = g16logpatterns["freq red. masses"].findall(block)
            frc_consts = g16logpatterns["freq frc consts"].findall(block)
            ir_intens = g16logpatterns["freq IR Inten"].findall(block)
            freq_modes = g16logpatterns["freq mode"].findall(block)
            freq_modes_reindex = []
            for i in range(0, len(freq_modes), 3 * self.__n_atom):
                freq_modes_reindex.append(freq_modes[i : i + 3 * self.__n_atom : 3])
                freq_modes_reindex.append(
                    freq_modes[i + 1 : i + 1 + 3 * self.__n_atom : 3]
                )
                freq_modes_reindex.append(
                    freq_modes[i + 2 : i + 2 + 3 * self.__n_atom : 3]
                )

            for freq, red_mass, frc_const, ir_inten, freq_mode in zip(
                chain.from_iterable(freqs),
                chain.from_iterable(red_masses),
                chain.from_iterable(frc_consts),
                chain.from_iterable(ir_intens),
                freq_modes_reindex,
            ):
                self.frequencies.append(
                    {
                        "is imaginary": float(freq) < 0,
                        "freq": float(freq) * atom_ureg.cm_1,
                        "reduced masses": float(red_mass) * atom_ureg.amu,
                        "force constants": (
                            float(frc_const) * atom_ureg.mdyne / atom_ureg.angstrom
                        ),
                        "IR intensities": (
                            float(ir_inten) * atom_ureg.kmol / atom_ureg.mol
                        ),
                        "normal coordinates": np.array(freq_mode, dtype=np.float32)
                        * atom_ureg.angstrom,
                    }
                )

        temperature_match = g16logpatterns["Temperature"].search(self.__block)
        if temperature_match:
            self.temperature = float(temperature_match.group(1))

    def _parse_state(self):
        if "opt" in self._route_params:
            matches = g16logpatterns["opt stat"].findall(self.__block)
            if matches:
                for key, val in matches:
                    self.status[key] = val == "YES"
        matches = g16logpatterns["termination"].findall(self.__block)
        if matches:
            self.status["termination"] = matches[0]
        matches = g16logpatterns["failure reason"].findall(self.__block)
        if matches:
            self.status["failure reason"] = matches[0]

    def _parse_sum_energy(self):
        mappings = {
            "zero-point Energies": "zero-point sum",
            "thermal Energies": "E sum",
            "thermal Enthalpies": "H sum",
            "thermal Free Energies": "G sum",
            ("Zero-point", ""): "zero-point correction",
            ("Thermal", " to Energy"): "TCE",
            ("Thermal", " to Enthalpy"): "TCH",
            ("Thermal", " to Gibbs Free Energy"): "TCG",
        }
        matches = g16logpatterns["sum energies"].findall(self.__block)
        if matches:
            for tag, val in matches:
                self.sum_energy[mappings[tag]] = round(
                    float(val) * atom_ureg.hartree / atom_ureg.particle, 6
                )
        matches = g16logpatterns["corrections"].findall(self.__block)
        if matches:
            for tag1, tag2, val in matches:
                self.sum_energy[mappings[(tag1, tag2)]] = round(
                    float(val) * atom_ureg.hartree / atom_ureg.particle, 6
                )

    def _parse_gradient(self):
        start = g16logpatterns["gradients start"].search(self.__block)
        end = g16logpatterns["gradients end"].search(self.__block)
        if start and end:
            gradients = g16logpatterns["gradients match"].findall(
                self.__block[start.start() : end.end()]
            )
            temp_gradients = [
                (float(fx), float(fy), float(fz)) for fx, fy, fz in gradients
            ]
            self.gradients = (
                np.array(temp_gradients, np.float32)
                * atom_ureg.hartree
                / atom_ureg.bohr
            )

    def is_error(self) -> bool:
        """
        Check if the current frame is an error frame
        """
        if "termination" in self.status and self.status["termination"] == "Error":
            return True
        if "SCF Done" in self.status and self.status["SCF Done"] == False:
            return True
        return False

    def _parse_tail(self):
        start = g16logpatterns["tail_start"].search(self.__block)
        end = g16logpatterns["tail_end"].search(self.__block)
        if start and end:
            tail = self.__block[start.end() : end.start()]
            tail = tail.replace("\n ", "")
            energies_match = g16logpatterns["tail_match"].findall(tail)
            for e, v in energies_match:
                if "HF" in e:
                    self.scf_energy = float(v) * atom_ureg.hartree / atom_ureg.particle
                if "MP2" in e:
                    self.mp2_energy = float(v) * atom_ureg.hartree / atom_ureg.particle
                if "MP3" in e:
                    self.mp2_energy = float(v) * atom_ureg.hartree / atom_ureg.particle
                if "MP4" in e:
                    self.mp2_energy = float(v) * atom_ureg.hartree / atom_ureg.particle
                if "CCSD" in e:
                    self.ccsd_energy = float(v) * atom_ureg.hartree / atom_ureg.particle
            thermal_match = g16logpatterns["tail_thermal_match"].findall(tail)
            for e, v in thermal_match:
                if "ZeroPoint" in e:
                    self.sum_energy["zero-point correction"] = float(v) * atom_ureg.hartree / atom_ureg.particle
                if "Thermal" in e:
                    self.sum_energy["TCE"] = float(v) * atom_ureg.hartree / atom_ureg.particle
                if "ETot" in e:
                    self.sum_energy["E sum"] = float(v) * atom_ureg.hartree / atom_ureg.particle
                if "HTot" in e:
                    self.sum_energy["H sum"] = float(v) * atom_ureg.hartree / atom_ureg.particle
                if "GTot" in e:
                    self.sum_energy["G sum"] = float(v) * atom_ureg.hartree / atom_ureg.particle
