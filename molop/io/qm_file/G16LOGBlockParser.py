import time
import math
from itertools import chain
from typing import Literal

import numpy as np

from molop.io.bases.molblock_base import QMBaseBlockParser
from molop.logger.logger import logger
from molop.unit import atom_ureg
from molop.utils import g16logpatterns, parameter_comment_parser


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
        version=None,
        parameter_comment=None,
        only_extract_structure=False,
    ):
        super().__init__(block, only_extract_structure)
        self._file_path = file_path
        self._charge = charge
        self._multiplicity = multiplicity
        self.__n_atom = n_atom
        self._version = version
        self._parameter_comment = parameter_comment
        (
            self._link0,
            self._route_params,
            self._dieze_tag,
            self._functional,
            self._basis_set,
        ) = parameter_comment_parser(self._parameter_comment)
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

    @property
    def functional(self) -> str:
        return self._functional

    @property
    def basis_set(self) -> str:
        return self._basis_set

    @property
    def _lower_route_params(self) -> dict:
        return {k.lower(): v for k, v in self._route_params.items()}

    def _parse_coords(self):
        matches = g16logpatterns["coords start"].search(self._block)
        if matches:
            block = self._block[matches.start() :]
        else:
            block = self._block
        coords = g16logpatterns["coords"].findall(
            block,
        )[: self.__n_atom]
        temp_coords = []
        for atom_num, x, y, z in coords:
            self._atoms.append(int(atom_num))
            temp_coords.append((float(x), float(y), float(z)))
        self._coords = np.array(temp_coords, np.float32) * atom_ureg.angstrom

    def _parse(self):
        self._parse_state()
        self._parse_energy()
        self._parse_mulliken_charges()
        self._parse_spin_density()
        self._parse_gradient()
        self._parse_spins()
        self._parse_orbitals()
        self._parse_sum_energy()
        self._parse_frequencies()
        self._wiberg_bond_order = self._parse_bond_order("wiberg")
        self._atom_atom_overlap_bond_order = self._parse_bond_order("atom_atom_overlap")
        self._mo_bond_order = self._parse_bond_order("mo")
        self._parse_nbo_charges()
        self._parse_nbo_bond_order()
        self._parse_dipole()

    def _parse_state(self):
        if "opt" in self._lower_route_params:
            matches = g16logpatterns["opt stat"].findall(self._block)
            if matches:
                for key, val in matches:
                    self._state[key] = val == "YES"
        matches = g16logpatterns["termination"].findall(self._block)
        if matches:
            self._state["termination"] = matches[0]
        matches = g16logpatterns["failure reason"].findall(self._block)
        if matches:
            self._state["failure reason"] = matches[0]

    def _parse_energy(self):
        if self.energy is None:
            energy = g16logpatterns["energy"].findall(self._block)
            if not energy:
                energy = g16logpatterns["energy alter"].findall(self._block)
                if energy:
                    energy = ["".join(energy[0].split("\n "))]
            if energy:
                self._energy = (
                    round(float(energy[0]), 6) * atom_ureg.hartree / atom_ureg.particle
                )
                self._state["SCF Done"] = True
            else:
                self._state["SCF Done"] = False

    def _parse_mulliken_charges(self):
        if g16logpatterns["mulliken start"].search(self._block):
            charges = g16logpatterns["mulliken match"].findall(
                self._block[
                    g16logpatterns["mulliken start"]
                    .search(self._block)
                    .start() : g16logpatterns["mulliken end"]
                    .search(self._block)
                    .end()
                ]
            )
            for charge in charges:
                self._partial_charges.append(float(charge))

    def _parse_spin_density(self):
        if g16logpatterns["spin density start"].search(self._block):
            spin_densities = g16logpatterns["spin density match"].findall(
                self._block[
                    g16logpatterns["spin density start"]
                    .search(self._block)
                    .start() : g16logpatterns["spin density end"]
                    .search(self._block)
                    .end()
                ]
            )
            for spin_density in spin_densities:
                self._spin_densities.append(float(spin_density))

    def _parse_gradient(self):
        if g16logpatterns["gradients start"].search(self._block):
            gradients = g16logpatterns["gradients match"].findall(
                self._block[
                    g16logpatterns["gradients start"]
                    .search(self._block)
                    .start() : g16logpatterns["gradients end"]
                    .search(self._block)
                    .end()
                ]
            )
            temp_gradients = []
            for fx, fy, fz in gradients:
                temp_gradients.append((float(fx), float(fy), float(fz)))
            self._gradients = (
                np.array(temp_gradients, np.float32)
                * atom_ureg.hartree
                / atom_ureg.bohr
            )

    def _parse_spins(self):
        s = g16logpatterns["spins"].findall(self._block)
        if s:
            self._spin_multiplicity = float(s[-1][1])
            self._spin_eigenvalue = float(s[-1][0])

    def _parse_orbitals(self):
        matches = g16logpatterns["orbital"].findall(self._block)
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
                            [energies[j : j + 10] for j in range(0, len(energies), 10)],
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
                            [energies[j : j + 10] for j in range(0, len(energies), 10)],
                        )
                    )
                )
                if occ_stat == "occ.":
                    homo_idx = len(temp_beta_orbitals) - 1
        if temp_alpha_orbitals and homo_idx:
            self._alpha_FMO_orbits = (
                np.array(temp_alpha_orbitals, dtype=np.float32)
                * atom_ureg.hartree
                / atom_ureg.particle
            )
            self._alpha_energy["homo"] = self._alpha_FMO_orbits[homo_idx]
            self._alpha_energy["lumo"] = self._alpha_FMO_orbits[homo_idx + 1]
            self._alpha_energy["gap"] = (
                round((self._alpha_energy["lumo"] - self._alpha_energy["homo"]).m, 6)
                * atom_ureg.hartree
                / atom_ureg.particle
            )
        if temp_beta_orbitals and homo_idx:
            self._beta_FMO_orbits = (
                np.array(temp_beta_orbitals, dtype=np.float32)
                * atom_ureg.hartree
                / atom_ureg.particle
            )
            self._beta_energy["homo"] = self._beta_FMO_orbits[homo_idx]
            self._beta_energy["lumo"] = self._beta_FMO_orbits[homo_idx + 1]
            self._beta_energy["gap"] = (
                round((self._beta_energy["lumo"] - self._beta_energy["homo"]).m, 6)
                * atom_ureg.hartree
                / atom_ureg.particle
            )

    def _parse_frequencies(self):
        if "freq" in self._lower_route_params:
            matches = g16logpatterns["freq start"].search(self._block)
            if matches:
                block = self._block[matches.end() :]
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
                    self._frequencies.append(
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

    def _parse_sum_energy(self):
        mappings = {
            "zero-point Energies": "zero-point gas",
            "thermal Energies": "E sum",
            "thermal Enthalpies": "H sum",
            "thermal Free Energies": "G sum",
            ("Zero-point", ""): "zero-point correction",
            ("Thermal", " to Energy"): "TCE",
            ("Thermal", " to Enthalpy"): "TCH",
            ("Thermal", " to Gibbs Free Energy"): "TCG",
        }
        matches = g16logpatterns["sum energies"].findall(self._block)
        if matches:
            for tag, val in matches:
                self._sum_energy[mappings[tag]] = round(
                    float(val) * atom_ureg.hartree / atom_ureg.particle, 6
                )
        matches = g16logpatterns["corrections"].findall(self._block)
        if matches:
            for tag1, tag2, val in matches:
                self._sum_energy[mappings[(tag1, tag2)]] = round(
                    float(val) * atom_ureg.hartree / atom_ureg.particle, 6
                )

    def _parse_bond_order(self, _type: Literal["wiberg", "atom_atom_overlap", "mo"]):
        start_matches = g16logpatterns[f"{_type}_start"].search(self._block)
        end_matches = g16logpatterns[f"{_type}_end"].search(self._block)
        if start_matches and end_matches:
            temp_block = self._block[start_matches.start() : end_matches.end()]
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
        return None

    def _parse_dipole(self):
        matches = g16logpatterns["dipole"].search(self._block)
        if matches:
            dipole = np.array(list(map(float, matches.groups())))
            self._dipole = dipole * atom_ureg.debye

    def _parse_nbo_charges(self):
        start_matches = g16logpatterns[f"nbo charge start"].search(self._block)
        end_matches = g16logpatterns[f"nbo charge end"].search(self._block)
        if start_matches and end_matches:
            temp_block = self._block[start_matches.start() : end_matches.end()]
            charges = list(
                map(float, g16logpatterns["nbo charge match"].findall(temp_block))
            )
            self._nbo_charges = charges

    def _parse_nbo_bond_order(self):
        start_matches = g16logpatterns[f"nbo summary start"].search(self._block)
        end_matches = g16logpatterns[f"nbo summary end"].search(self._block)
        if start_matches and end_matches:
            temp_block = self._block[start_matches.start() : end_matches.end()]
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
            self._nbo_bond_order = [
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

    def is_error(self) -> bool:
        """
        Check if the current frame is an error frame
        """
        if "termination" in self.state and self.state["termination"] == "Error":
            return True
        if "SCF Done" in self.state and self.state["SCF Done"] == False:
            return True
        return False
