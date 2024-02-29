import time
from itertools import chain
from typing import Literal

import numpy as np

from molop.io.bases.molblock_base import QMBaseBlockParser
from molop.logger.logger import logger
from molop.unit import atom_ureg
from molop.utils import parameter_comment_parser
from molop.utils import patterns


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
        matches = patterns["coords start"].search(self._block)
        if matches:
            block = self._block[matches.start() :]
        else:
            block = self._block
        coords = patterns["coords"].findall(
            block,
        )[: self.__n_atom]
        temp_coords = []
        for atom_num, x, y, z in coords:
            self._atoms.append(int(atom_num))
            temp_coords.append((float(x), float(y), float(z)))
        self._coords = np.array(temp_coords) * atom_ureg.angstrom

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
        # self._parse_hessian()
        # self._parse_nbo()

    def _parse_state(self):
        if "opt" in self._lower_route_params:
            matches = patterns["opt stat"].findall(self._block)
            if matches:
                for key, val in matches:
                    self._state[key] = val == "YES"
        matches = patterns["termination"].findall(self._block)
        if matches:
            self._state["termination"] = matches[0]
        matches = patterns["failure reason"].findall(self._block)
        if matches:
            self._state["failure reason"] = matches[0]

    def _parse_energy(self):
        if self.energy is None:
            energy = patterns["energy"].findall(self._block)
            if not energy:
                energy = patterns["energy alter"].findall(self._block)
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
        if patterns["mulliken start"].search(self._block):
            charges = patterns["mulliken match"].findall(
                self._block[
                    patterns["mulliken start"]
                    .search(self._block)
                    .start() : patterns["mulliken end"]
                    .search(self._block)
                    .end()
                ]
            )
            for charge in charges:
                self._partial_charges.append(float(charge))

    def _parse_spin_density(self):
        if patterns["spin density start"].search(self._block):
            spin_densities = patterns["spin density match"].findall(
                self._block[
                    patterns["spin density start"]
                    .search(self._block)
                    .start() : patterns["spin density end"]
                    .search(self._block)
                    .end()
                ]
            )
            for spin_density in spin_densities:
                self._spin_densities.append(float(spin_density))

    def _parse_gradient(self):
        if patterns["gradients start"].search(self._block):
            gradients = patterns["gradients match"].findall(
                self._block[
                    patterns["gradients start"]
                    .search(self._block)
                    .start() : patterns["gradients end"]
                    .search(self._block)
                    .end()
                ]
            )
            for fx, fy, fz in gradients:
                self._gradients.append(
                    (
                        float(fx) * atom_ureg.hartree / atom_ureg.bohr,
                        float(fy) * atom_ureg.hartree / atom_ureg.bohr,
                        float(fz) * atom_ureg.hartree / atom_ureg.bohr,
                    )
                )

    def _parse_spins(self):
        s = patterns["spins"].findall(self._block)
        if s:
            self._spin_multiplicity = float(s[-1][1])
            self._spin_eigenvalue = float(s[-1][0])

    def _parse_orbitals(self):
        matches = patterns["orbital"].findall(self._block)
        for orbital_type, occ_stat, energies in matches:
            if orbital_type == "Alpha":
                self._alpha_FMO_orbits.extend(
                    (
                        round(e, 6) * atom_ureg.hartree / atom_ureg.particle
                        for e in map(
                            float,
                            [energies[j : j + 10] for j in range(0, len(energies), 10)],
                        )
                    )
                )
                if occ_stat == "occ.":
                    homo_idx = len(self._alpha_FMO_orbits) - 1
            elif orbital_type == "Beta":
                self._beta_FMO_orbits.extend(
                    (
                        round(e, 6) * atom_ureg.hartree / atom_ureg.particle
                        for e in map(
                            float,
                            [energies[j : j + 10] for j in range(0, len(energies), 10)],
                        )
                    )
                )
                if occ_stat == "occ.":
                    homo_idx = len(self._beta_FMO_orbits) - 1
        if self._alpha_FMO_orbits:
            self._alpha_energy["homo"] = self._alpha_FMO_orbits[homo_idx]
            self._alpha_energy["lumo"] = self._alpha_FMO_orbits[homo_idx + 1]
            self._alpha_energy["gap"] = (
                round((self._alpha_energy["lumo"] - self._alpha_energy["homo"]).m, 6)
                * atom_ureg.hartree
                / atom_ureg.particle
            )
        if self._beta_FMO_orbits:
            self._beta_energy["homo"] = self._beta_FMO_orbits[homo_idx]
            self._beta_energy["lumo"] = self._beta_FMO_orbits[homo_idx + 1]
            self._beta_energy["gap"] = (
                round((self._beta_energy["lumo"] - self._beta_energy["homo"]).m, 6)
                * atom_ureg.hartree
                / atom_ureg.particle
            )

    def _parse_frequencies(self):
        if "freq" in self._lower_route_params:
            matches = patterns["freq start"].search(self._block)
            if matches:
                block = self._block[matches.end() :]
                freqs = patterns["freq"].findall(block)
                red_masses = patterns["freq red. masses"].findall(block)
                frc_consts = patterns["freq frc consts"].findall(block)
                ir_intens = patterns["freq IR Inten"].findall(block)
                freq_modes = patterns["freq mode"].findall(block)
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
                            "normal coordinates": np.array(freq_mode).astype(np.float16)
                            * atom_ureg.angstrom,
                        }
                    )

    def _parse_sum_energy(self):
        mappings = {
            "zero-point Energies": "zero-point gas",
            "thermal Energies": "E gas",
            "thermal Enthalpies": "H gas",
            "thermal Free Energies": "G gas",
            ("Zero-point", ""): "zero-point correction",
            ("Thermal", " to Energy"): "TCE",
            ("Thermal", " to Enthalpy"): "TCH",
            ("Thermal", " to Gibbs Free Energy"): "TCG",
        }
        matches = patterns["sum energies"].findall(self._block)
        if matches:
            for tag, val in matches:
                self._sum_energy[mappings[tag]] = round(
                    float(val) * atom_ureg.hartree / atom_ureg.particle, 6
                )
        matches = patterns["corrections"].findall(self._block)
        if matches:
            for tag1, tag2, val in matches:
                self._sum_energy[mappings[(tag1, tag2)]] = round(
                    float(val) * atom_ureg.hartree / atom_ureg.particle, 6
                )

    def _parse_hessian(self):
        lines = self._block.splitlines()
        hess_lines = []
        append_line = False
        for line in reversed(lines):
            if r"\\@" in line or line.startswith(" @") or line.startswith(r" \@"):
                append_line = True
            if append_line:
                hess_lines.append(line.strip("\n").strip(" "))
            if "NImag" in line:
                append_line = "end"
                break
        if append_line == "end":
            hess_str = "".join(hess_lines[::-1]).split(r"\\")[-3]
            hess_val = [float(val) for val in hess_str.split(",")]
            n = self.__n_atom * 3
            if len(hess_val) != n * (n + 1) // 2:
                logger.warning(
                    f"The number of elements {len(hess_val)} in the Hessian matrix is not consistent with the number of atoms {self.__n_atom} in {self._file_path}"
                )
            else:
                self._hessian = hess_val

    def _parse_nbo(self):
        lines = self._block.splitlines()
        for idx, line in enumerate(lines):
            if "Summary of Natural Population Analysis" in line:
                for i in range(idx + 6, idx + 6 + self.__n_atom):
                    _, _, natural_charge, core, valence, rydberg, _ = lines[i].split()
                    self._nbo_analysis.append(
                        {
                            "natural_charge": float(natural_charge),
                            "core": float(core),
                            "valence": float(valence),
                            "rydberg": float(rydberg),
                        }
                    )
                break
