"""Electronic-state and multi-reference result data classes."""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from typing import Any, ClassVar, overload

import numpy as np
from pint._typing import UnitLike
from pint.facets.plain import PlainQuantity
from pydantic import Field

from molop.io.base_models.Bases import (
    BaseDataClassWithUnit,
    PropertyBundle,
    PropertyColumnValue,
    PropertyScalarValue,
    PropertyTable,
    PropertyTransition,
    Spectrum,
)
from molop.unit import atom_ureg
from molop.utils.types import PintArray3

from .qm_requests import ActiveSpace


class ElectronicConfiguration(BaseDataClassWithUnit):
    label: str | None = Field(default=None, description="Configuration label")
    coefficient: float | None = Field(default=None, description="CI or configuration coefficient")
    weight: float | None = Field(default=None, description="Configuration weight")
    occupation: list[float] = Field(default_factory=list, description="Orbital occupation pattern")
    orbital_indices: list[int] = Field(default_factory=list, description="0-based orbital indices")
    raw: str = Field(default="", description="Raw configuration text")


class ElectronicState(BaseDataClassWithUnit):
    default_units: ClassVar[dict[str, UnitLike]] = {
        "energy": atom_ureg.hartree,
        "excitation_energy": atom_ureg.eV,
        "transition_dipole": atom_ureg.debye,
    }
    set_default_units: ClassVar[bool] = True

    state_index: int | None = Field(default=None, description="0-based electronic state index")
    root: int | None = Field(default=None, description="Program root index")
    label: str | None = Field(default=None, description="Electronic-state label")
    multiplicity: int | None = Field(default=None, description="Spin multiplicity")
    spin: float | None = Field(default=None, description="Spin quantum number")
    irrep: str | None = Field(default=None, description="Irreducible representation")
    method: str | None = Field(default=None, description="Method used for this state")
    energy: PlainQuantity | None = Field(default=None, description="State total energy")
    excitation_energy: PlainQuantity | None = Field(default=None, description="Excitation energy")
    oscillator_strength: float | None = Field(default=None, description="Oscillator strength")
    transition_dipole: PintArray3 | None = Field(
        default=None, description="Transition dipole vector"
    )
    configurations: list[ElectronicConfiguration] = Field(
        default_factory=list, description="Dominant electronic configurations"
    )
    properties: dict[str, Any] = Field(
        default_factory=dict, description="Additional state-specific properties"
    )
    source: str | None = Field(default=None, description="Program/source section")


class ElectronicStates(BaseDataClassWithUnit, Sequence[ElectronicState]):
    states: list[ElectronicState] = Field(
        default_factory=list, description="Electronic states ordered by source output"
    )

    def __iter__(self) -> Iterator[ElectronicState]:  # type: ignore[override]
        return iter(self.states)

    def __len__(self) -> int:
        return len(self.states)

    @overload
    def __getitem__(self, index: int) -> ElectronicState: ...

    @overload
    def __getitem__(self, index: slice) -> list[ElectronicState]: ...

    def __getitem__(self, index: int | slice) -> ElectronicState | list[ElectronicState]:
        return self.states[index]

    @property
    def ground_state(self) -> ElectronicState | None:
        return self.states[0] if self.states else None

    @property
    def excited_states(self) -> list[ElectronicState]:
        return self.states[1:]

    def to_state_table(self) -> PropertyTable:
        columns: dict[str, PropertyColumnValue] = {
            "state_index": [
                idx if state.state_index is None else state.state_index
                for idx, state in enumerate(self.states)
            ],
            "root": [state.root for state in self.states],
            "label": [state.label for state in self.states],
            "multiplicity": [state.multiplicity for state in self.states],
            "spin": [state.spin for state in self.states],
            "irrep": [state.irrep for state in self.states],
            "method": [state.method for state in self.states],
            "oscillator_strength": [state.oscillator_strength for state in self.states],
        }

        energies = [state.energy for state in self.states]
        if any(energy is not None for energy in energies):
            columns["energy"] = (
                np.array(
                    [
                        np.nan if energy is None else energy.to(atom_ureg.hartree).magnitude
                        for energy in energies
                    ]
                )
                * atom_ureg.hartree
            )

        excitation_energies = [state.excitation_energy for state in self.states]
        if any(excitation_energy is not None for excitation_energy in excitation_energies):
            columns["excitation_energy"] = (
                np.array(
                    [
                        np.nan
                        if excitation_energy is None
                        else excitation_energy.to(atom_ureg.eV).magnitude
                        for excitation_energy in excitation_energies
                    ]
                )
                * atom_ureg.eV
            )

        return PropertyTable(
            columns=columns,
            row_labels=[
                state.label or str(state.root) if state.root is not None else str(idx)
                for idx, state in enumerate(self.states)
            ],
            metadata={"source": "ElectronicStates", "index_base": 0},
        )

    def to_transitions(self) -> list[PropertyTransition]:
        transitions: list[PropertyTransition] = []
        ground_state = self.ground_state
        initial_state = ground_state.state_index if ground_state else 0
        for idx, state in enumerate(self.states):
            if state is ground_state:
                continue
            final_state = state.state_index if state.state_index is not None else idx
            transitions.append(
                PropertyTransition(
                    label=state.label,
                    initial_state=initial_state,
                    final_state=final_state,
                    energy=state.excitation_energy,
                    oscillator_strength=state.oscillator_strength,
                    transition_dipole=state.transition_dipole,
                    properties={
                        "root": state.root,
                        "multiplicity": state.multiplicity,
                        "spin": state.spin,
                        "irrep": state.irrep,
                        "method": state.method,
                        "source": state.source,
                    },
                )
            )
        return transitions

    def to_spectrum(self, label: str = "Electronic transitions") -> Spectrum:
        transitions = self.to_transitions()
        return Spectrum(
            label=label,
            x_label="excitation energy",
            y_label="oscillator strength",
            transitions=transitions,
            metadata={"source": "ElectronicStates"},
        )

    def to_property_bundle(self) -> PropertyBundle:
        transitions = self.to_transitions()
        return PropertyBundle(
            tables={"electronic_states": self.to_state_table()} if self.states else {},
            spectra={"electronic_transitions": self.to_spectrum()} if transitions else {},
            transitions={"electronic": transitions} if transitions else {},
            metadata={"source": "ElectronicStates"},
        )


class MultireferenceResult(BaseDataClassWithUnit):
    method: str | None = Field(default=None, description="Multi-reference method")
    reference_method: str | None = Field(default=None, description="Reference method")
    ci_type: str | None = Field(default=None, description="Configuration interaction type")
    active_space: ActiveSpace | None = Field(default=None, description="Active-space summary")
    electronic_states: ElectronicStates | None = Field(
        default=None, description="State/root resolved multi-reference results"
    )
    corrections: dict[str, PlainQuantity | float | str | None] = Field(
        default_factory=dict, description="Energy corrections or named post-CI corrections"
    )
    diagnostics: list[str] = Field(default_factory=list, description="Parser or method diagnostics")
    properties: dict[str, Any] = Field(
        default_factory=dict, description="Additional multi-reference properties"
    )

    def to_property_bundle(self) -> PropertyBundle:
        scalar_properties: dict[str, PropertyScalarValue] = {
            f"{key}_correction": value
            for key, value in self.corrections.items()
            if value is not None
        }
        bundle = PropertyBundle(
            scalar_properties=scalar_properties,
            metadata={
                "source": "MultireferenceResult",
                "method": self.method,
                "reference_method": self.reference_method,
                "ci_type": self.ci_type,
                "diagnostics": self.diagnostics,
            },
        )
        if self.electronic_states is not None:
            states_bundle = self.electronic_states.to_property_bundle()
            bundle.tables.update(states_bundle.tables)
            bundle.spectra.update(states_bundle.spectra)
            bundle.transitions.update(states_bundle.transitions)
        return bundle


__all__ = [
    "ElectronicConfiguration",
    "ElectronicState",
    "ElectronicStates",
    "MultireferenceResult",
]
