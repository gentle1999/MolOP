"""
Author: TMJ
Date: 2024-10-19 09:57:26
LastEditors: TMJ
LastEditTime: 2025-12-16 01:12:17
Description: 请填写简介
"""

from pint import UnitRegistry, set_application_registry
from pint._typing import UnitLike
from pint.facets.plain import PlainQuantity


si_ureg = UnitRegistry(system="SI")
atom_ureg = UnitRegistry(system="atomic")
set_application_registry(si_ureg)
set_application_registry(atom_ureg)


def unit_transform(value: PlainQuantity | None, unit: UnitLike) -> PlainQuantity | None:
    """
    Transform the unit of a quantity.

    Parameters:
        value (Optional[PlainQuantity]): The quantity to be transformed.
        unit (UnitLike): The target unit.

    Returns:
        Optional[PlainQuantity]: The transformed quantity.
    """
    if value is None:
        return None
    if value.units == unit:
        return value
    return value.to(unit)
