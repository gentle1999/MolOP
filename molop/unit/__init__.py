"""
Author: TMJ
Date: 2024-01-01 21:40:35
LastEditors: TMJ
LastEditTime: 2024-10-18 19:10:30
Description: 请填写简介
"""

from typing import Union

from pint import UnitRegistry, set_application_registry
from pint.facets.plain import PlainQuantity, PlainUnit

si_ureg = UnitRegistry(system="SI")
atom_ureg = UnitRegistry(system="atomic")

set_application_registry(si_ureg)
set_application_registry(atom_ureg)


def unit_transform(value: Union[PlainQuantity, None], unit: PlainUnit):
    if value is None:
        return None
    if value.units == unit:
        return value
    return value.to(unit)
