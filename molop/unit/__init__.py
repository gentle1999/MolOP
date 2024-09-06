"""
Author: TMJ
Date: 2024-01-01 21:40:35
LastEditors: TMJ
LastEditTime: 2024-09-06 11:48:52
Description: 请填写简介
"""

from typing import Union

from pint import UnitRegistry, set_application_registry
from pint.facets.plain import PlainQuantity

si_ureg = UnitRegistry(system="SI")
atom_ureg = UnitRegistry(system="atomic")

set_application_registry(si_ureg)
set_application_registry(atom_ureg)


def unit_transform(value: Union[PlainQuantity, None], unit):
    if value is None:
        return None
    return value.to(unit)
