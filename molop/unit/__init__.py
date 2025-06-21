"""
Author: TMJ
Date: 2024-10-19 09:57:26
LastEditors: TMJ
LastEditTime: 2025-06-16 12:32:15
Description: 请填写简介
"""

from typing import Union

from pint import UnitRegistry, set_application_registry
from pint.facets.numpy import NumpyRegistry
from pint.facets.plain import PlainQuantity, PlainUnit

si_ureg = UnitRegistry(system="SI")
atom_ureg = UnitRegistry(system="atomic")
numpy_ureg = UnitRegistry(system="atomic", force_ndarray=True)
set_application_registry(si_ureg)
set_application_registry(atom_ureg)
set_application_registry(numpy_ureg)


def unit_transform(
    value: Union[PlainQuantity, None], unit: Union[str, PlainUnit]
) -> Union[PlainQuantity, None]:
    """
    Transform the unit of a quantity.

    Parameters:
        value (Union[PlainQuantity, None]): The quantity to be transformed.
        unit (Union[str, PlainUnit]): The target unit.

    Returns:
        Union[PlainQuantity, None]: The transformed quantity.
    """
    if value is None:
        return None
    if value.units == unit:
        return value
    return value.to(unit)
