"""
Author: TMJ
Date: 2024-01-01 21:40:35
LastEditors: TMJ
LastEditTime: 2024-01-01 21:42:59
Description: 请填写简介
"""
from pint import UnitRegistry, set_application_registry

si_ureg = UnitRegistry(system="SI")
atom_ureg = UnitRegistry(system="atomic")

set_application_registry(si_ureg)
set_application_registry(atom_ureg)
