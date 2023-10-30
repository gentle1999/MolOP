'''
Author: TMJ
Date: 2023-03-22 19:55:04
LastEditors: TMJ
LastEditTime: 2023-03-24 19:28:36
Description: 请填写简介
'''
class ArgumentError(Exception):
    'Error when Argument loss'

class MolError(Exception):
    'Mol inherent errors'

class ReactError(Exception):
    'Reaction inherent errors'