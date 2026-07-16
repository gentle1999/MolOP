"""
Author: TMJ
Date: 2023-10-30 14:05:04
LastEditors: TMJ
LastEditTime: 2025-12-14 16:45:09
Description: 请填写简介
"""

from dataclasses import dataclass
from typing import Annotated, Any, TypeAlias, TypeVar

import numpy as np
import numpy.typing as npt
from openbabel import pybel
from pint.facets.numpy.quantity import NumpyQuantity
from pydantic import AfterValidator
from rdkit import Chem


RdMol: TypeAlias = Chem.rdchem.Mol
RWMol: TypeAlias = Chem.rdchem.RWMol
OMol: TypeAlias = pybel.Molecule
RdConformer: TypeAlias = Chem.rdchem.Conformer
DType = TypeVar("DType", bound=np.generic)
ShapeDimension: TypeAlias = int | str | None


@dataclass(frozen=True)
class ArrayShape:
    """Pydantic validator for fixed, wildcard, or repeated array dimensions.

    Integer dimensions are exact, ``None`` accepts any size, and equal string
    labels must resolve to the same size. For example, ``ArrayShape("N", 3)``
    accepts Cartesian coordinate arrays and ``ArrayShape("N", "N")`` accepts
    square matrices.
    """

    dimensions: tuple[ShapeDimension, ...]

    def __init__(self, *dimensions: ShapeDimension) -> None:
        if not dimensions:
            raise ValueError("ArrayShape requires at least one dimension")
        if any(isinstance(dimension, int) and dimension < 0 for dimension in dimensions):
            raise ValueError("ArrayShape integer dimensions must be non-negative")
        if any(isinstance(dimension, str) and not dimension for dimension in dimensions):
            raise ValueError("ArrayShape dimension labels must not be empty")
        object.__setattr__(self, "dimensions", tuple(dimensions))

    def __call__(self, value: Any) -> Any:
        actual_shape = getattr(value, "shape", None)
        if actual_shape is None:
            raise ValueError("value must expose a shape attribute")
        actual = tuple(int(size) for size in actual_shape)
        if len(actual) != len(self.dimensions):
            raise ValueError(f"value must have shape {self}; got {actual}")

        bound_dimensions: dict[str, int] = {}
        for expected, size in zip(self.dimensions, actual, strict=True):
            if expected is None:
                continue
            if isinstance(expected, int):
                if size != expected:
                    raise ValueError(f"value must have shape {self}; got {actual}")
                continue
            if expected in bound_dimensions and bound_dimensions[expected] != size:
                raise ValueError(f"value must have shape {self}; got {actual}")
            bound_dimensions[expected] = size
        return value

    def __str__(self) -> str:
        return f"({', '.join(str(item) for item in self.dimensions)})"


@dataclass(frozen=True)
class ArrayRank:
    """Pydantic validator for arrays whose exact dimensions are intentionally open."""

    minimum: int = 1
    maximum: int | None = None

    def __post_init__(self) -> None:
        if self.minimum < 0:
            raise ValueError("ArrayRank minimum must be non-negative")
        if self.maximum is not None and self.maximum < self.minimum:
            raise ValueError("ArrayRank maximum must be greater than or equal to minimum")

    def __call__(self, value: Any) -> Any:
        actual_shape = getattr(value, "shape", None)
        if actual_shape is None:
            raise ValueError("value must expose a shape attribute")
        rank = len(actual_shape)
        if rank < self.minimum or (self.maximum is not None and rank > self.maximum):
            expected = (
                f"at least {self.minimum}"
                if self.maximum is None
                else f"between {self.minimum} and {self.maximum}"
            )
            raise ValueError(f"value must have {expected} dimensions; got {rank}")
        return value


Array: TypeAlias = Annotated[npt.NDArray[DType], AfterValidator(ArrayRank())]
ArrayN: TypeAlias = Annotated[npt.NDArray[DType], AfterValidator(ArrayShape("N"))]
ArrayNx3: TypeAlias = Annotated[npt.NDArray[DType], AfterValidator(ArrayShape("N", 3))]
Array3: TypeAlias = Annotated[npt.NDArray[DType], AfterValidator(ArrayShape(3))]
Array3x3: TypeAlias = Annotated[npt.NDArray[DType], AfterValidator(ArrayShape(3, 3))]
Array4x4: TypeAlias = Annotated[npt.NDArray[DType], AfterValidator(ArrayShape(4, 4))]
SquareArray: TypeAlias = Annotated[npt.NDArray[DType], AfterValidator(ArrayShape("N", "N"))]

PintArray: TypeAlias = Annotated[NumpyQuantity, AfterValidator(ArrayRank())]
PintArrayN: TypeAlias = Annotated[NumpyQuantity, AfterValidator(ArrayShape("N"))]
PintArrayNx3: TypeAlias = Annotated[NumpyQuantity, AfterValidator(ArrayShape("N", 3))]
PintArray3: TypeAlias = Annotated[NumpyQuantity, AfterValidator(ArrayShape(3))]
PintArray6: TypeAlias = Annotated[NumpyQuantity, AfterValidator(ArrayShape(6))]
PintArray10: TypeAlias = Annotated[NumpyQuantity, AfterValidator(ArrayShape(10))]
PintArray15: TypeAlias = Annotated[NumpyQuantity, AfterValidator(ArrayShape(15))]
PintArray3x3: TypeAlias = Annotated[NumpyQuantity, AfterValidator(ArrayShape(3, 3))]
PintSquareMatrix: TypeAlias = Annotated[NumpyQuantity, AfterValidator(ArrayShape("N", "N"))]
PintArray6Or3x3: TypeAlias = PintArray6 | PintArray3x3

# Backward-compatible aliases. Their old definitions did not enforce shape.
arrayNx3: TypeAlias = ArrayNx3[DType]
arrayN: TypeAlias = ArrayN[DType]
