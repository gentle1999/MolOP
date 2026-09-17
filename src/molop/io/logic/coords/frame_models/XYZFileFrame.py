"""
Author: TMJ
Date: 2025-07-28 23:05:56
LastEditors: TMJ
LastEditTime: 2026-06-19 00:57:34
Description: 请填写简介
"""

from typing import cast

from pydantic import Field, model_validator
from typing_extensions import Self

from molop.io.base_models.ChemFileFrame import BaseCoordsFrame, _HasCoords
from molop.io.base_models.DataClasses import Comment
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.io.logic.coords.frame_models._coords_renderers import render_xyz_frame


class XYZFileFrameMixin:
    comment: str = Field(default="", description="comment")

    @model_validator(mode="after")
    def sync_common_comments(self) -> Self:
        typed_self = cast(BaseCoordsFrame, self)
        if self.comment:
            typed_self.comments.ensure(
                Comment(
                    text=self.comment,
                    source_format="xyz",
                )
            )
        elif common_comment := typed_self.comments.first_text(kind="comment"):
            self.comment = common_comment
        return self

    def _render(self, comment: str | None = None, **kwargs) -> str:
        """Render the XYZ file frame as a string.

        Args:
            comment (str | None): The comment to use. Defaults to None.

        Returns:
            str: The rendered XYZ file frame.
        """
        typed_self = cast(_HasCoords, self)
        comments = cast(BaseCoordsFrame, self).comments
        stored_comment = self.comment or comments.first_text(kind="comment")
        return render_xyz_frame(typed_self, comment=comment, stored_comment=stored_comment)


class XYZFileFrameMemory(
    MemoryStorageMixin, XYZFileFrameMixin, BaseCoordsFrame["XYZFileFrameMemory"]
): ...


class XYZFileFrameDisk(
    DiskStorageMixin, XYZFileFrameMixin, BaseCoordsFrame["XYZFileFrameDisk"]
): ...
