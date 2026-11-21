"""mesh-uv-weld: 网格顶点焊接工具包（出题用初始失败环境）。

公开契约见 :mod:`weld.vertices`。
"""

from .vertices import (
    EPS,
    WeldError,
    compact,
    key_for,
    merge_uv_islands,
    preserve_seams,
    remove_degenerate,
    validate,
    weld_batch,
    weld_positions,
)

__all__ = [
    "EPS",
    "WeldError",
    "key_for",
    "weld_positions",
    "preserve_seams",
    "remove_degenerate",
    "compact",
    "merge_uv_islands",
    "validate",
    "weld_batch",
]
