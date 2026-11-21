"""mesh-uv-weld 的验收测试：10 个缺陷各 1 条 + 5 条组合/边界。

全部用例都直接调用 `weld.vertices` 的桩函数，因此初始状态下 15 条全红，且没有
collection / import / syntax error。测试只依赖标准库 math 与 pytest。
"""

import math

import pytest

from weld.vertices import (
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


# --------------------------------------------------------------------------- #
# 网格辅助：只构造数据，不参与实现
# --------------------------------------------------------------------------- #


def make_mesh(positions, indices, uvs=None, normals=None):
    """按固定四字段构造网格，缺失的 UV/法线用确定值补齐。"""
    if uvs is None:
        uvs = [[0.0, 0.0] for _ in positions]
    if normals is None:
        normals = [[0.0, 0.0, 1.0] for _ in positions]
    return {
        "positions": [list(p) for p in positions],
        "uvs": [list(uv) for uv in uvs],
        "normals": [list(n) for n in normals],
        "indices": [list(tri) for tri in indices],
    }


def triangle_mesh(x=0.0):
    """XY 平面上的单个三角形，整体沿 X 平移 x。"""
    return make_mesh(
        [[x, 0.0, 0.0], [x + 1.0, 0.0, 0.0], [x, 1.0, 0.0]],
        [[0, 1, 2]],
        uvs=[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
    )


# --------------------------------------------------------------------------- #
# 缺陷 1：位置量化用 round 而不是 floor(x / epsilon)
# --------------------------------------------------------------------------- #


def test_key_for_uses_floor_on_cell_boundary():
    assert key_for((0.999999, 0.0, 0.0), 1.0) == (0, 0, 0)
    assert key_for((1.0, 0.0, 0.0), 1.0) == (1, 0, 0)
    assert key_for((-0.000001, 0.0, 0.0), 1.0) == (-1, 0, 0)


# --------------------------------------------------------------------------- #
# 缺陷 3：合并后索引不重映射
# --------------------------------------------------------------------------- #


def test_weld_positions_remaps_indices_after_merge():
    m = make_mesh(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [9.0, 9.0, 9.0],
        ],
        [[0, 1, 2], [3, 2, 4]],
    )
    out = weld_positions(m, EPS)

    assert out["positions"] == [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
        [9.0, 9.0, 9.0],
    ]
    assert out["indices"] == [[0, 1, 2], [0, 2, 3]]
    assert out["remap"] == [0, 1, 2, 0, 3, 4]


# --------------------------------------------------------------------------- #
# 缺陷 2：合并不看 UV，UV 缝被焊死
# --------------------------------------------------------------------------- #


def test_preserve_seams_keeps_vertices_with_different_uv():
    m = make_mesh(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        [[0, 1, 2], [3, 4, 5]],
        uvs=[
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [0.0, 0.0],
            [0.5, 0.0],
            [1.0, 1.0],
        ],
    )
    out = preserve_seams(m, EPS)

    # 0 与 3 位置、UV 都相同 -> 焊接；1/4、2/5 位置相同但 UV 不同 -> 保留缝
    assert out["positions"] == [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ]
    assert out["uvs"] == [
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [0.5, 0.0],
        [1.0, 1.0],
    ]
    assert out["indices"] == [[0, 1, 2], [0, 3, 4]]
    assert out["remap"] == [0, 1, 2, 0, 3, 4]


# --------------------------------------------------------------------------- #
# 缺陷 4：退化三角形不剔除，法线为 NaN
# --------------------------------------------------------------------------- #


def test_remove_degenerate_drops_zero_area_triangles():
    m = make_mesh(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        [[0, 1, 2], [0, 1, 3]],
    )
    out = remove_degenerate(m)

    assert out["indices"] == [[0, 1, 3]]
    # 本函数只删面，不删顶点；未引用顶点交给 compact。
    assert len(out["positions"]) == 4
    assert len(out["uvs"]) == 4
    assert len(out["normals"]) == 4


# --------------------------------------------------------------------------- #
# 缺陷 5：未引用顶点不删，文件变大
# --------------------------------------------------------------------------- #


def test_compact_drops_unreferenced_vertices():
    m = make_mesh(
        [
            [0.0, 0.0, 0.0],
            [9.0, 9.0, 9.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [8.0, 8.0, 8.0],
        ],
        [[0, 2, 3]],
    )
    out = compact(m)

    assert out["positions"] == [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ]
    assert out["indices"] == [[0, 1, 2]]
    assert out["remap"] == [0, -1, 1, 2, -1]


# --------------------------------------------------------------------------- #
# 缺陷 6：UV 岛判定用包围盒相交，误合并
# --------------------------------------------------------------------------- #


def test_merge_uv_islands_ignores_overlapping_bboxes():
    m = make_mesh(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            # 下面三个顶点构成另一个拓扑岛，位置/UV 与第一个三角形完全相同，
            # 包围盒也完全重叠，但没有任何共享索引。
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        [[0, 1, 2], [3, 2, 1], [4, 5, 6]],
        uvs=[
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
        ],
    )
    out = merge_uv_islands(m, EPS)

    # 岛内 3 与 0 焊接；第三个三角形是独立岛，即使包围盒重叠也不合并。
    assert out["positions"] == [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ]
    assert out["indices"] == [[0, 1, 2], [0, 2, 1], [3, 4, 5]]
    assert out["remap"] == [0, 1, 2, 0, 3, 4, 5]


# --------------------------------------------------------------------------- #
# 缺陷 7：索引越界不校验
# --------------------------------------------------------------------------- #


def test_validate_reports_out_of_range_index():
    m = make_mesh(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        [[0, 1, 9]],
    )
    problems = validate(m)

    assert isinstance(problems, list)
    assert problems
    assert any("index" in problem.lower() for problem in problems)


# --------------------------------------------------------------------------- #
# 缺陷 8：NaN 坐标不报错
# --------------------------------------------------------------------------- #


def test_validate_reports_nan_position():
    m = make_mesh(
        [[0.0, 0.0, 0.0], [math.nan, 0.0, 0.0], [0.0, 1.0, 0.0]],
        [[0, 1, 2]],
    )
    problems = validate(m)

    assert problems
    assert any("nan" in problem.lower() for problem in problems)


# --------------------------------------------------------------------------- #
# 缺陷 9：批量共用 remap 表，网格之间串
# --------------------------------------------------------------------------- #


def test_weld_batch_keeps_remaps_independent_and_rolls_back():
    mesh_a = triangle_mesh(0.0)
    # 前两个顶点在 epsilon=0.5 下同格，但输入三角形本身不退化为零面积。
    mesh_b = make_mesh(
        [[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [0.0, 1.0, 0.0]],
        [[0, 1, 2]],
        uvs=[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
    )
    config = {"epsilon": 0.5}

    out = weld_batch([mesh_a, mesh_b], config)

    assert len(out["meshes"]) == 2
    assert out["remaps"] == [[0, 1, 2], [0, 0, 1]]
    assert out["meshes"][0]["positions"] == triangle_mesh(0.0)["positions"]
    assert out["meshes"][1]["positions"] == [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    assert out["meshes"][1]["indices"] == [[0, 0, 1]]

    bad = make_mesh(
        [[10.0, 0.0, 0.0], [11.0, 0.0, 0.0], [10.0, 1.0, 0.0]],
        [[0, 1, 9]],
    )
    with pytest.raises(WeldError):
        weld_batch([mesh_a, bad], config)
    # 整批失败时不能原地改写已经处理过的网格。
    assert mesh_a == triangle_mesh(0.0)


# --------------------------------------------------------------------------- #
# 缺陷 10：epsilon <= 0 不校验，除零
# --------------------------------------------------------------------------- #


def test_key_for_rejects_non_positive_epsilon():
    with pytest.raises(WeldError):
        key_for((0.0, 0.0, 0.0), -1.0)


# --------------------------------------------------------------------------- #
# 组合 / 边界 1：单三角形焊接后不变
# --------------------------------------------------------------------------- #


def test_single_triangle_weld_is_unchanged():
    m = triangle_mesh(0.0)
    out = weld_positions(m, EPS)

    assert out == {
        "positions": [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        "indices": [[0, 1, 2]],
        "remap": [0, 1, 2],
    }


# --------------------------------------------------------------------------- #
# 组合 / 边界 2：空网格返回空
# --------------------------------------------------------------------------- #


def test_empty_mesh_welds_to_empty():
    m = make_mesh([], [])
    out = weld_positions(m, EPS)

    assert out == {"positions": [], "indices": [], "remap": []}


# --------------------------------------------------------------------------- #
# 组合 / 边界 3：epsilon = 0 抛 WeldError
# --------------------------------------------------------------------------- #


def test_weld_positions_rejects_zero_epsilon():
    with pytest.raises(WeldError):
        weld_positions(triangle_mesh(0.0), 0.0)


# --------------------------------------------------------------------------- #
# 组合 / 边界 4：焊接后三角形数不增
# --------------------------------------------------------------------------- #


def test_welding_does_not_increase_triangle_count():
    m = make_mesh(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
        ],
        [[0, 1, 2], [3, 2, 4]],
    )
    out = weld_positions(m, EPS)

    assert len(out["indices"]) == len(m["indices"]) == 2
    assert all(0 <= index < len(out["positions"]) for tri in out["indices"] for index in tri)


# --------------------------------------------------------------------------- #
# 组合 / 边界 5：批量全成功时顺序不变
# --------------------------------------------------------------------------- #


def test_weld_batch_preserves_mesh_order():
    meshes = [triangle_mesh(0.0), triangle_mesh(10.0), triangle_mesh(20.0)]
    out = weld_batch(meshes, {"epsilon": EPS})

    assert len(out["meshes"]) == 3
    assert out["remaps"] == [[0, 1, 2], [0, 1, 2], [0, 1, 2]]
    assert [result["positions"][0][0] for result in out["meshes"]] == [
        0.0,
        10.0,
        20.0,
    ]
