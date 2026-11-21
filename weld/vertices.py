"""网格顶点焊接核心。

约定
----
* 网格是普通 dict：``{"positions": [[x, y, z], ...], "uvs": [[u, v], ...],
  "normals": [[nx, ny, nz], ...], "indices": [[i, j, k], ...]}``。
  ``positions``、``uvs``、``normals`` 按顶点下标一一对齐；``indices`` 的每个
  三元组都指向这三个数组中的同一顶点下标。
* 位置比较用调用方传入的十进制步长 ``epsilon``：每个分量按
  ``floor(x / epsilon)`` 量化，量化键相同才视为同一位置。``epsilon`` 必须
  是有限正数；``epsilon <= 0`` 一律抛 :class:`WeldError`，绝不能除零。
* UV 缝必须保留：位置相同但 UV 不同的顶点不能合并。UV 各分量比较的全局
  容差是 ``1e-6``。
* 几何面积、长度、坐标有限性检查的全局容差统一为 :data:`EPS`（``1e-6``）。
* UV 岛是三角形按公共边（共享同一对顶点索引）得到的连通分量；UV 缝在网格
  里已经表现为重复顶点，因此缝两侧不共享索引、属于不同岛。禁止用 UV 包围盒
  是否相交来判断岛，否则包围盒重叠但拓扑分离的岛会被误合并。
* 只用标准库 ``math``，不使用随机数，相同输入必须得到相同输出。
* 非法输入（缺少字段、数组长度不一致、下标越界、非有限坐标、非正 epsilon）
  抛 :class:`WeldError`；``validate`` 是例外，它把问题收集成字符串列表返回。
"""

import math

__all__ = [
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

#: 全局几何/UV 比较容差。
EPS = 1e-6

#: 网格的四个必需字段。
_FIELDS = ("positions", "uvs", "normals", "indices")


class WeldError(Exception):
    """顶点焊接中的非法输入、非法参数或批量操作失败时抛出。

    约定：
    * ``epsilon <= 0``、``epsilon`` 非有限、网格字段缺失或数组长度不一致时抛出；
    * 顶点下标越界、坐标/法线/UV 含 NaN 或无穷大时抛出；
    * 批量焊接中任一网格非法时抛出，且不能返回或写出任何部分结果（整批回滚）。
    """


# --------------------------------------------------------------------------- #
# 内部校验辅助
# --------------------------------------------------------------------------- #


def _check_epsilon(epsilon):
    """校验量化步长并返回 float。"""
    try:
        value = float(epsilon)
    except (TypeError, ValueError):
        raise WeldError("epsilon must be a finite positive number")
    if not math.isfinite(value) or value <= 0.0:
        raise WeldError("epsilon must be a finite positive number")
    return value


def _check_number(value, what):
    """把单个分量转成有限 float，否则抛 WeldError。"""
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise WeldError("%s is not a number" % what)
    if not math.isfinite(number):
        raise WeldError("%s is not finite" % what)
    return number


def _check_vector(vector, size, what):
    """校验定长数值向量，返回 float 列表。"""
    if not isinstance(vector, (list, tuple)) or len(vector) != size:
        raise WeldError("%s must have %d components" % (what, size))
    return [_check_number(component, "%s[%d]" % (what, i))
            for i, component in enumerate(vector)]


def _check_mesh(mesh):
    """校验网格结构与数值，返回 ``(positions, uvs, normals, triangles)``。

    ``triangles`` 是 ``tuple[int, int, int]`` 的列表。返回的顶点数组仍是输入
    对象，实现只读不写。
    """
    if not isinstance(mesh, dict):
        raise WeldError("mesh must be a dict")
    arrays = {}
    for name in _FIELDS:
        if name not in mesh:
            raise WeldError("mesh is missing field %r" % name)
        if not isinstance(mesh[name], (list, tuple)):
            raise WeldError("mesh field %r must be a list" % name)
        arrays[name] = mesh[name]

    positions = arrays["positions"]
    uvs = arrays["uvs"]
    normals = arrays["normals"]
    count = len(positions)
    if len(uvs) != count:
        raise WeldError("uvs length %d does not match positions length %d"
                        % (len(uvs), count))
    if len(normals) != count:
        raise WeldError("normals length %d does not match positions length %d"
                        % (len(normals), count))

    for i, position in enumerate(positions):
        _check_vector(position, 3, "positions[%d]" % i)
    for i, uv in enumerate(uvs):
        _check_vector(uv, 2, "uvs[%d]" % i)
    for i, normal in enumerate(normals):
        _check_vector(normal, 3, "normals[%d]" % i)

    triangles = []
    for triangle_index, triangle in enumerate(arrays["indices"]):
        if not isinstance(triangle, (list, tuple)) or len(triangle) != 3:
            raise WeldError("indices[%d] is not a triangle" % triangle_index)
        parsed = []
        for corner, index in enumerate(triangle):
            if not isinstance(index, int):
                raise WeldError("indices[%d][%d] is not an integer index"
                                % (triangle_index, corner))
            if index < 0 or index >= count:
                raise WeldError("indices[%d][%d] index %d out of range of %d vertices"
                                % (triangle_index, corner, index, count))
            parsed.append(index)
        triangles.append(tuple(parsed))

    return positions, uvs, normals, triangles


def _quantize(value, epsilon):
    """安全的 ``floor(value / epsilon)``。"""
    try:
        return math.floor(value / epsilon)
    except (OverflowError, ValueError):
        raise WeldError("position cannot be quantized with the given epsilon")


def _position_key(position, epsilon):
    """把已校验的三维位置量化成整数键（向下取整）。"""
    return (
        _quantize(float(position[0]), epsilon),
        _quantize(float(position[1]), epsilon),
        _quantize(float(position[2]), epsilon),
    )


def _triangle_area(positions, triangle):
    """三维三角形面积 ``0.5 * |(b - a) x (c - a)|``。"""
    ax, ay, az = (float(c) for c in positions[triangle[0]])
    bx, by, bz = (float(c) for c in positions[triangle[1]])
    cx, cy, cz = (float(c) for c in positions[triangle[2]])
    ux, uy, uz = bx - ax, by - ay, bz - az
    vx, vy, vz = cx - ax, cy - ay, cz - az
    nx = uy * vz - uz * vy
    ny = uz * vx - ux * vz
    nz = ux * vy - uy * vx
    return 0.5 * math.sqrt(nx * nx + ny * ny + nz * nz)


def _copy_vectors(vectors):
    """逐个复制向量，避免输出与输入共享可变列表。"""
    return [list(vector) for vector in vectors]


def _find_uv_group(groups, u, v):
    """在已登记的 UV 代表中找容差范围内匹配的，返回其新下标或 ``-1``。"""
    for candidate_u, candidate_v, new_index in groups:
        if (abs(u - candidate_u) <= EPS
                and abs(v - candidate_v) <= EPS):
            return new_index
    return -1


# --------------------------------------------------------------------------- #
# 公开接口
# --------------------------------------------------------------------------- #


def key_for(pos, epsilon):
    """把一个三维位置量化成整数键。

    参数
    ----
    pos : sequence[float]
        长度为 3 的坐标 ``[x, y, z]``。
    epsilon : float
        位置量化步长，必须是有限正数。

    返回
    ----
    tuple[int, int, int]
        ``(floor(x / epsilon), floor(y / epsilon), floor(z / epsilon))``。
        必须使用向下取整而不是 ``round``：例如 ``epsilon=1.0`` 时，
        ``0.999999`` 属于键 ``0``，而 ``1.0`` 才属于键 ``1``；负坐标同样
        向下取整（``-0.000001`` 属于键 ``-1``）。

    异常
    ----
    WeldError
        ``epsilon`` 不是有限正数；``pos`` 不是长度为 3 的有限坐标。
    """
    step = _check_epsilon(epsilon)
    values = _check_vector(pos, 3, "pos")
    return tuple(_quantize(value, step) for value in values)


def weld_positions(mesh, epsilon):
    """只按位置焊接顶点，忽略 UV（用于基础焊接/缺陷诊断）。

    参数
    ----
    mesh : dict
        见模块说明。本函数不读取 UV 或法线来决定是否合并。
    epsilon : float
        位置量化步长，必须大于 0。

    返回
    ----
    dict
        ``{"positions": [[x, y, z], ...], "indices": [[i, j, k], ...],
        "remap": [new_index, ...]}``。

        * 位置量化键相同的顶点合并为一个，保留输入顺序中**第一次出现**的
          代表位置；
        * ``indices`` 中每个旧下标都必须替换为 ``remap[旧下标]``，不能留下
          仍指向已合并顶点的索引；
        * ``remap`` 与输入顶点等长，``remap[i]`` 是旧顶点 ``i`` 的新下标；
        * 空网格返回 ``{"positions": [], "indices": [], "remap": []}``。

    异常
    ----
    WeldError
        网格缺字段/长度不一致、``epsilon <= 0``、下标越界或坐标非有限。
    """
    step = _check_epsilon(epsilon)
    positions, _uvs, _normals, triangles = _check_mesh(mesh)

    representatives = {}
    merged_positions = []
    remap = []
    for old_index, position in enumerate(positions):
        key = _position_key(position, step)
        new_index = representatives.get(key)
        if new_index is None:
            new_index = len(merged_positions)
            representatives[key] = new_index
            merged_positions.append(list(position))
        remap.append(new_index)

    indices = [[remap[index] for index in triangle] for triangle in triangles]
    return {"positions": merged_positions, "indices": indices, "remap": remap}


def preserve_seams(mesh, epsilon):
    """按“位置 + UV”焊接顶点，并把 UV 缝保留下来。

    位置量化键相同**且** UV 两个分量都在 ``1e-6`` 内相同，才允许合并。
    只要 UV 有一个分量相差超过 ``1e-6``，即使位置完全相同也必须保留为两个
    顶点，否则缝会被焊死、贴图会整体撕裂。

    参数
    ----
    mesh : dict
        见模块说明，``uvs`` 必须与 ``positions`` 等长。
    epsilon : float
        位置量化步长，必须大于 0。

    返回
    ----
    dict
        ``{"positions": [...], "uvs": [...], "indices": [...],
        "remap": [...]}``。合并键为 ``(位置量化键, UV 容差分组)``，每个键保留
        第一次出现的顶点；``indices`` 必须按 ``remap`` 改写；``remap`` 与输入
        顶点等长。空网格返回各字段均为 ``[]`` 的同构字典。

    异常
    ----
    WeldError
        网格字段缺失/长度不一致、``epsilon <= 0``、下标越界、UV 数量不匹配
        或坐标非有限。
    """
    step = _check_epsilon(epsilon)
    positions, uvs, _normals, triangles = _check_mesh(mesh)

    # 每个位置键 -> [(u, v, 新下标), ...]，逐个与代表 UV 做容差比较。
    groups = {}
    merged_positions = []
    merged_uvs = []
    remap = []
    for old_index, position in enumerate(positions):
        key = _position_key(position, step)
        u, v = float(uvs[old_index][0]), float(uvs[old_index][1])
        candidates = groups.setdefault(key, [])
        new_index = _find_uv_group(candidates, u, v)
        if new_index < 0:
            new_index = len(merged_positions)
            candidates.append((u, v, new_index))
            merged_positions.append(list(position))
            merged_uvs.append([u, v])
        remap.append(new_index)

    indices = [[remap[index] for index in triangle] for triangle in triangles]
    return {
        "positions": merged_positions,
        "uvs": merged_uvs,
        "indices": indices,
        "remap": remap,
    }


def remove_degenerate(mesh):
    """删除退化三角形，保留其余三角形和全部顶点。

    退化判定：三角形三维面积 ``<= EPS``（``1e-6``），包括三个顶点共线、
    两个位置重合或三个位置重合。删除后不能重新计算法线、也不能留下索引
    指向已删除的面；但本函数**不删除顶点**（未引用顶点由 :func:`compact`
    负责），所以 ``positions``/``uvs``/``normals`` 原样保留。

    参数
    ----
    mesh : dict
        见模块说明。

    返回
    ----
    dict
        ``{"positions": [...], "uvs": [...], "normals": [...],
        "indices": [...]}``。顶点数组与输入等长、同序，只有 ``indices`` 被
        过滤。空网格返回四个空列表。

    异常
    ----
    WeldError
        网格字段缺失/长度不一致、下标越界或存在非有限坐标。
    """
    positions, uvs, normals, triangles = _check_mesh(mesh)
    kept = [list(triangle) for triangle in triangles
            if _triangle_area(positions, triangle) > EPS]
    return {
        "positions": _copy_vectors(positions),
        "uvs": _copy_vectors(uvs),
        "normals": _copy_vectors(normals),
        "indices": kept,
    }


def compact(mesh):
    """删除没有被任何三角形引用的顶点，并重映射索引。

    参数
    ----
    mesh : dict
        见模块说明。

    返回
    ----
    dict
        ``{"positions": [...], "uvs": [...], "normals": [...],
        "indices": [...], "remap": [...]}``。

        * 只保留下标在 ``indices`` 中出现过的顶点，并保持它们第一次被引用
          的先后顺序（按旧下标升序更便于复现）；
        * ``indices`` 必须改写为新下标；
        * ``remap`` 与输入顶点等长，保留顶点映射到新下标，删除顶点映射为
          ``-1``；
        * 没有三角形时所有顶点都未引用，顶点数组和 ``indices`` 为空，
          ``remap`` 全为 ``-1``。

    异常
    ----
    WeldError
        网格字段缺失/长度不一致、下标越界或存在非有限坐标。
    """
    positions, uvs, normals, triangles = _check_mesh(mesh)

    referenced = sorted({index for triangle in triangles for index in triangle})
    old_to_new = {old: new for new, old in enumerate(referenced)}
    remap = [old_to_new.get(old, -1) for old in range(len(positions))]

    indices = [[old_to_new[index] for index in triangle]
               for triangle in triangles]
    return {
        "positions": [list(positions[old]) for old in referenced],
        "uvs": [list(uvs[old]) for old in referenced],
        "normals": [list(normals[old]) for old in referenced],
        "indices": indices,
        "remap": remap,
    }


def merge_uv_islands(mesh, epsilon):
    """在同一个 UV 岛内部焊接顶点，跨岛绝不合并。

    UV 岛判定只看三角形拓扑：两个三角形共享一条边（一条边由同一对顶点索引
    构成）就属于同一个岛；不共享边的三角形属于不同岛。**不能用 UV 包围盒
    相交判断岛**：两个拓扑分离的岛即使 UV 包围盒完全重叠，也不能合并。

    岛内合并键是 ``(位置量化键, UV 容差分组)``，即只有位置在 ``epsilon``
    量化下相同、UV 在 ``1e-6`` 内相同的重复顶点才焊接。不同岛即使位置和 UV
    都相同也必须保持为不同顶点，否则跨岛数据会串在一起。

    参数
    ----
    mesh : dict
        见模块说明。
    epsilon : float
        位置量化步长，必须大于 0。

    返回
    ----
    dict
        ``{"positions": [...], "uvs": [...], "indices": [...],
        "remap": [...]}``。代表顶点按输入顺序第一次出现；``indices`` 按
        ``remap`` 改写；``remap`` 与输入顶点等长。空网格返回空字段。

    异常
    ----
    WeldError
        网格字段缺失/长度不一致、``epsilon <= 0``、下标越界、UV 数量不匹配
        或坐标非有限。
    """
    step = _check_epsilon(epsilon)
    positions, uvs, _normals, triangles = _check_mesh(mesh)

    # 并查集：共享一条无向边的三角形属于同一个 UV 岛。
    parent = list(range(len(triangles)))

    def find(item):
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    def union(left, right):
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[max(left_root, right_root)] = min(left_root, right_root)

    edge_owner = {}
    for triangle_index, triangle in enumerate(triangles):
        a, b, c = triangle
        for edge in ((a, b), (b, c), (c, a)):
            edge_key = (min(edge), max(edge))
            if edge_key in edge_owner:
                union(triangle_index, edge_owner[edge_key])
            else:
                edge_owner[edge_key] = triangle_index

    # 每个顶点所属岛的签名。只由顶点相连的三角形（可能与多个岛相连）的根组成，
    # 因此仅共享顶点的两个岛不会混为一谈。
    vertex_islands = [set() for _ in positions]
    for triangle_index, triangle in enumerate(triangles):
        root = find(triangle_index)
        for vertex in triangle:
            vertex_islands[vertex].add(root)

    groups = {}
    merged_positions = []
    merged_uvs = []
    remap = []
    for old_index, position in enumerate(positions):
        islands = vertex_islands[old_index]
        if islands:
            signature = ("islands", tuple(sorted(islands)))
        else:
            # 未被任何三角形引用的顶点不属于任何岛，各自独立处理。
            signature = ("unreferenced", old_index)
        key = (signature, _position_key(position, step))
        u, v = float(uvs[old_index][0]), float(uvs[old_index][1])
        candidates = groups.setdefault(key, [])
        new_index = _find_uv_group(candidates, u, v)
        if new_index < 0:
            new_index = len(merged_positions)
            candidates.append((u, v, new_index))
            merged_positions.append(list(position))
            merged_uvs.append([u, v])
        remap.append(new_index)

    indices = [[remap[index] for index in triangle] for triangle in triangles]
    return {
        "positions": merged_positions,
        "uvs": merged_uvs,
        "indices": indices,
        "remap": remap,
    }


def validate(mesh):
    """检查网格问题并返回人类可读的字符串列表；无问题返回 ``[]``。

    这是诊断函数：遇到下列几何/结构问题时**不要抛异常**，而是把每个问题
    收集到返回列表。至少检查：

    * ``indices`` 中下标越界（问题字符串包含 ``index`` 和具体下标）；
    * ``positions``/``uvs``/``normals`` 中存在 NaN（问题字符串包含 ``nan``）
      或无穷大；
    * 三角形面积 ``<= EPS`` 的退化面（问题字符串包含 ``degenerate``）；
    * 缺少字段、数组长度不一致、``indices`` 元素不是三元组。

    返回列表的顺序应稳定（按字段和下标升序），便于测试复现。输入网格本身
    不能被修改。

    参数
    ----
    mesh : dict
        见模块说明。

    返回
    ----
    list[str]
        问题描述列表；没有问题时为 ``[]``。

    异常
    ----
    无。只要还能遍历输入就返回问题列表；无法遍历的极端输入可抛
    :class:`WeldError`。
    """
    problems = []
    if not isinstance(mesh, dict):
        return ["mesh is not a dict"]

    arrays = {}
    for name in _FIELDS:
        if name not in mesh:
            problems.append("mesh is missing field %s" % name)
        elif not isinstance(mesh[name], (list, tuple)):
            problems.append("mesh field %s is not a list" % name)
        else:
            arrays[name] = list(mesh[name])

    vertex_count = len(arrays["positions"]) if "positions" in arrays else None
    for name in ("uvs", "normals"):
        if name in arrays and vertex_count is not None:
            if len(arrays[name]) != vertex_count:
                problems.append(
                    "%s length %d does not match positions length %d"
                    % (name, len(arrays[name]), vertex_count))

    def component_issues(vector, size, label):
        """返回 (问题列表, 有限 float 列表或 None)。"""
        if not isinstance(vector, (list, tuple)) or len(vector) != size:
            return ["%s must have %d components" % (label, size)], None
        issues = []
        values = []
        for k, component in enumerate(vector):
            try:
                number = float(component)
            except (TypeError, ValueError):
                issues.append("%s[%d] is not a number" % (label, k))
                values.append(None)
                continue
            values.append(number)
            if math.isnan(number):
                issues.append("%s[%d] contains nan" % (label, k))
            elif math.isinf(number):
                issues.append("%s[%d] contains infinite value" % (label, k))
        if any(value is None or not math.isfinite(value) for value in values):
            return issues, None
        return issues, values

    position_values = []
    if "positions" in arrays:
        for i, vector in enumerate(arrays["positions"]):
            issues, values = component_issues(vector, 3, "positions[%d]" % i)
            problems.extend(issues)
            position_values.append(values)
    if "uvs" in arrays:
        for i, vector in enumerate(arrays["uvs"]):
            issues, _values = component_issues(vector, 2, "uvs[%d]" % i)
            problems.extend(issues)
    if "normals" in arrays:
        for i, vector in enumerate(arrays["normals"]):
            issues, _values = component_issues(vector, 3, "normals[%d]" % i)
            problems.extend(issues)

    valid_triangles = []
    if "indices" in arrays:
        for triangle_index, triangle in enumerate(arrays["indices"]):
            if not isinstance(triangle, (list, tuple)) or len(triangle) != 3:
                problems.append("indices[%d] is not a triangle" % triangle_index)
                continue
            parsed = []
            usable = True
            for corner, index in enumerate(triangle):
                if not isinstance(index, int):
                    problems.append("indices[%d][%d] is not an integer index"
                                    % (triangle_index, corner))
                    usable = False
                    continue
                if vertex_count is not None and not (0 <= index < vertex_count):
                    problems.append(
                        "indices[%d][%d] index %d out of range of %d vertices"
                        % (triangle_index, corner, index, vertex_count))
                    usable = False
                    continue
                parsed.append(index)
            if usable:
                valid_triangles.append((triangle_index, tuple(parsed)))

    for triangle_index, triangle in valid_triangles:
        corners = []
        usable = True
        for vertex in triangle:
            if vertex >= len(position_values) or position_values[vertex] is None:
                usable = False
                break
            corners.append(position_values[vertex])
        if not usable:
            continue
        area = _triangle_area(corners, (0, 1, 2))
        if area <= EPS:
            problems.append("indices[%d] is degenerate (area %.6g <= %.6g)"
                            % (triangle_index, area, EPS))

    return problems


def weld_batch(meshes, config):
    """批量焊接多个网格；每个网格使用独立的 remap 表。

    参数
    ----
    meshes : sequence[dict]
        待焊接网格，见模块说明。
    config : dict
        ``{"epsilon": float}`` 是必需项，``epsilon`` 必须为有限正数；
        可选 ``"preserve_seams": bool``，默认 ``False``。为 ``True`` 时对每个
        网格调用 :func:`preserve_seams`，否则调用 :func:`weld_positions`。
        其他未知键忽略。

    返回
    ----
    dict
        ``{"meshes": [result0, result1, ...], "remaps": [remap0, remap1, ...]}``，
        严格保持输入网格顺序。每个 ``result`` 是该网格的焊接结果，``remaps[i]``
        必须与 ``meshes[i]`` 的输入顶点一一对应；**不同网格不得共用同一个
        remap 字典/列表**，否则第二个网格的索引会串到第一个网格上。

    异常
    ----
    WeldError
        ``config`` 缺 ``epsilon``、``epsilon <= 0``、任一网格字段非法、下标
        越界或坐标非有限。实现必须先验证再计算：任一网格失败时整批失败，
        不返回部分结果，也不得原地修改传入的网格（整体回滚）。
    """
    if not isinstance(config, dict):
        raise WeldError("config must be a dict")
    if "epsilon" not in config:
        raise WeldError("config is missing 'epsilon'")
    step = _check_epsilon(config["epsilon"])
    use_seams = bool(config.get("preserve_seams", False))

    if isinstance(meshes, (list, tuple)):
        batch = list(meshes)
    else:
        try:
            batch = list(meshes)
        except TypeError:
            raise WeldError("meshes must be a sequence of meshes")

    # 先验证整批，任一网格非法就整体失败，不产生部分结果。
    for mesh in batch:
        _check_mesh(mesh)

    results = []
    remaps = []
    for mesh in batch:
        if use_seams:
            result = preserve_seams(mesh, step)
        else:
            result = weld_positions(mesh, step)
        results.append(result)
        # 每个网格各自持有一份 remap，绝不与其他网格共享。
        remaps.append(list(result["remap"]))

    return {"meshes": results, "remaps": remaps}
