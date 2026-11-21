# mesh-uv-weld

## 题面

模型导入后做顶点焊接，UV 缝被焊死、贴图整体撕裂：位置量化用 round 而不是 floor(x/eps)，边界被错分；合并时不看 UV，UV 缝被焊死；合并后索引不重映射，模型撕裂；退化三角形不剔除，法线为 NaN；未引用顶点不删，文件变大；UV 岛判定用包围盒相交，误合并；索引越界不校验，越界读；NaN 坐标不报错，量化键为 NaN；批量共用 remap 表，网格之间串；epsilon<=0 不校验，除零。weld/vertices.py 是核心，tests/test_vertices.py 15个红。要修的是：按 docstring 写清的契约把 `weld/vertices.py` 里 8 个 `raise NotImplementedError` 的桩逐个实现，把上面 10 类缺陷全部修掉，让 `python -m pytest tests -q` 从 15 failed 变成 15 passed。

## 运行方式

```bash
cd mesh-uv-weld
python -m pytest tests -q
```

初始状态（未改动任何实现）必须恰好是：

```
15 failed
```

15 条用例全部红，且**没有** collection error / import error / syntax error——失败原因只能是
`NotImplementedError`。修完之后同样的命令应变成 `15 passed`。要求 Python 3，只装 pytest，
离线可跑，不使用随机数，重复运行结果一致。

## 字段表

| 字段 | 取值 |
| --- | --- |
| 所属领域 | 游戏、图形与多媒体 |
| 子方向 | 图形渲染 |
| 任务类型 | Bugfix |
| 语言/框架 | Python 3, pytest |
| 环境可复现等级 | 无外部依赖 |

## 约束

- **tests 不可修改。** `tests/test_vertices.py` 是验收标准，只能改 `weld/vertices.py`（以及
  `weld/__init__.py` 的导出）；不要为了让用例通过而改动、删除、跳过或重命名测试。
- 只用标准库（`math`）加 pytest，不引入 numpy 等第三方包；不使用随机数，结果可复现。
- 网格数据结构固定为 `{"positions": [[x, y, z], ...], "uvs": [[u, v], ...],
  "normals": [[nx, ny, nz], ...], "indices": [[i, j, k], ...]}`，四个字段按顶点下标对齐。
- 位置量化步长 `epsilon` 必须为正；位置量化用 `floor(x / epsilon)`，不能用 `round`。
- UV 缝不能焊死：位置相同但 UV 不同的顶点必须保留。UV 各分量和几何面积/长度的比较容差
  统一为 `1e-6`。
- 非法输入一律抛 `WeldError`；`validate` 例外，它把问题收集成字符串列表返回。
- 批量焊接中任一网格非法时整批失败、整体回滚，不能返回部分结果或原地改写输入网格。

## 文件结构

```
mesh-uv-weld/
├── README.md
├── conftest.py              # 把仓库根加入 sys.path，使 tests 可直接 import weld
├── .gitignore
├── weld/
│   ├── __init__.py
│   └── vertices.py          # 8 个桩 + WeldError，核心待实现文件
└── tests/
    └── test_vertices.py     # 15 条验收用例（不可修改）
```
