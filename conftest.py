"""把仓库根加入 sys.path，保证 `python -m pytest tests -q` 可直接 import weld。

不修改 tests/ 的内容；tests 只按路径发现用例。
"""

import os
import sys

_ROOT = os.path.dirname(os.path.abspath(__file__))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
