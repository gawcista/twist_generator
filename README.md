# TwistGenerator

Language: **English** | [中文](#中文文档)

TwistGenerator is a small Python toolkit for building twisted bilayer and aligned bilayer slab structures from VASP-style POSCAR/CONTCAR inputs. The default documentation is shown in English; the Chinese version is available at the bottom of this page.

## Core Distance Convention

The bilayer spacing parameter `d` is the target Cartesian distance, in Angstrom, between the closest inner surfaces of the two atomic slabs.

The final c-axis length is built from the physical slab thicknesses, the target interlayer gap, and the external vacuum already present in the input structures:

```text
vacuum_A = c_A - thickness_A
vacuum_B = c_B - thickness_B
vacuum_final = max(vacuum_A, vacuum_B)

c_final = thickness_A + d + thickness_B + vacuum_final
```

This means:

- `d` is not the final c-axis length.
- `d` is not the center-to-center distance between layers.
- By default, the final external vacuum is the larger external vacuum from the two input slabs.
- Pass `vacuum=<number>` to override the final external vacuum explicitly.

The code identifies the slab thickness with the largest periodic z-gap method. This is more robust for slab structures whose atoms are wrapped across the fractional z boundary.

## Basic Setup

Load a single-layer POSCAR:

```python
from twist_generator.structure import atomic_structure

lattice = atomic_structure("POSCAR")
```

If Python imports an older installed package instead of this checkout, run from the parent directory or set:

```bash
export PYTHONPATH=/home/gawcista/scripts
```

## `generator` Module

Use `twist_generator.generator` when you want ready-made structure-generation helpers that write output files directly.

### Hexagonal 210 Twist

```python
from twist_generator.structure import atomic_structure
from twist_generator.generator import generate_hexagonal_210

lattice = atomic_structure("POSCAR")

generate_hexagonal_210(
    m=7,
    n=3,
    lattice=lattice,
    d=3.35,
    filename="twisted_hex_210.vasp",
)
```

### Square Twists

```python
from twist_generator.generator import generate_square_100, generate_square_110

generate_square_100(
    m=5,
    n=2,
    lattice=lattice,
    d=3.35,
    filename="twisted_square_100.vasp",
)

generate_square_110(
    m=5,
    n=2,
    lattice=lattice,
    d=3.35,
    filename="twisted_square_110.vasp",
)
```

### Aligned Heterobilayer

```python
from twist_generator.structure import atomic_structure
from twist_generator.generator import generate_aligned

lattice_a = atomic_structure("POSCAR_A")
lattice_b = atomic_structure("POSCAR_B")

generate_aligned(
    latticeA=lattice_a,
    latticeB=lattice_b,
    m=4,
    n=5,
    d=3.35,
    vacuum=18.0,
    translation=[0.0, 0.0, 0.0],
    filename="aligned_4_5.vasp",
)
```

Available helpers:

```python
generate_bilayer(P1, P2, lattice, d=4, filename='CONTCAR', vacuum=None, **kargs)
generate_square_100(m, n, **kargs)
generate_square_110(m, n, **kargs)
generate_hexagonal_210(m, n, **kargs)
generate_aligned(latticeA, latticeB, m, n, d=4, filename='CONTCAR', vacuum=None, **kargs)
```

## `build` Module

Use `twist_generator.build` when you want the generated structure object back instead of writing a file immediately.

### Twisted Bilayer From Explicit Matrices

```python
from twist_generator.build import build_bilayer

P1 = [[7, -3, 0], [3, 4, 0], [0, 0, 1]]
P2 = [[7, -4, 0], [4, 3, 0], [0, 0, 1]]

bilayer = build_bilayer(
    P1=P1,
    P2=P2,
    d=3.35,
    lattice=lattice,
    use_axis='A',
    translation=[0.0, 0.0, 0.0],
    vacuum=None,
)

bilayer.print_POSCAR("CONTCAR")
```

`use_axis='A'` keeps layer A's in-plane basis and strains/rebases layer B onto it. `use_axis='B'` does the opposite.

### Aligned Bilayer From Two Materials

```python
from twist_generator.build import build_bilayer_aligned

bilayer = build_bilayer_aligned(
    P1=[[4, 0, 0], [0, 4, 0], [0, 0, 1]],
    P2=[[5, 0, 0], [0, 5, 0], [0, 0, 1]],
    d=3.35,
    lattA=lattice_a,
    lattB=lattice_b,
    use_axis='A',
    translation=[0.25, 0.0, 0.0],
    vacuum=18.0,
)

bilayer.print_POSCAR("aligned.vasp")
```

Available builders:

```python
build_bilayer(P1, P2, d, lattice, use_axis='A', translation=[0., 0., 0.], vacuum=None)
build_bilayer_aligned(P1, P2, d, lattA, lattB, use_axis='A', translation=[0., 0., 0.], vacuum=None)
build_trilayer_ABC_hex(P1, P2, P3, d, lattA, lattB, lattC, autocorrection=True)
```

`build_trilayer_ABC_hex` is present for legacy trilayer workflows. The current spacing tests focus on bilayers.

## `search` Module

Use `twist_generator.search` to print candidate commensurate cells before generating a structure. The search functions print tab-separated rows. For twist searches, the columns are generally:

```text
m    n    angle(deg)    atom_count
```

For aligned searches, the columns are:

```text
m    n    supercell_length    strain_percent
```

### Search Hexagonal Candidates

```python
from twist_generator.structure import atomic_structure
from twist_generator.search import search_hexagonal_210

lattice = atomic_structure("POSCAR")

search_hexagonal_210(
    Mcut=20,
    lattice=lattice,
    writeflag=False,
)
```

After choosing a printed `(m, n)` pair, generate the file with `generator`:

```python
from twist_generator.generator import generate_hexagonal_210

generate_hexagonal_210(
    m=13,
    n=6,
    lattice=lattice,
    d=3.35,
    filename="hex_13_6.vasp",
)
```

### Search Square Candidates

```python
from twist_generator.search import search_square_100, search_square_210

search_square_100(Mcut=20, lattice=lattice, writeflag=False)
search_square_210(Mcut=20, lattice=lattice, writeflag=False)
```

### Search Aligned Heterobilayers

```python
from twist_generator.search import search_aligned

search_aligned(
    Mcut=30,
    latticeA=lattice_a,
    latticeB=lattice_b,
    strain_tol=0.001,
    writeflag=False,
)
```

For the most predictable workflow, use `writeflag=False`, select a candidate from the printed table, then call the matching `generator` helper manually.

## Input Requirements

The spacing reconstruction assumes the input structures are 2D slabs:

- The c axis is aligned with z, approximately `[0, 0, c]`.
- The a and b lattice vectors do not contain z components.
- There is a unique and clear vacuum region along c.
- The input is not a 3D bulk-like periodic structure.

If these requirements are not met, `combine_lattice()` raises `ValueError` instead of silently producing an incorrect structure.

## Output Checks

A valid generated bilayer should satisfy:

```text
layer_A_top to layer_B_bottom = d
periodic external vacuum = vacuum_final
```

If the spacing looks unusual, check:

- whether the input POSCAR c-axis length is reasonable;
- whether the single-layer structure has a unique vacuum region;
- whether `d` was meant as surface spacing, not center spacing;
- whether you should pass `vacuum=<number>` explicitly.

## Tests

Run the unit tests:

```bash
python -m unittest discover -s tests -v
```

Or force this checkout onto `PYTHONPATH`:

```bash
PYTHONPATH=/home/gawcista/scripts python -m unittest discover -s tests -v
```

The tests cover:

- toy slab spacing and external vacuum;
- heterobilayer vacuum selection;
- wrapped fractional z coordinates;
- explicit `vacuum` override;
- aligned bilayer x/y translation;
- invalid negative parameters, tilted c axis, and bulk-like z distributions.

---

<a id="中文文档"></a>

[English](#twistgenerator) | **中文**

# TwistGenerator

TwistGenerator 是一个用于从 VASP POSCAR/CONTCAR 输入生成二维材料转角双层和 aligned 双层 slab 结构的 Python 工具。默认文档显示英文；中文版本在此处展开查看。

## 核心层间距约定

双层结构的层间距参数 `d` 表示上下两层原子 slab 最近表面之间的目标距离，单位为 Angstrom。

最终双层晶胞的 c 轴长度由两层真实厚度、目标层间距和输入结构原本的外部真空共同决定：

```text
vacuum_A = c_A - thickness_A
vacuum_B = c_B - thickness_B
vacuum_final = max(vacuum_A, vacuum_B)

c_final = thickness_A + d + thickness_B + vacuum_final
```

这意味着：

- `d` 不是最终 c 轴长度。
- `d` 不是两层几何中心之间的距离。
- 默认最终外部真空会取两层输入 slab 中较大的原始外部真空。
- 可以传入 `vacuum=<number>` 显式覆盖最终外部真空。

代码使用最大 periodic z-gap 方法识别 slab 厚度和外部真空。对于原子跨过 fractional z 边界的 slab 结构，这比简单使用 `z_max - z_min` 更稳健。

## 基本设置

读取单层 POSCAR：

```python
from twist_generator.structure import atomic_structure

lattice = atomic_structure("POSCAR")
```

如果 Python 导入了已安装的旧版本，而不是当前源码目录，可以从父目录运行，或设置：

```bash
export PYTHONPATH=/home/gawcista/scripts
```

## `generator` 模块

当你希望直接生成并写出结构文件时，使用 `twist_generator.generator`。

### 六角晶格 210 转角

```python
from twist_generator.structure import atomic_structure
from twist_generator.generator import generate_hexagonal_210

lattice = atomic_structure("POSCAR")

generate_hexagonal_210(
    m=7,
    n=3,
    lattice=lattice,
    d=3.35,
    filename="twisted_hex_210.vasp",
)
```

### 方形晶格转角

```python
from twist_generator.generator import generate_square_100, generate_square_110

generate_square_100(
    m=5,
    n=2,
    lattice=lattice,
    d=3.35,
    filename="twisted_square_100.vasp",
)

generate_square_110(
    m=5,
    n=2,
    lattice=lattice,
    d=3.35,
    filename="twisted_square_110.vasp",
)
```

### Aligned 异质双层

```python
from twist_generator.structure import atomic_structure
from twist_generator.generator import generate_aligned

lattice_a = atomic_structure("POSCAR_A")
lattice_b = atomic_structure("POSCAR_B")

generate_aligned(
    latticeA=lattice_a,
    latticeB=lattice_b,
    m=4,
    n=5,
    d=3.35,
    vacuum=18.0,
    translation=[0.0, 0.0, 0.0],
    filename="aligned_4_5.vasp",
)
```

可用 helper：

```python
generate_bilayer(P1, P2, lattice, d=4, filename='CONTCAR', vacuum=None, **kargs)
generate_square_100(m, n, **kargs)
generate_square_110(m, n, **kargs)
generate_hexagonal_210(m, n, **kargs)
generate_aligned(latticeA, latticeB, m, n, d=4, filename='CONTCAR', vacuum=None, **kargs)
```

## `build` 模块

当你希望拿到结构对象而不是立即写文件时，使用 `twist_generator.build`。

### 用显式变换矩阵构造转角双层

```python
from twist_generator.build import build_bilayer

P1 = [[7, -3, 0], [3, 4, 0], [0, 0, 1]]
P2 = [[7, -4, 0], [4, 3, 0], [0, 0, 1]]

bilayer = build_bilayer(
    P1=P1,
    P2=P2,
    d=3.35,
    lattice=lattice,
    use_axis='A',
    translation=[0.0, 0.0, 0.0],
    vacuum=None,
)

bilayer.print_POSCAR("CONTCAR")
```

`use_axis='A'` 表示保留 A 层的面内晶格基矢，并将 B 层匹配到 A 层；`use_axis='B'` 则相反。

### 两种材料的 aligned 双层

```python
from twist_generator.build import build_bilayer_aligned

bilayer = build_bilayer_aligned(
    P1=[[4, 0, 0], [0, 4, 0], [0, 0, 1]],
    P2=[[5, 0, 0], [0, 5, 0], [0, 0, 1]],
    d=3.35,
    lattA=lattice_a,
    lattB=lattice_b,
    use_axis='A',
    translation=[0.25, 0.0, 0.0],
    vacuum=18.0,
)

bilayer.print_POSCAR("aligned.vasp")
```

可用 builder：

```python
build_bilayer(P1, P2, d, lattice, use_axis='A', translation=[0., 0., 0.], vacuum=None)
build_bilayer_aligned(P1, P2, d, lattA, lattB, use_axis='A', translation=[0., 0., 0.], vacuum=None)
build_trilayer_ABC_hex(P1, P2, P3, d, lattA, lattB, lattC, autocorrection=True)
```

`build_trilayer_ABC_hex` 是保留的三层结构 workflow；当前层间距测试主要覆盖双层结构。

## `search` 模块

使用 `twist_generator.search` 可以在真正生成结构之前打印 commensurate cell 候选。转角搜索通常输出：

```text
m    n    angle(deg)    atom_count
```

aligned 搜索输出：

```text
m    n    supercell_length    strain_percent
```

### 搜索六角晶格候选

```python
from twist_generator.structure import atomic_structure
from twist_generator.search import search_hexagonal_210

lattice = atomic_structure("POSCAR")

search_hexagonal_210(
    Mcut=20,
    lattice=lattice,
    writeflag=False,
)
```

选定打印出的 `(m, n)` 后，再用 `generator` 生成结构：

```python
from twist_generator.generator import generate_hexagonal_210

generate_hexagonal_210(
    m=13,
    n=6,
    lattice=lattice,
    d=3.35,
    filename="hex_13_6.vasp",
)
```

### 搜索方形晶格候选

```python
from twist_generator.search import search_square_100, search_square_210

search_square_100(Mcut=20, lattice=lattice, writeflag=False)
search_square_210(Mcut=20, lattice=lattice, writeflag=False)
```

### 搜索 aligned 异质双层

```python
from twist_generator.search import search_aligned

search_aligned(
    Mcut=30,
    latticeA=lattice_a,
    latticeB=lattice_b,
    strain_tol=0.001,
    writeflag=False,
)
```

最稳妥的流程是使用 `writeflag=False` 只打印候选，选定后再手动调用对应的 `generator` helper。

## 输入结构要求

当前层间距重建逻辑假设输入是二维 slab：

- c 轴沿 z 方向，近似 `[0, 0, c]`。
- a/b 晶格矢量没有 z 分量。
- c 方向存在唯一且清晰的外部真空区域。
- 输入不是 3D bulk-like 周期结构。

如果不满足这些条件，`combine_lattice()` 会抛出 `ValueError`，避免静默生成错误结构。

## 输出结构检查

生成的双层结构应满足：

```text
layer_A_top 到 layer_B_bottom 的距离 = d
周期边界方向的外部真空 = vacuum_final
```

如果层间距看起来不正常，优先检查：

- 输入 POSCAR 的 c 轴长度是否合理；
- 单层结构是否有唯一的真空区；
- `d` 是否被当作表面间距，而不是中心距；
- 是否需要显式传入 `vacuum=<number>`。

## 测试

运行单元测试：

```bash
python -m unittest discover -s tests -v
```

或强制使用当前源码目录：

```bash
PYTHONPATH=/home/gawcista/scripts python -m unittest discover -s tests -v
```

测试覆盖：

- toy slab 的层间距和外部真空；
- 异质双层自动选择较大外部真空；
- wrapped fractional z 坐标；
- 显式 `vacuum` 覆盖；
- aligned 双层 x/y 平移；
- 负参数、倾斜 c 轴和 bulk-like z 分布等异常输入。
