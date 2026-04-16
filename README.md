# TwistGenerator

TwistGenerator 是一个用于生成二维材料转角双层结构的 Python 工具。当前版本主要面向 VASP POSCAR/CONTCAR 输入，支持按给定变换矩阵构造上下两层 supercell，并输出新的 POSCAR 结构。

## 核心概念

双层结构的层间距参数 `d` 表示两层原子 slab 最近表面之间的距离，单位为 Angstrom。

最终双层晶胞的 c 轴长度由两层真实厚度、目标层间距和输入结构原本的外部真空共同决定：

```text
vacuum_A = c_A - thickness_A
vacuum_B = c_B - thickness_B
vacuum_final = max(vacuum_A, vacuum_B)

c_final = thickness_A + d + thickness_B + vacuum_final
```

也就是说：

- `d` 不是最终 c 轴长度。
- `d` 不是两层几何中心之间的距离。
- 默认外部真空会自动取两层输入结构中较大的原始真空。
- 如需手动控制最终外部真空，可传入 `vacuum=<number>`。

代码使用最大 periodic z-gap 方法识别 slab 厚度和外部真空。该方法适合有明确真空层的二维 slab 结构，也能处理原子坐标跨过 fractional z 边界的情况。

## 基本用法

从 POSCAR 读取单层结构：

```python
from twist_generator.structure import atomic_structure
from twist_generator.generator import generate_hexagonal_210

lattice = atomic_structure("POSCAR")

generate_hexagonal_210(
    m=7,
    n=3,
    lattice=lattice,
    d=3.35,
    filename="twisted_bilayer.vasp",
)
```

手动指定最终外部真空：

```python
generate_hexagonal_210(
    m=7,
    n=3,
    lattice=lattice,
    d=3.35,
    vacuum=18.0,
    filename="twisted_bilayer_vac18.vasp",
)
```

直接使用变换矩阵构造双层：

```python
from twist_generator.build import build_bilayer

P1 = [[7, -3, 0], [3, 4, 0], [0, 0, 1]]
P2 = [[7, -4, 0], [4, 3, 0], [0, 0, 1]]

bilayer = build_bilayer(P1, P2, d=3.35, lattice=lattice)
bilayer.print_POSCAR("CONTCAR")
```

## 常用 API

```python
combine_lattice(latticeA, latticeB, d, vacuum=None)
```

合并两个已经匹配好 a/b 晶格的 slab。`d` 为两层表面间距，`vacuum=None` 时自动取输入结构中较大的外部真空。

```python
build_bilayer(P1, P2, d, lattice, use_axis='A', translation=[0., 0., 0.], vacuum=None)
```

对同一个单层结构分别应用 `P1` 和 `P2`，构造转角双层。`translation` 作用于第二层的 fractional 坐标；z 方向层间距仍由 `d` 控制。

```python
build_bilayer_aligned(P1, P2, d, lattA, lattB, use_axis='A', translation=[0., 0., 0.], vacuum=None)
```

用于构造两种不同材料的 aligned 双层。第二层会应用 `translation`。

## 输入结构要求

当前层间距重建逻辑假设输入是二维 slab：

- c 轴基本沿 z 方向，即第三个晶格矢量形如 `[0, 0, c]`。
- a/b 晶格矢量没有 z 分量。
- c 方向存在唯一、清晰的外部真空区域。
- 输入不是 3D bulk-like 周期结构。

如果不满足这些条件，`combine_lattice()` 会抛出 `ValueError`，避免静默生成错误结构。

## 输出结构检查

构造完成后，双层结构应满足：

```text
layer_A_top 到 layer_B_bottom 的距离 = d
周期边界方向的外部真空 = vacuum_final
```

如果感觉层间距异常，优先确认：

- 输入 POSCAR 的 c 轴长度是否合理。
- 单层结构是否有足够且唯一的真空区。
- 你设置的 `d` 是否是表面间距，而不是中心距。
- 是否需要显式传入 `vacuum` 覆盖自动真空。

## 测试

运行单元测试：

```bash
python -m unittest discover -s tests -v
```

如果当前目录下的源码没有被优先导入，可从父目录运行，或临时设置：

```bash
PYTHONPATH=/home/gawcista/scripts python -m unittest discover -s tests -v
```

测试覆盖：

- 单层 toy slab 的表面间距和外部真空。
- 异质双层自动取较大外部真空。
- fractional z 跨边界的 slab。
- 手动 `vacuum` 覆盖。
- aligned 双层的 x/y 平移。
- 负参数、倾斜 c 轴、bulk-like 分布等异常输入。
