# [NanoStudio: A Python Package for Automated Construction of Twisted Two-Dimensional Materials and Multi-Walled Nanotubes]()

Here, we present NanoStudio, an open-source Python package for the automated construction of low-dimensional nanostructures based on honeycomb lattices. The package provides a unified implementation for generating graphene with arbitrary chirality, multilayer graphene with user-defined stacking sequences, twisted graphene with arbitrary commensurate twist angles, quasi-one-dimensional nanoribbons, and multi-walled nanotubes. 

![alt text](image.png)

## Authors

This package was primarily written by Chenglong Qin (clqin@xhu.edu.cn).

## License

The model is released under the MIT License.

# Parameter Documentation for Nanostructure Generation Code
This object-oriented code generates **graphene, twisted bilayer graphene (TBG), and carbon nanotubes (CNTs)** based on the ASE library. Below is a full English summary of classes, parameters, methods and usage examples.

## 1. Base Class: `Nano`
Abstract base class providing common parameters, coordinate transformation, mathematical tools and basic structure functions for all carbon nanostructures.

### 1.1 Constructor `__init__()`
| Parameter | Type | Description | Unit & Note |
| ---- | ---- | ---- | ---- |
| `n` | int | Chiral index *n* | Standard index for carbon nanostructures |
| `m` | int | Chiral index *m* | Standard index for carbon nanostructures |
| `period` | list | Supercell repetition `[x, y, z]` | Lattice replication along three axes |
| `pbc` | list[bool] | Periodic boundary conditions `[x, y, z]` | `True`: enable PBC; `False`: add vacuum layer |
| `bond_length` | float | C-C bond length | Unit: Å, default = 1.42 Å |
| `atom_type` | list | Element symbol of atoms | e.g. `['C']`, `['C', 'B']` |
| `vacuum` | float | Vacuum layer thickness | Unit: Å, applied to non-periodic directions |

### 1.2 General Methods
- `rotation_a1a2(theta)`: 2D rotation matrix under oblique lattice basis **(a1, a2)**. Input `theta` in **radians**.
- `a1a2_to_ij(theta)`: Transformation matrix from oblique coordinates to Cartesian coordinates.
- `rotation_ij(theta)`: Standard 2D rotation matrix in Cartesian system. Input `theta` in **radians**.
- `extended_gcd(a, b)`: Extended Euclidean algorithm. Solve $ax+by=\gcd(a,b)$, returns `(gcd, x, y)`.
- `_get_structures()`: Generate ASE `Atoms` object with lattice, atomic positions, PBC and vacuum settings.
- `_get_properties()`: Calculate structural properties including chiral angle, bond length and PBC.
- `show_properties()`: Print structural parameters in console.

---

## 2. Subclass: `Graphene(Nano)`
For constructing **single-layer and multi-layer stacked graphene** (AA / AB / arbitrary stacking).

### 2.1 Constructor `__init__()`
| Parameter | Type | Description | Default Value |
| ---- | ---- | ---- | ---- |
| `n` | int | Chiral index *n* | — |
| `m` | int | Chiral index *m* | — |
| `p` | int | Component *p* of translation vector **T** | — |
| `q` | int | Component *q* of translation vector **T** | — |
| `xy_period` | list | Supercell replication along x & y `[nx, ny]` | `[1, 1]` |
| `xy_pbc` | list[bool] | PBC for x & y directions | `[True, True]` |
| `bond_length` | float | C-C bond length | `1.42` Å |
| `atom_type` | list | Atomic element | `['C']` |
| `vacuum` | float | Vacuum layer thickness | `7.5` Å |

> Internal logic: Z direction is automatically set to non-periodic with no replication.

### 2.2 Main Method `get_structures()`
| Parameter | Type | Description | Restriction |
| ---- | ---- | ---- | ---- |
| `layer_num` | int | Number of graphene layers | Must match length of `translations` |
| `translations` | list[tuple] | In-plane offset for each layer (fractional coordinates) | Format: `[(x1,y1), (x2,y2), ...]` |
| `layer_spacing` | float | Vertical interlayer distance | Unit: Å, typical value = 3.4 Å |

### 2.3 Private Methods
- `_build_coords()`: Generate fractional coordinates for single-layer graphene (A / B sublattice).
- `_get_properties()`: Calculate lattice matrix, total atom number, chiral angle and translation vectors.

---

## 3. Subclass: `TwistGraphene(Nano)`
For constructing **twisted multi-layer graphene (TBG)** with rotational misalignment between layers.

### 3.1 Constructor `__init__()`
| Parameter | Type | Description | Default Value |
| ---- | ---- | ---- | ---- |
| `n` | int | Chiral index *n* | — |
| `m` | int | Chiral index *m* | — |
| `xy_period` | tuple | Supercell replication `(nx, ny)` | `(1, 1)` |
| `xy_pbc` | list[bool] | PBC for x & y directions | `[True, True]` |
| `bond_length` | float | C-C bond length | `1.42` Å |
| `atom_type` | list | Atomic element | `['C']` |
| `vacuum` | float | Vacuum layer thickness | `7.5` Å |

### 3.2 Main Method `get_structures()`
| Parameter | Type | Description | Restriction |
| ---- | ---- | ---- | ---- |
| `layer_num` | int | Total number of layers | Length must match `translations` & `layer_sequence` |
| `translations` | list[tuple] | In-plane offset for each layer | Fractional coordinates |
| `layer_sequence` | list[int] | Layer type flag | `0`: original layer; `1`: rotated layer |
| `layer_spacing` | float | Vertical interlayer distance | Unit: Å |

### 3.3 Notes
- Twist angle is automatically calculated from `(n, m)`.
- Layer rotation is realized via Cartesian rotation matrix.
- `_get_single_layer()`: Generate coordinates for original or rotated graphene monolayer.

---

## 4. Subclass: `Nanotube(Nano)`
For constructing **single-walled and multi-walled carbon nanotubes (CNTs)**; supports multiple tubes in one model.

### 4.1 Constructor `__init__()`
| Parameter | Type | Description | Default Value |
| ---- | ---- | ---- | ---- |
| `n` | list[int] | List of chiral index *n* for each tube | e.g. `[15, 6]` for two tubes |
| `m` | list[int] | List of chiral index *m* for each tube | Same length as `n` |
| `z_period` | list[int] | Axial replication along tube axis (z) | One value per tube |
| `z_pbc` | bool | PBC along z-axis (tube axis) | `True` |
| `bond_length` | float | C-C bond length | `1.42` Å |
| `atom_type` | list | Atomic element | `['C']` |
| `vacuum` | float | Radial vacuum thickness | `7.5` Å |

> Rule: PBC is disabled for radial x & y directions by default.

### 4.2 Main Method `get_structures()`
| Parameter | Type | Description |
| ---- | ---- | ---- |
| `translation` | list[tuple] | In-plane offset for each nanotube | Length equals total tube number |

### 4.3 Calculated Structural Properties
- `chirality (n, m)`: Chiral indices
- `diameter`: Tube diameter (Å)
- `length`: Axial period length (Å)
- `chiral_theta`: Chiral angle (°)
- `total_atom_num`: Total number of atoms

### 4.4 Coordinate Logic
- 2D graphene plane is rolled into cylindrical coordinates to form CNT geometry.
- Atomic positions for A / B sublattices are calculated via rotation angle and axial shift.

---

## 5. Default Global Parameters
| Physical Quantity | Default Value | Remarks |
| ---- | ---- | ---- |
| C-C bond length | 1.42 Å | Standard value for carbon materials |
| Vacuum layer | 7.5 Å | Used for all non-periodic directions |
| Graphene interlayer distance | 3.4 Å | Typical van der Waals gap for graphite |
| Default atom type | `['C']` | Pure carbon structure |

---

## 6. Code Usage Examples
### 6.1 Multilayer Graphene
```python
model = Graphene(4, 2, -2, 6, xy_period=[1, 1], xy_pbc=[True, True])
struct = model.get_structures(layer_num=3, translations=[(0, 0), (1/3, 1/3), (-1/3, 2/3)])
```
- Chirality: $n=4, m=2$; Translation vector: $p=-2, q=6$
- 3 layers with different in-plane shifts for varied stacking configurations

### 6.2 Twisted Graphene
```python
model = TwistGraphene(4, 2, xy_period=[2, 2], xy_pbc=[True, True])
struct = model.get_structures(layer_num=3, translations=[(0, 0), (0, 0), (0, 0)], layer_sequence=[1, 0, 1])
```
- Chirality: $n=4, m=2$; 2×2 supercell in x-y plane
- 3 layers: rotated layer → original layer → rotated layer

### 6.3 Multi-Wall Carbon Nanotubes
```python
model = Nanotube([15, 6], [15, 6], z_period=[5, 5], atom_type=['C'], z_pbc=False)
struct = model.get_structures(translation=[(0, 0), (0, 0)])
```
- Two armchair CNTs: $(15,15)$ and $(6,6)$
- 5 times axial replication; PBC disabled along tube axis; no in-plane offset





























表格Physical QuantityDefault ValueRemarksC-C bond length1.42 ÅStandard value for carbon materialsVacuum layer7.5 ÅUsed for all non-periodic directionsGraphene interlayer distance3.4 ÅTypical van der Waals gap for graphiteDefault atom type['C']Pure carbon structure6. Code Usage Examples6.1 Multilayer Graphenepython运行model = Graphene(4, 2, -2, 6, xy_period=[1, 1], xy_pbc=[True, True])
struct = model.get_structures(layer_num=3, translations=[(0, 0), (1/3, 1/3), (-1/3, 2/3)])

Chirality: \(n=4, m=2\); Translation vector: \(p=-2, q=6\)
3 layers with different in-plane shifts for varied stacking configurations
6.2 Twisted Graphenepython运行model = TwistGraphene(4, 2, xy_period=[2, 2], xy_pbc=[True, True])
struct = model.get_structures(layer_num=3, translations=[(0, 0), (0, 0), (0, 0)], layer_sequence=[1, 0, 1])

Chirality: \(n=4, m=2\); 2×2 supercell in x-y plane
3 layers: rotated layer → original layer → rotated layer
6.3 Multi-Wall Carbon Nanotubespython运行model = Nanotube([15, 6], [15, 6], z_period=[5, 5], atom_type=['C'], z_pbc=False)
struct = model.get_structures(translation=[(0, 0), (0, 0)])

Two armchair CNTs: \((15,15)\) and \((6,6)\)
5 times axial replication; PBC disabled along tube axis; no in-plane offset