# LEGO Manual

Complete reference for LEGO v2.0 -- Lattice Elemental Geometry Operations.

## Table of Contents

- [Parameter File Format](#parameter-file-format)
- [Crystal Structures](#crystal-structures)
- [POSCAR Format](#poscar-format)
- [Rotation](#rotation)
- [Box Sizing](#box-sizing)
- [Cut Planes](#cut-planes)
- [Output Formats](#output-formats)
- [Python API](#python-api)
- [Building on macOS with OpenMP](#building-on-macos-with-openmp)
- [Migration from Old LEGO](#migration-from-old-lego)

---

## Building on macOS with OpenMP

The default Apple Clang compiler on macOS does not include OpenMP support.
LEGO will build and run without it (the atom counting pass runs single-threaded),
which is fine for systems up to a few million atoms.

For larger systems, install GCC via Homebrew to get OpenMP:

```bash
brew install gcc
make CC=gcc-15 OMP_FLAG="-fopenmp"
```

Replace `gcc-15` with whatever version Homebrew installed. You can check with:

```bash
ls /opt/homebrew/bin/gcc-*
```

To verify OpenMP is active, the program prints `OpenMP enabled: N threads available`
at startup.

---

## Parameter File Format

LEGO reads a plain-text parameter file where each line contains a keyword
followed by its value(s), separated by whitespace. Comments start with `#`.
Keywords are case-insensitive.

```
# This is a comment
structure    fcc
lattice_const 3.615
box_x        50
```

### Complete Keyword Reference

#### Crystal Definition

| Keyword | Arguments | Default | Description |
|---------|-----------|---------|-------------|
| `structure` | name | (required) | Built-in structure name (`sc`, `fcc`, `bcc`, `hex`, `hcp`, `L12`, `B2`, `D011`) or `poscar` for a custom file |
| `structure_file` | path | | Path to a POSCAR file. Required when `structure poscar` |
| `lattice_const` | float | from POSCAR | Override lattice parameter a (Angstrom). Scales the entire lattice uniformly relative to the first lattice vector length |
| `lattice_const_b` | float | from POSCAR | Override lattice parameter b. Only the second lattice vector is rescaled |
| `lattice_const_c` | float | from POSCAR | Override lattice parameter c. Only the third lattice vector is rescaled |
| `mass` | float | from POSCAR | Override atomic mass for type 0 (AMU) |
| `mass1` | float | from POSCAR | Override atomic mass for type 1 (AMU) |
| `atomtype` | int | -1 (no override) | Force all atoms to this type index |

#### Box Dimensions

| Keyword | Arguments | Default | Description |
|---------|-----------|---------|-------------|
| `box_x` | float | (required) | Target box dimension in x (Angstrom) |
| `box_y` | float | (required) | Target box dimension in y (Angstrom) |
| `box_z` | float | (required) | Target box dimension in z (Angstrom) |
| `box_mode` | string | `periodic` | `periodic` rounds box to the nearest integer multiple of the minimum repeat distance. `exact` uses the target dimensions as given |

#### Rotation

| Keyword | Arguments | Default | Description |
|---------|-----------|---------|-------------|
| `new_x` | 3 floats | | Crystallographic direction for the x-axis (Miller indices) |
| `new_y` | 3 floats | | Crystallographic direction for the y-axis |
| `new_z` | 3 floats | | Crystallographic direction for the z-axis |
| `axis1` | 3 floats | 1 0 0 | Axis for the first axis-angle rotation |
| `angle1` | float | 0 | Angle for the first rotation (degrees) |
| `axis2` | 3 floats | 1 0 0 | Axis for the second axis-angle rotation |
| `angle2` | float | 0 | Angle for the second rotation (degrees) |
| `rotmat1_x` | 3 floats | | First row of an explicit 3x3 rotation matrix |
| `rotmat1_y` | 3 floats | | Second row |
| `rotmat1_z` | 3 floats | | Third row |

#### Cut Planes

| Keyword | Arguments | Default | Description |
|---------|-----------|---------|-------------|
| `cut_norm1` | 3 floats | 0 0 0 | Normal vector of cut plane 1 (auto-normalized) |
| `plane_distance1` | float | | Signed distance from origin to plane 1: d = dot(point, normal) |
| `plane_point1` | 3 floats | | A point on plane 1 (alternative to `plane_distance1`) |
| `cutdir1` | float | 0 | `+1`: remove atoms with dot(pos, norm) >= d. `-1`: remove atoms with dot(pos, norm) <= d. `0`: no cutting |
| `cut_norm2` | 3 floats | 0 0 0 | Normal vector of cut plane 2 |
| `plane_distance2` | float | | Distance for plane 2 |
| `plane_point2` | 3 floats | | Point on plane 2 |
| `cutdir2` | float | 0 | Cut direction for plane 2 |

#### Output

| Keyword | Arguments | Default | Description |
|---------|-----------|---------|-------------|
| `outfile` | string | (required) | Output filename |
| `output_format` | string | `imd` | `imd` or `lammps` |
| `compress` | string | `no` | Gzip compress LAMMPS output (`yes`/`true`/`1` or `no`/`false`/`0`) |
| `startnr` | int | 1 | Starting atom ID number |

#### Miscellaneous

| Keyword | Arguments | Default | Description |
|---------|-----------|---------|-------------|
| `epsilon` | float | 1.0e-6 | Tolerance added to boundary checks. Prevents atoms sitting exactly on box faces |
| `shiftvec` | 3 floats | 0 0 0 | Translate all atoms by this vector (pre-rotation Cartesian coordinates) |

#### Wyckoff Overrides (Legacy)

| Keyword | Arguments | Description |
|---------|-----------|-------------|
| `pos_8d_init` | 3 floats | Override 8d Wyckoff position for D0_11 |
| `pos_4c_init` | 3 floats | Override 4c Wyckoff position |
| `pos_4c_init1` | 3 floats | Override second 4c Wyckoff position |

#### Deprecated Keywords

| Keyword | Notes |
|---------|-------|
| `lower` | Ignored. Loop bounds are computed automatically |
| `upper` | Ignored. Loop bounds are computed automatically |

---

## Crystal Structures

### Built-in Structures

LEGO ships with POSCAR files for common crystal structures in the `structures/`
directory:

**sc** (simple cubic): 1 atom per cell at (0, 0, 0). Unit lattice vectors.

**fcc** (face-centered cubic): 4 atoms per conventional cell at (0,0,0),
(0.5,0.5,0), (0.5,0,0.5), (0,0.5,0.5). Unit lattice vectors.

**bcc** (body-centered cubic): 2 atoms per conventional cell at (0,0,0)
and (0.5,0.5,0.5). Unit lattice vectors.

**hex / hcp** (hexagonal): 2 atoms per cell with hexagonal lattice vectors
a1 = (1, 0, 0), a2 = (-0.5, sqrt(3)/2, 0), a3 = (0, 0, c/a). The c/a
ratio is set in the POSCAR and can be overridden with `lattice_const_c`.

**L12** (e.g. Ni3Al): 4 atoms per cell. Type 0 (Al) at corners, type 1 (Ni)
at face centers. Use `mass` and `mass1` to set atomic masses.

**B2** (e.g. NiAl): 2 atoms per cell. Type 0 at (0,0,0), type 1 at
(0.5,0.5,0.5).

**D011** (e.g. Fe3C): Orthorhombic structure with multiple Wyckoff positions.
The Wyckoff positions can be overridden with `pos_8d_init`, `pos_4c_init`,
and `pos_4c_init1`.

### Custom Structures

Any structure can be used by providing a POSCAR file:

```
structure       poscar
structure_file  /path/to/my_structure.poscar
lattice_const   5.43
```

POSCAR files can be generated from CIF or other formats using the Python
wrapper's conversion functions (see [Python API](#python-api)).

---

## POSCAR Format

LEGO reads VASP POSCAR files (both VASP 4 and VASP 5+ formats). The file
format is:

```
Comment line
scale_factor
  a1x a1y a1z       # lattice vector a1
  a2x a2y a2z       # lattice vector a2
  a3x a3y a3z       # lattice vector a3
Element1 Element2    # element symbols (VASP 5+ only)
N1 N2               # atom counts per type
Direct               # or Cartesian
  x1 y1 z1          # fractional (or Cartesian) coordinates
  x2 y2 z2
  ...
```

Key points:

- The scale factor multiplies all lattice vectors. Set to 1.0 if vectors are
  already in Angstrom.
- The element symbols line is optional (VASP 4 compatibility). If absent,
  types are named T1, T2, etc., and masses default to 1.0 AMU.
- Coordinates can be `Direct` (fractional) or `Cartesian`.
- The `Selective dynamics` line, if present, is skipped.
- LEGO recognizes common element symbols and assigns standard atomic masses
  automatically (H through Pb). Unknown elements get mass 1.0 AMU.
- The `lattice_const` parameter overrides the lattice scale by rescaling
  all vectors proportionally to match the specified value for |a1|.

---

## Rotation

LEGO supports three rotation methods. They are applied in priority order:
explicit matrix > crystallographic directions > axis/angle. The final rotation
is R_total = R2 * R1, where R2 comes from axis2/angle2 (default: identity).

### Method 1: Crystallographic Directions (new_x / new_z)

The most common method. Specify two of the three crystallographic directions
that should align with the Cartesian axes. The third direction is computed
via cross product.

To create a (111)-oriented FCC crystal with x along [1-10]:

```
new_x    1 -1  0
new_z    1  1  1
```

LEGO normalizes the two given vectors and computes the missing one:

- Given `new_x` and `new_z`: `new_y = cross(new_z, new_x)`
- Given `new_x` and `new_y`: `new_z = cross(new_x, new_y)`
- Given `new_y` and `new_z`: `new_x = cross(new_y, new_z)`

The rotation matrix R1 has the normalized direction vectors as its rows:

```
R1 = | nx_x  nx_y  nx_z |
     | ny_x  ny_y  ny_z |
     | nz_x  nz_y  nz_z |
```

The two specified directions must not be parallel.

### Method 2: Axis/Angle

For arbitrary rotations not aligned with crystallographic directions:

```
axis1    0 0 1
angle1   45.0
```

This rotates by 45 degrees around the z-axis. A second rotation can be
chained:

```
axis2    1 0 0
angle2   30.0
```

The combined rotation is R_total = R2 * R1 (R1 applied first).

### Method 3: Explicit Matrix

For full control, specify all nine entries of the rotation matrix:

```
rotmat1_x   0.707107  -0.707107  0.0
rotmat1_y   0.707107   0.707107  0.0
rotmat1_z   0.0        0.0       1.0
```

The determinant must be +/-1 (to machine precision). Setting any of the three
`rotmat1_*` keywords activates this mode and overrides directions and
axis/angle.

### How Rotation Affects the Lattice

The rotated lattice matrix is M = R_total * L^T, where L contains the
(possibly rescaled) lattice vectors as rows. The position of atom (i,j,k,b)
is:

```
pos = M * (i + bx, j + by, k + bz)^T + R_total * shiftvec
```

where (bx, by, bz) are the fractional basis coordinates of atom b.

---

## Box Sizing

### Periodic Mode (default)

In periodic mode, LEGO adjusts each box dimension to the largest integer
multiple of the minimum repeat distance that does not exceed the target size.
This ensures the box is commensurate with the (rotated) lattice so that
periodic boundary conditions are artifact-free.

The minimum repeat distance along each Cartesian axis is the smallest
translation L such that M^{-1} * (L, 0, 0) (for x-axis, etc.) yields
integer lattice coordinates. LEGO finds this by searching integer multiples.

The actual box dimension is:

```
box_actual = n * repeat_dist - epsilon
```

where n = floor(box_target / repeat_dist) and epsilon prevents atoms from
sitting exactly on the boundary.

### Exact Mode

In exact mode, box dimensions are used as specified without adjustment:

```
box_mode    exact
```

This is useful for non-periodic simulations or when the box size is already
known to be commensurate.

### Loop Bounds

LEGO automatically computes the range of lattice indices (i, j, k) needed
to fill the box. It transforms all eight box corners through M^{-1} to find
the fractional coordinate extremes, then adds a safety margin. The old
`lower`/`upper` parameters are deprecated and ignored.

---

## Cut Planes

Two independent cut planes can remove atoms from the generated crystal.
Each plane is defined by:

1. A **normal vector** (`cut_norm1` / `cut_norm2`), which is automatically
   normalized.
2. A **distance** from the origin (`plane_distance1`) or a **point** on the
   plane (`plane_point1`). The plane distance is d = dot(point, normal).
3. A **cut direction** (`cutdir1`):
   - `+1`: remove atoms where dot(pos, norm) >= d (far side of origin)
   - `-1`: remove atoms where dot(pos, norm) <= d (near side of origin)
   - `0`: disable this cut plane

### Example: Removing a Half-Space

Remove all atoms with x > 25:

```
cut_norm1        1.0 0.0 0.0
plane_distance1  25.0
cutdir1          1
```

### Example: Creating a Slab

Remove atoms below z=10 and above z=40:

```
cut_norm1        0.0 0.0 1.0
plane_distance1  10.0
cutdir1         -1

cut_norm2        0.0 0.0 1.0
plane_distance2  40.0
cutdir2          1
```

### Example: Angled Cut Using a Point

Define the plane by a point rather than a distance:

```
cut_norm1        1.0 1.0 0.0
plane_point1     25.0 25.0 0.0
cutdir1          1
```

---

## Output Formats

### IMD Format

The IMD format uses an ASCII header followed by atom data:

```
#F A 1 1 1 3 0 0
#C number type mass x y z
#X box_x 0 0
#Y 0 box_y 0
#Z 0 0 box_z
##
#E
id type mass x y z
```

- The box is orthogonal (off-diagonal terms are zero).
- Types are 0-based integers.
- Masses are in IMD internal units: AMU * 0.000103650 (conversion factor
  from AMU to eV ps^2 / Angstrom^2).
- Positions are in Angstrom.

### LAMMPS Data Format

The LAMMPS output uses the `atomic` atom style:

```
LEGO generated structure: comment

N atoms
M atom types

0.0 box_x xlo xhi
0.0 box_y ylo yhi
0.0 box_z zlo zhi

Masses

1 mass_in_AMU
2 mass_in_AMU

Atoms # atomic

id type x y z
```

- Types are 1-based integers (shifted from internal 0-based).
- Masses are in AMU.
- Positions are in Angstrom.
- Gzip compression is available: set `compress yes` in the parameter file.
  The output filename gets `.gz` appended automatically if not already present.

---

## Python API

The Python wrapper in `python/lego.py` provides a `Lego` class for scripting.

### Installation

No installation required. Add the `python/` directory to your Python path or
import directly:

```python
import sys
sys.path.insert(0, "/path/to/lego/python")
from lego import Lego
```

### Lego Class

```python
lego = Lego(binary=None)
```

Initializes the wrapper. If `binary` is not given, it searches:
1. `PATH` for a `lego` executable
2. The parent directory of `lego.py`
3. A `build/` subdirectory of the parent

### Running LEGO

**With a parameter file:**

```python
result = lego.run("Cu_fcc.param")
```

**With keyword arguments:**

```python
result = lego.run(
    structure="bcc",
    lattice_const=3.1652,
    mass=183.84,
    box_x=50, box_y=50, box_z=50,
    new_x="1 -1 0", new_z="1 1 0",
    output_format="imd",
    outfile="W_bcc_110.imd"
)
```

A temporary parameter file is generated, passed to the binary, and deleted
after the run.

### LegoResult

The `run()` method returns a `LegoResult` object:

| Attribute | Type | Description |
|-----------|------|-------------|
| `total_atoms` | int | Number of atoms generated |
| `box` | tuple(float, float, float) | Actual box dimensions |
| `repeat_dist` | tuple(float, float, float) | Minimum repeat distances |
| `success` | bool | Whether the run exited with code 0 |
| `outfile` | str | Output filename |
| `stdout` | str | Raw standard output |
| `stderr` | str | Raw standard error |
| `returncode` | int | Process exit code |

### Structure Conversion

**CIF to POSCAR** (requires pymatgen):

```python
Lego.cif_to_poscar("input.cif", "output.poscar")
```

**Any format to POSCAR** (requires ASE):

```python
Lego.ase_to_poscar("input.xyz", "output.poscar")
Lego.ase_to_poscar("input.cfg", "output.poscar", format="cfg")
```

### Command-Line Usage

The wrapper can also be run as a script:

```bash
python python/lego.py Cu_fcc.param
python python/lego.py --binary ./lego Cu_fcc.param
python python/lego.py --cif-to-poscar structure.cif structure.poscar
```

---

## Migration from Old LEGO

LEGO v2.0 is a rewrite of the original `mk_config` code by Erik Bitzek
(2000-2011). The parameter file format is backward-compatible with the
following changes:

### Removed/Deprecated Parameters

| Old Parameter | Status | Notes |
|---------------|--------|-------|
| `lower` | Deprecated | Loop bounds are now computed automatically from box geometry. A warning is printed if encountered |
| `upper` | Deprecated | Same as `lower` |

### New Parameters

| Parameter | Description |
|-----------|-------------|
| `structure_file` | Path to custom POSCAR file (use with `structure poscar`) |
| `box_mode` | Choose `periodic` (default) or `exact` box sizing |
| `output_format` | Choose `imd` (default) or `lammps` |
| `compress` | Gzip LAMMPS output |
| `lattice_const_b` | Independent b-axis lattice constant |
| `lattice_const_c` | Independent c-axis lattice constant |
| `mass1` | Mass override for type 1 (in addition to `mass` for type 0) |
| `rotmat1_x/y/z` | Explicit rotation matrix specification |
| `plane_point1/2` | Alternative to `plane_distance` using a point on the plane |

### Behavioral Changes

- **Automatic loop bounds**: The `lower`/`upper` parameters are no longer
  needed. LEGO computes the minimum (i,j,k) range by transforming box
  corners through M^{-1} and adding a safety margin.

- **POSCAR-based structures**: All built-in structures are now defined as
  POSCAR files in the `structures/` directory rather than being hardcoded.
  This makes it easy to add new structures by dropping in a POSCAR file.

- **Mass handling**: Masses are read from POSCAR element symbols (with a
  built-in lookup table for common elements) rather than requiring manual
  specification. The `mass`/`mass1` parameters override the automatic values.

- **IMD mass units**: The old code expected masses in IMD internal units.
  The new code expects AMU and converts automatically (factor: 0.000103650).
  If migrating old parameter files, divide old mass values by 0.000103650
  to get AMU, or simply use standard AMU values.

- **LAMMPS output**: LAMMPS data file output is new in v2.0. Set
  `output_format lammps` to use it.

- **OpenMP parallelization**: The atom counting pass is parallelized with
  OpenMP for faster execution on large systems. Build with OpenMP support
  for best performance.

### Migrating a Parameter File

A minimal old-style parameter file:

```
structure    bcc
lower        -200
upper         200
epsilon      1.0e-6
outfile      W.bcc
startnr      1
mass         0.01905605250
lattice_const 3.1652
box_x        50
box_y        50
box_z        50
cut_norm1    1.0 0.0 0.0
plane_distance1 0
cutdir1      0
cut_norm2    1.0 0.0 0.0
plane_distance2 0
cutdir2      0
```

Equivalent v2.0 parameter file:

```
structure       bcc
lattice_const   3.1652
mass            183.84

box_x           50
box_y           50
box_z           50

box_mode        periodic
output_format   imd
outfile         W.imd
```

Key differences:
- `lower`/`upper` removed (automatic)
- `mass` changed from IMD units (0.019) to AMU (183.84)
- Cut planes with `cutdir 0` can be omitted entirely
- `epsilon` and `startnr` use sensible defaults and can be omitted
