# LEGO - Lattice Elemental Geometry Operations

Crystal structure generator for molecular dynamics simulations. Creates
arbitrarily rotated and cut crystal structures in IMD and LAMMPS data formats.

Originally by Erik Bitzek (2000-2011). Rewritten 2026.

## Building

Requires a C compiler, zlib, and optionally OpenMP.

### Linux

```bash
make
```

GCC with OpenMP is typically available out of the box.

### macOS

The default Apple Clang does not support OpenMP. LEGO builds without it
(single-threaded counting pass), which is fine for smaller systems:

```bash
make
```

For OpenMP support (recommended for large systems), install GCC via Homebrew:

```bash
brew install gcc
make CC=gcc-15 OMP_FLAG="-fopenmp"
```

Replace `gcc-15` with the version Homebrew installed (check with
`ls /opt/homebrew/bin/gcc-*`).

### Custom Compiler

```bash
make CC=icc CFLAGS="-O3 -xHost"
```

The build produces a `lego` binary in the project root. The `structures/`
directory must be accessible at runtime (see [Structure Search Path](#structure-search-path)).

## Quick Start

1. Create a parameter file `Cu_fcc.param`:

```
structure       fcc
lattice_const   3.615
mass            63.546

new_x           1 -1 0
new_z           1  1 1

box_x           50
box_y           50
box_z           50

box_mode        periodic
output_format   imd
outfile         Cu_fcc_111.imd
```

2. Run LEGO:

```bash
./lego Cu_fcc.param
```

This generates a (111)-oriented FCC copper crystal with periodic box dimensions
of approximately 50 Angstrom in each direction.

## Parameter Reference

Parameters are specified one per line in the form `keyword value`. Lines
starting with `#` are comments.

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `structure` | Crystal structure name (see [Built-in Structures](#built-in-structures)) or `poscar` for custom |
| `outfile` | Output filename |
| `box_x` | Target box size in x (Angstrom) |
| `box_y` | Target box size in y (Angstrom) |
| `box_z` | Target box size in z (Angstrom) |

### Crystal Definition

| Parameter | Default | Description |
|-----------|---------|-------------|
| `structure` | (required) | Built-in name or `poscar` |
| `structure_file` | | Path to POSCAR file (required when `structure poscar`) |
| `lattice_const` | from POSCAR | Override lattice constant a (Angstrom) |
| `lattice_const_b` | from POSCAR | Override lattice constant b |
| `lattice_const_c` | from POSCAR | Override lattice constant c |
| `mass` | from POSCAR | Override mass for type 0 (AMU) |
| `mass1` | from POSCAR | Override mass for type 1 (AMU) |

### Rotation

Three methods are available (highest priority first):

| Parameter | Description |
|-----------|-------------|
| `rotmat1_x` | Explicit rotation matrix row 1 (3 values) |
| `rotmat1_y` | Explicit rotation matrix row 2 (3 values) |
| `rotmat1_z` | Explicit rotation matrix row 3 (3 values) |
| `new_x` | Crystallographic direction for x-axis (3 Miller indices) |
| `new_y` | Crystallographic direction for y-axis (3 Miller indices) |
| `new_z` | Crystallographic direction for z-axis (3 Miller indices) |
| `axis1` | Rotation axis for first rotation (3 values) |
| `angle1` | Rotation angle for first rotation (degrees) |
| `axis2` | Rotation axis for second rotation (3 values) |
| `angle2` | Rotation angle for second rotation (degrees) |

When using `new_x`/`new_z`, specify any two of the three directions. The
third is computed automatically via cross product.

### Box and Output

| Parameter | Default | Description |
|-----------|---------|-------------|
| `box_mode` | `periodic` | `periodic` (round to repeat distance) or `exact` |
| `output_format` | `imd` | `imd` or `lammps` |
| `compress` | `no` | Gzip LAMMPS output (`yes`/`no`) |
| `epsilon` | `1.0e-6` | Tolerance for boundary checks |
| `startnr` | `1` | Starting atom ID number |
| `atomtype` | (none) | Override all atom types to this value |

### Cut Planes

Two independent cut planes are available. Each requires a normal vector, a
distance or point, and a direction:

| Parameter | Description |
|-----------|-------------|
| `cut_norm1` | Normal vector of cut plane 1 (3 values) |
| `plane_distance1` | Signed distance from origin to plane 1 |
| `plane_point1` | Point on plane 1 (alternative to distance) |
| `cutdir1` | Cut direction: `+1` remove far side, `-1` remove near side, `0` disable |
| `cut_norm2` | Normal vector of cut plane 2 |
| `plane_distance2` | Signed distance from origin to plane 2 |
| `plane_point2` | Point on plane 2 |
| `cutdir2` | Cut direction for plane 2 |

### Miscellaneous

| Parameter | Description |
|-----------|-------------|
| `shiftvec` | Shift all atoms by this vector (3 values, pre-rotation Cartesian) |

## Built-in Structures

| Name | POSCAR File | Description |
|------|-------------|-------------|
| `sc` | `sc.poscar` | Simple cubic |
| `fcc` | `fcc.poscar` | Face-centered cubic |
| `bcc` | `bcc.poscar` | Body-centered cubic |
| `hex` / `hcp` | `hex.poscar` | Hexagonal close-packed |
| `L12` | `L12_Ni3Al.poscar` | L1_2 ordered (e.g. Ni3Al) |
| `B2` | `B2_NiAl.poscar` | B2 ordered (e.g. NiAl) |
| `D011` | `D011_Fe3C.poscar` | D0_11 ordered (e.g. Fe3C) |

Custom structures can be used by setting `structure poscar` and providing a
`structure_file` path to any valid VASP POSCAR file.

### Structure Search Path

LEGO searches for built-in POSCAR files in this order:

1. `$LEGO_STRUCTURES` environment variable
2. `structures/` relative to the current working directory
3. `structures/` relative to the executable location

## Output Formats

### IMD

The default format. Produces plain text files with an orthogonal box header:

```
#F A 1 1 1 3 0 0
#C number type mass x y z
#X box_x 0 0
#Y 0 box_y 0
#Z 0 0 box_z
##
#E
id type mass x y z
...
```

Masses are converted to IMD internal units (eV ps^2 / Angstrom^2).

### LAMMPS

Produces LAMMPS data files in `atomic` style:

```
LEGO generated structure: ...

N atoms
M atom types

0.0 box_x xlo xhi
0.0 box_y ylo yhi
0.0 box_z zlo zhi

Masses

1 mass1
...

Atoms # atomic

id type x y z
...
```

LAMMPS types are 1-based. Masses are in AMU. Gzip compression is available
via the `compress` parameter.

## Python Wrapper

A Python wrapper is provided in `python/lego.py`:

```python
from lego import Lego

# Initialize (finds the binary automatically)
lego = Lego()

# Run with a parameter file
result = lego.run("params.param")
print(result.total_atoms, result.box)

# Run with keyword arguments (generates a temp parameter file)
result = lego.run(
    structure="fcc",
    lattice_const=3.615,
    mass=63.546,
    box_x=100, box_y=100, box_z=100,
    new_x="1 1 0", new_z="0 0 1",
    output_format="lammps",
    outfile="Cu_fcc.lmp"
)

# Convert CIF to POSCAR (requires pymatgen)
Lego.cif_to_poscar("structure.cif", "structure.poscar")

# Convert any format to POSCAR (requires ASE)
Lego.ase_to_poscar("structure.xyz", "structure.poscar")
```

The `LegoResult` object provides:
- `total_atoms` -- number of atoms generated
- `box` -- tuple of (box_x, box_y, box_z)
- `repeat_dist` -- tuple of minimum repeat distances
- `success` -- whether the run succeeded
- `stdout` / `stderr` -- raw output

## Examples

See the `examples/` directory for ready-to-run parameter files:

- `W_bcc_100.param` -- BCC tungsten, default (100) orientation
- `W_bcc_110.param` -- BCC tungsten, (110) orientation
- `Cu_fcc_111.param` -- FCC copper, (111) orientation
- `Cu_fcc_lammps.param` -- FCC copper, LAMMPS output
- `Ni3Al_L12.param` -- L1_2 Ni3Al ordered alloy
- `NiAl_B2.param` -- B2 NiAl ordered alloy
- `custom_poscar.param` -- Using a custom POSCAR file

## License

GPLv3. See [LICENSE](LICENSE).

## Citation

If you use LEGO in your research, please cite it. See [CITATION.cff](CITATION.cff)
for metadata.
