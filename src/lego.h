/* LEGO - Lattice Elemental Geometry Operations
 * Crystal structure generator for molecular dynamics simulations
 *
 * Originally by Erik Bitzek (2000-2011)
 * Rewritten 2026
 *
 * Licensed under GPLv3
 */

#ifndef LEGO_H
#define LEGO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* ========== Constants ========== */

#define MAX_BASIS    64   /* max atoms in a conventional cell basis */
#define MAX_TYPES     8   /* max distinct atom types */
#define EPS_DEFAULT   1.0e-6
#define MAX_REPEAT_SEARCH 1000  /* max iterations for repeat distance search */
#define INT_TOL       1.0e-4    /* tolerance for "is this an integer?" checks */

/* IMD mass unit conversion: AMU -> eV*ps^2/Angstrom^2 */
#define AMU_TO_IMD    0.000103650

#ifndef STRUCTURES_DIR
#define STRUCTURES_DIR "structures"
#endif

/* ========== Vector / Matrix Types ========== */

typedef struct {
    double x, y, z;
} vec3;

typedef double mat3[3][3]; /* mat3[row][col] */

/* ========== Enums ========== */

typedef enum { MODE_PERIODIC, MODE_EXACT } BoxMode;
typedef enum { FMT_IMD, FMT_LAMMPS } OutputFormat;

/* ========== Crystal Structure ========== */

typedef struct {
    char comment[256];
    double scale;
    mat3 lattice;                    /* rows = lattice vectors a1, a2, a3 */
    int ntypes;
    char species[MAX_TYPES][8];      /* element symbols */
    int counts[MAX_TYPES];           /* atoms per type */
    double masses[MAX_TYPES];        /* atomic masses in AMU */
    int natoms;                      /* sum of counts[] */
    vec3 basis[MAX_BASIS];           /* fractional coordinates */
    int basis_type[MAX_BASIS];       /* type index (0-based) */
} Crystal;

/* ========== Output Writer Interface ========== */

typedef struct Config Config; /* forward declaration */

typedef struct {
    void *ctx;
    int (*write_header)(void *ctx, const Config *cfg);
    int (*write_atom)(void *ctx, int id, int type, double mass, vec3 pos);
    int (*write_footer)(void *ctx, const Config *cfg);
    void (*close)(void *ctx);
} OutputWriter;

/* ========== Configuration ========== */

struct Config {
    /* Crystal definition */
    char structure_name[256];
    char structure_file[1024];
    Crystal crystal;

    /* Box parameters */
    BoxMode box_mode;
    double box_target[3];            /* target box dimensions (user input) */
    double box_actual[3];            /* actual box after periodic adjustment */
    double repeat_dist[3];           /* min repeat distances per axis */

    /* Rotation: crystallographic direction mapping */
    vec3 new_x, new_y, new_z;
    int has_new_x, has_new_y, has_new_z;

    /* Rotation: axis/angle */
    vec3 axis1, axis2;
    double angle1, angle2;

    /* Rotation: explicit matrix */
    mat3 rotmat1;
    int rotmat1_given;

    /* Cut planes */
    vec3 cut_norm1, cut_norm2;
    double plane_dist1, plane_dist2;
    int plane_dist1_given, plane_dist2_given;
    vec3 plane_point1, plane_point2;
    int plane_point1_given, plane_point2_given;
    double cutdir1, cutdir2;

    /* Shift vector (in Cartesian, pre-rotation) */
    vec3 shiftvec;
    vec3 shiftvec_rotated;

    /* Output */
    char outfile[1024];
    OutputFormat output_format;
    int startnr;
    int compress;
    int atomtype_override;           /* -1 = no override */

    /* Lattice constant overrides */
    double lattice_const;            /* 0 = use POSCAR */
    double lattice_const_b;          /* 0 = use POSCAR */
    double lattice_const_c;          /* 0 = use POSCAR */

    /* Mass overrides: -1.0 = use POSCAR/default */
    double mass_override[MAX_TYPES];

    /* Wyckoff overrides (legacy D011 support) */
    vec3 pos_8d_init, pos_4c_init, pos_4c_init1;
    int has_wyckoff;

    /* Misc */
    double epsilon;

    /* Derived quantities (computed by lattice_setup) */
    mat3 R_total;                    /* combined rotation matrix */
    mat3 M;                          /* rotated lattice vectors */
    mat3 Minv;                       /* inverse of M */
    int ijk_min[3];                  /* loop bounds */
    int ijk_max[3];

    /* Runtime */
    long total_atoms;
};

/* ========== Function Prototypes ========== */

/* math3d.c */
double   vec3_len(vec3 v);
vec3     vec3_add(vec3 a, vec3 b);
vec3     vec3_sub(vec3 a, vec3 b);
vec3     vec3_scale(vec3 v, double s);
vec3     vec3_normalize(vec3 v);
double   vec3_dot(vec3 a, vec3 b);
vec3     vec3_cross(vec3 a, vec3 b);
void     mat3_identity(mat3 m);
void     mat3_zero(mat3 m);
void     mat3_copy(mat3 dst, const mat3 src);
double   mat3_det(const mat3 m);
int      mat3_inverse(const mat3 m, mat3 inv);
vec3     mat3_mul_vec(const mat3 m, vec3 v);
void     mat3_mul_mat(const mat3 a, const mat3 b, mat3 result);
void     mat3_from_axis_angle(mat3 m, vec3 axis, double angle);
void     mat3_transpose(const mat3 m, mat3 result);
void     mat3_print(const mat3 m, const char *label);

/* params.c */
void config_init(Config *cfg);
int  params_read(const char *filename, Config *cfg);

/* poscar.c */
int  poscar_read(const char *filename, Crystal *cryst);
int  crystal_load(Config *cfg);

/* lattice.c */
int  lattice_setup(Config *cfg);
void lattice_print_info(const Config *cfg);

/* output_imd.c */
OutputWriter imd_open(const char *filename);

/* output_lammps.c */
OutputWriter lammps_open(const char *filename, int compress);

/* lego.c */
vec3 compute_atom_pos(const Config *cfg, int i, int j, int k, int b);
int  atom_in_box(const Config *cfg, vec3 pos);
int  atom_cut_away(const Config *cfg, vec3 pos);
long count_atoms(const Config *cfg);
void generate_atoms(const Config *cfg, OutputWriter *writer);

#endif /* LEGO_H */
