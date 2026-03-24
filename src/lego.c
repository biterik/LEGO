/* LEGO - Lattice Elemental Geometry Operations
 *
 * Crystal structure generator for molecular dynamics simulations.
 * Creates arbitrarily rotated and cut crystal structures in
 * IMD or LAMMPS data format.
 *
 * Originally by Erik Bitzek (2000-2011)
 * Rewritten 2026
 *
 * Licensed under GPLv3
 */

#include "lego.h"
#include <stdarg.h>

/* Defined in poscar.c */
extern void poscar_set_argv0(const char *argv0);

/* ========== Atom Position and Filtering ========== */

vec3 compute_atom_pos(const Config *cfg, int i, int j, int k, int b)
{
    /* Fractional coordinates: lattice point + basis */
    double fi = (double)i + cfg->crystal.basis[b].x;
    double fj = (double)j + cfg->crystal.basis[b].y;
    double fk = (double)k + cfg->crystal.basis[b].z;

    /* Transform to Cartesian via M: pos = M * (fi, fj, fk)^T */
    vec3 pos;
    pos.x = cfg->M[0][0] * fi + cfg->M[0][1] * fj + cfg->M[0][2] * fk;
    pos.y = cfg->M[1][0] * fi + cfg->M[1][1] * fj + cfg->M[1][2] * fk;
    pos.z = cfg->M[2][0] * fi + cfg->M[2][1] * fj + cfg->M[2][2] * fk;

    /* Apply rotated shift vector */
    pos.x += cfg->shiftvec_rotated.x;
    pos.y += cfg->shiftvec_rotated.y;
    pos.z += cfg->shiftvec_rotated.z;

    return pos;
}

int atom_in_box(const Config *cfg, vec3 pos)
{
    double eps = cfg->epsilon;
    return (pos.x + eps >= 0.0 && pos.x + eps < cfg->box_actual[0] &&
            pos.y + eps >= 0.0 && pos.y + eps < cfg->box_actual[1] &&
            pos.z + eps >= 0.0 && pos.z + eps < cfg->box_actual[2]);
}

int atom_cut_away(const Config *cfg, vec3 pos)
{
    /* Check cut plane 1 */
    if (cfg->cutdir1 != 0.0) {
        double pn1 = vec3_dot(pos, cfg->cut_norm1);
        if (cfg->cutdir1 < 0.0 && pn1 <= cfg->plane_dist1) return 1;
        if (cfg->cutdir1 > 0.0 && pn1 >= cfg->plane_dist1) return 1;
    }

    /* Check cut plane 2 */
    if (cfg->cutdir2 != 0.0) {
        double pn2 = vec3_dot(pos, cfg->cut_norm2);
        if (cfg->cutdir2 < 0.0 && pn2 <= cfg->plane_dist2) return 1;
        if (cfg->cutdir2 > 0.0 && pn2 >= cfg->plane_dist2) return 1;
    }

    return 0;
}

/* ========== Atom Counting (Pass 1) ========== */

long count_atoms(const Config *cfg)
{
    long count = 0;
    int imin = cfg->ijk_min[0], imax = cfg->ijk_max[0];
    int jmin = cfg->ijk_min[1], jmax = cfg->ijk_max[1];
    int kmin = cfg->ijk_min[2], kmax = cfg->ijk_max[2];
    int natoms = cfg->crystal.natoms;

    #ifdef _OPENMP
    #pragma omp parallel for collapse(3) reduction(+:count) schedule(dynamic, 64)
    #endif
    for (int i = imin; i <= imax; i++)
        for (int j = jmin; j <= jmax; j++)
            for (int k = kmin; k <= kmax; k++) {
                for (int b = 0; b < natoms; b++) {
                    vec3 pos = compute_atom_pos(cfg, i, j, k, b);
                    if (atom_in_box(cfg, pos) && !atom_cut_away(cfg, pos))
                        count++;
                }
            }

    return count;
}

/* ========== Atom Generation (Pass 2) ========== */

void generate_atoms(const Config *cfg, OutputWriter *writer)
{
    int imin = cfg->ijk_min[0], imax = cfg->ijk_max[0];
    int jmin = cfg->ijk_min[1], jmax = cfg->ijk_max[1];
    int kmin = cfg->ijk_min[2], kmax = cfg->ijk_max[2];
    int natoms = cfg->crystal.natoms;
    int nr = cfg->startnr;

    for (int i = imin; i <= imax; i++)
        for (int j = jmin; j <= jmax; j++)
            for (int k = kmin; k <= kmax; k++) {
                for (int b = 0; b < natoms; b++) {
                    vec3 pos = compute_atom_pos(cfg, i, j, k, b);

                    if (!atom_in_box(cfg, pos))
                        continue;
                    if (atom_cut_away(cfg, pos))
                        continue;

                    int type = cfg->crystal.basis_type[b];
                    if (cfg->atomtype_override >= 0)
                        type = cfg->atomtype_override;

                    double mass = cfg->crystal.masses[
                        cfg->crystal.basis_type[b]];

                    writer->write_atom(writer->ctx, nr, type, mass, pos);
                    nr++;
                }
            }
}

/* ========== Usage ========== */

static void usage(const char *progname)
{
    printf("LEGO - Lattice Elemental Geometry Operations\n");
    printf("Crystal structure generator for MD simulations\n\n");
    printf("Usage: %s <parameter_file>\n\n", progname);
    printf("Built-in structures: sc, fcc, bcc, hex, hcp, L12, B2, D011\n");
    printf("Or provide custom structures via POSCAR files.\n\n");
    printf("See README.md for parameter file format.\n");
}

/* ========== Main ========== */

int main(int argc, char *argv[])
{
    Config cfg;

    if (argc < 2) {
        usage(argv[0]);
        exit(1);
    }

    /* Initialize */
    config_init(&cfg);
    poscar_set_argv0(argv[0]);

    printf("========================================\n");
    printf("  LEGO - Lattice Elemental Geometry Operations\n");
    printf("========================================\n\n");

    #ifdef _OPENMP
    printf("OpenMP enabled: %d threads available\n\n",
           omp_get_max_threads());
    #endif

    /* Read parameter file */
    printf("Reading parameters from: %s\n", argv[1]);
    params_read(argv[1], &cfg);

    /* Load crystal structure */
    crystal_load(&cfg);

    /* Set up lattice: rotation, M, M^{-1}, bounds, periodic sizing */
    lattice_setup(&cfg);

    /* Pass 1: Count atoms */
    printf("Counting atoms...\n");
    cfg.total_atoms = count_atoms(&cfg);
    printf("Total atoms: %ld\n\n", cfg.total_atoms);

    if (cfg.total_atoms == 0) {
        fprintf(stderr, "WARNING: No atoms generated. Check box dimensions "
                "and crystal parameters.\n");
        return 1;
    }

    /* Open output writer */
    OutputWriter writer;
    if (cfg.output_format == FMT_IMD) {
        printf("Output format: IMD\n");
        writer = imd_open(cfg.outfile);
    } else {
        printf("Output format: LAMMPS data\n");
        writer = lammps_open(cfg.outfile, cfg.compress);
    }
    printf("Output file: %s\n\n", cfg.outfile);

    /* Write header */
    writer.write_header(writer.ctx, &cfg);

    /* Pass 2: Generate atoms */
    printf("Generating atoms...\n");
    generate_atoms(&cfg, &writer);

    /* Finalize */
    writer.write_footer(writer.ctx, &cfg);
    writer.close(writer.ctx);

    /* Summary */
    printf("\n========================================\n");
    printf("  Summary\n");
    printf("========================================\n");
    printf("Structure: %s\n", cfg.crystal.comment);
    printf("Atoms: %ld\n", cfg.total_atoms);
    printf("Box: %.8f x %.8f x %.8f\n",
           cfg.box_actual[0], cfg.box_actual[1], cfg.box_actual[2]);
    printf("Min repeat: %.8f  %.8f  %.8f\n",
           cfg.repeat_dist[0], cfg.repeat_dist[1], cfg.repeat_dist[2]);
    printf("Output: %s (%s)\n", cfg.outfile,
           cfg.output_format == FMT_IMD ? "IMD" : "LAMMPS");
    printf("========================================\n");

    return 0;
}
