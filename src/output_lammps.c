/* LEGO - output_lammps.c
 * LAMMPS data format output writer.
 *
 * LAMMPS data format (atomic style):
 *   LEGO generated structure
 *
 *   N atoms
 *   M atom types
 *
 *   0.0 box_x xlo xhi
 *   0.0 box_y ylo yhi
 *   0.0 box_z zlo zhi
 *
 *   Masses
 *
 *   1 mass1
 *   2 mass2
 *
 *   Atoms # atomic
 *
 *   id type x y z
 *   ...
 */

#include "lego.h"
#include <zlib.h>

typedef struct {
    FILE *fp;       /* for plain text */
    gzFile gz;      /* for gzip */
    int compressed;
} LammpsCtx;

/* Wrapper for writing to either FILE* or gzFile */
static void lmp_printf(LammpsCtx *c, const char *fmt, ...)
{
    va_list args;
    va_start(args, fmt);

    if (c->compressed) {
        char buf[1024];
        vsnprintf(buf, sizeof(buf), fmt, args);
        gzputs(c->gz, buf);
    } else {
        vfprintf(c->fp, fmt, args);
    }

    va_end(args);
}

static int lammps_write_header(void *ctx, const Config *cfg)
{
    LammpsCtx *c = (LammpsCtx *)ctx;

    lmp_printf(c, "LEGO generated structure: %s\n\n", cfg->crystal.comment);
    lmp_printf(c, "%ld atoms\n", cfg->total_atoms);
    lmp_printf(c, "%d atom types\n\n", cfg->crystal.ntypes);

    lmp_printf(c, "0.0 %.16f xlo xhi\n", cfg->box_actual[0]);
    lmp_printf(c, "0.0 %.16f ylo yhi\n", cfg->box_actual[1]);
    lmp_printf(c, "0.0 %.16f zlo zhi\n\n", cfg->box_actual[2]);

    lmp_printf(c, "Masses\n\n");
    for (int t = 0; t < cfg->crystal.ntypes; t++) {
        /* LAMMPS types are 1-based, masses in AMU */
        lmp_printf(c, "%d %.6f\n", t + 1, cfg->crystal.masses[t]);
    }

    lmp_printf(c, "\nAtoms # atomic\n\n");

    return 0;
}

static int lammps_write_atom(void *ctx, int id, int type, double mass,
                             vec3 pos)
{
    LammpsCtx *c = (LammpsCtx *)ctx;
    (void)mass; /* mass is in the Masses section, not per atom */

    /* LAMMPS types are 1-based */
    lmp_printf(c, "%d %d %.16f %.16f %.16f\n",
               id, type + 1, pos.x, pos.y, pos.z);

    return 0;
}

static int lammps_write_footer(void *ctx, const Config *cfg)
{
    (void)ctx;
    (void)cfg;
    return 0;
}

static void lammps_close(void *ctx)
{
    LammpsCtx *c = (LammpsCtx *)ctx;
    if (c->compressed) {
        if (c->gz) gzclose(c->gz);
    } else {
        if (c->fp) fclose(c->fp);
    }
    free(c);
}

OutputWriter lammps_open(const char *filename, int compress)
{
    OutputWriter w;
    LammpsCtx *c = (LammpsCtx *)calloc(1, sizeof(LammpsCtx));

    if (!c) {
        fprintf(stderr, "Out of memory\n");
        exit(2);
    }

    c->compressed = compress;

    if (compress) {
        /* Append .gz if not already present */
        char gz_name[1088];
        int len = strlen(filename);
        if (len >= 3 && strcmp(filename + len - 3, ".gz") == 0)
            strncpy(gz_name, filename, sizeof(gz_name) - 1);
        else
            snprintf(gz_name, sizeof(gz_name), "%s.gz", filename);

        c->gz = gzopen(gz_name, "wb");
        if (!c->gz) {
            fprintf(stderr, "Can't open gzip output file: %s\n", gz_name);
            exit(2);
        }
        printf("Writing compressed LAMMPS data file: %s\n", gz_name);
    } else {
        c->fp = fopen(filename, "w");
        if (!c->fp) {
            fprintf(stderr, "Can't open output file: %s\n", filename);
            exit(2);
        }
    }

    w.ctx = c;
    w.write_header = lammps_write_header;
    w.write_atom = lammps_write_atom;
    w.write_footer = lammps_write_footer;
    w.close = lammps_close;

    return w;
}
