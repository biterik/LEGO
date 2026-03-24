/* LEGO - output_imd.c
 * IMD format output writer.
 *
 * IMD format:
 *   #F A 1 1 1 3 0 0
 *   #C number type mass x y z
 *   #X box_x 0 0
 *   #Y 0 box_y 0
 *   #Z 0 0 box_z
 *   ##
 *   #E
 *   id type mass x y z
 *   ...
 */

#include "lego.h"

typedef struct {
    FILE *fp;
} ImdCtx;

static int imd_write_header(void *ctx, const Config *cfg)
{
    ImdCtx *c = (ImdCtx *)ctx;

    fprintf(c->fp, "#F A 1 1 1 3 0 0\n");
    fprintf(c->fp, "#C number type mass x y z\n");
    fprintf(c->fp, "#X %.16f 0 0\n", cfg->box_actual[0]);
    fprintf(c->fp, "#Y 0 %.16f 0\n", cfg->box_actual[1]);
    fprintf(c->fp, "#Z 0 0 %.16f\n", cfg->box_actual[2]);
    fprintf(c->fp, "##\n");
    fprintf(c->fp, "#E\n");

    return 0;
}

static int imd_write_atom(void *ctx, int id, int type, double mass, vec3 pos)
{
    ImdCtx *c = (ImdCtx *)ctx;

    /* IMD uses 0-based types (already 0-based internally) */
    /* Mass is converted to IMD internal units */
    double imd_mass = mass * AMU_TO_IMD;

    fprintf(c->fp, "%d %d %.16f %.16f %.16f %.16f\n",
            id, type, imd_mass, pos.x, pos.y, pos.z);

    return 0;
}

static int imd_write_footer(void *ctx, const Config *cfg)
{
    (void)ctx;
    (void)cfg;
    /* IMD has no footer */
    return 0;
}

static void imd_close(void *ctx)
{
    ImdCtx *c = (ImdCtx *)ctx;
    if (c->fp)
        fclose(c->fp);
    free(c);
}

OutputWriter imd_open(const char *filename)
{
    OutputWriter w;
    ImdCtx *c = (ImdCtx *)calloc(1, sizeof(ImdCtx));

    if (!c) {
        fprintf(stderr, "Out of memory\n");
        exit(2);
    }

    c->fp = fopen(filename, "w");
    if (!c->fp) {
        fprintf(stderr, "Can't open output file: %s\n", filename);
        exit(2);
    }

    w.ctx = c;
    w.write_header = imd_write_header;
    w.write_atom = imd_write_atom;
    w.write_footer = imd_write_footer;
    w.close = imd_close;

    return w;
}
