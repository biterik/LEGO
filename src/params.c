/* LEGO - params.c
 * Parameter file parser.
 * Adapted from mk_config_param.c by Erik Bitzek.
 */

#include "lego.h"

/* ========== Helpers ========== */

/* Read a double from the current token stream. Returns 0 on success. */
static int read_double(double *val, const char *name, int line)
{
    char *str = strtok(NULL, " \t\n");
    if (!str) {
        fprintf(stderr, "Missing value for '%s' at line %d\n", name, line);
        return -1;
    }
    *val = atof(str);
    return 0;
}

/* Read an int from the current token stream. */
static int read_int(int *val, const char *name, int line)
{
    char *str = strtok(NULL, " \t\n");
    if (!str) {
        fprintf(stderr, "Missing value for '%s' at line %d\n", name, line);
        return -1;
    }
    *val = atoi(str);
    return 0;
}

/* Read a string from the current token stream. */
static int read_str(char *val, int maxlen, const char *name, int line)
{
    char *str = strtok(NULL, " \t\n");
    if (!str) {
        fprintf(stderr, "Missing value for '%s' at line %d\n", name, line);
        return -1;
    }
    strncpy(val, str, maxlen - 1);
    val[maxlen - 1] = '\0';
    return 0;
}

/* Read a vec3 (3 doubles) from the current token stream. */
static int read_vec3(vec3 *v, const char *name, int line)
{
    if (read_double(&v->x, name, line) ||
        read_double(&v->y, name, line) ||
        read_double(&v->z, name, line))
        return -1;
    return 0;
}

/* ========== Config Initialization ========== */

void config_init(Config *cfg)
{
    memset(cfg, 0, sizeof(Config));

    cfg->box_mode = MODE_PERIODIC;
    cfg->output_format = FMT_IMD;
    cfg->epsilon = EPS_DEFAULT;
    cfg->startnr = 1;
    cfg->compress = 0;
    cfg->atomtype_override = -1;

    cfg->lattice_const = 0.0;
    cfg->lattice_const_b = 0.0;
    cfg->lattice_const_c = 0.0;

    for (int i = 0; i < MAX_TYPES; i++)
        cfg->mass_override[i] = -1.0;

    /* Default rotation axes */
    cfg->axis1.x = 1.0;
    cfg->axis2.x = 1.0;

    /* Cut planes: default = no cutting */
    cfg->cutdir1 = 0.0;
    cfg->cutdir2 = 0.0;

    /* Sentinel for plane distances */
    cfg->plane_dist1 = -999999.9999;
    cfg->plane_dist2 = -999999.9999;
    cfg->plane_dist1_given = 0;
    cfg->plane_dist2_given = 0;
    cfg->plane_point1_given = 0;
    cfg->plane_point2_given = 0;

    cfg->has_new_x = 0;
    cfg->has_new_y = 0;
    cfg->has_new_z = 0;
    cfg->rotmat1_given = 0;
    cfg->has_wyckoff = 0;
}

/* ========== Parameter File Reader ========== */

int params_read(const char *filename, Config *cfg)
{
    FILE *infile;
    char buffer[1024];
    char *token;
    char *res;
    int curline = 0;

    infile = fopen(filename, "r");
    if (!infile) {
        fprintf(stderr, "Can't open parameter file: %s\n", filename);
        exit(2);
    }

    do {
        res = fgets(buffer, sizeof(buffer), infile);
        if (!res) break;
        curline++;

        token = strtok(buffer, " \t\n");
        if (!token) continue;          /* blank line */
        if (token[0] == '#') continue; /* comment */

        if (strcasecmp(token, "outfile") == 0) {
            read_str(cfg->outfile, sizeof(cfg->outfile), "outfile", curline);
        }
        else if (strcasecmp(token, "structure") == 0) {
            read_str(cfg->structure_name, sizeof(cfg->structure_name),
                     "structure", curline);
        }
        else if (strcasecmp(token, "structure_file") == 0) {
            read_str(cfg->structure_file, sizeof(cfg->structure_file),
                     "structure_file", curline);
        }
        else if (strcasecmp(token, "output_format") == 0) {
            char fmt[64];
            read_str(fmt, sizeof(fmt), "output_format", curline);
            if (strcasecmp(fmt, "imd") == 0)
                cfg->output_format = FMT_IMD;
            else if (strcasecmp(fmt, "lammps") == 0)
                cfg->output_format = FMT_LAMMPS;
            else {
                fprintf(stderr, "Unknown output_format '%s' at line %d. "
                        "Use 'imd' or 'lammps'.\n", fmt, curline);
                exit(2);
            }
        }
        else if (strcasecmp(token, "box_mode") == 0) {
            char mode[64];
            read_str(mode, sizeof(mode), "box_mode", curline);
            if (strcasecmp(mode, "periodic") == 0)
                cfg->box_mode = MODE_PERIODIC;
            else if (strcasecmp(mode, "exact") == 0)
                cfg->box_mode = MODE_EXACT;
            else {
                fprintf(stderr, "Unknown box_mode '%s' at line %d. "
                        "Use 'periodic' or 'exact'.\n", mode, curline);
                exit(2);
            }
        }
        else if (strcasecmp(token, "compress") == 0) {
            char val[64];
            read_str(val, sizeof(val), "compress", curline);
            cfg->compress = (strcasecmp(val, "yes") == 0 ||
                             strcasecmp(val, "1") == 0 ||
                             strcasecmp(val, "true") == 0);
        }
        else if (strcasecmp(token, "lower") == 0) {
            int dummy;
            read_int(&dummy, "lower", curline);
            fprintf(stderr, "WARNING: 'lower' is deprecated and ignored "
                    "(bounds computed automatically).\n");
        }
        else if (strcasecmp(token, "upper") == 0) {
            int dummy;
            read_int(&dummy, "upper", curline);
            fprintf(stderr, "WARNING: 'upper' is deprecated and ignored "
                    "(bounds computed automatically).\n");
        }
        else if (strcasecmp(token, "epsilon") == 0) {
            read_double(&cfg->epsilon, "epsilon", curline);
        }
        else if (strcasecmp(token, "startnr") == 0) {
            read_int(&cfg->startnr, "startnr", curline);
        }
        else if (strcasecmp(token, "atomtype") == 0) {
            read_int(&cfg->atomtype_override, "atomtype", curline);
        }
        else if (strcasecmp(token, "mass") == 0) {
            read_double(&cfg->mass_override[0], "mass", curline);
        }
        else if (strcasecmp(token, "mass1") == 0) {
            read_double(&cfg->mass_override[1], "mass1", curline);
        }
        else if (strcasecmp(token, "lattice_const") == 0) {
            read_double(&cfg->lattice_const, "lattice_const", curline);
        }
        else if (strcasecmp(token, "lattice_const_b") == 0) {
            read_double(&cfg->lattice_const_b, "lattice_const_b", curline);
        }
        else if (strcasecmp(token, "lattice_const_c") == 0) {
            read_double(&cfg->lattice_const_c, "lattice_const_c", curline);
        }
        else if (strcasecmp(token, "box_x") == 0) {
            read_double(&cfg->box_target[0], "box_x", curline);
        }
        else if (strcasecmp(token, "box_y") == 0) {
            read_double(&cfg->box_target[1], "box_y", curline);
        }
        else if (strcasecmp(token, "box_z") == 0) {
            read_double(&cfg->box_target[2], "box_z", curline);
        }
        else if (strcasecmp(token, "shiftvec") == 0) {
            read_vec3(&cfg->shiftvec, "shiftvec", curline);
        }
        else if (strcasecmp(token, "new_x") == 0) {
            read_vec3(&cfg->new_x, "new_x", curline);
            cfg->has_new_x = 1;
        }
        else if (strcasecmp(token, "new_y") == 0) {
            read_vec3(&cfg->new_y, "new_y", curline);
            cfg->has_new_y = 1;
        }
        else if (strcasecmp(token, "new_z") == 0) {
            read_vec3(&cfg->new_z, "new_z", curline);
            cfg->has_new_z = 1;
        }
        else if (strcasecmp(token, "axis1") == 0) {
            read_vec3(&cfg->axis1, "axis1", curline);
        }
        else if (strcasecmp(token, "axis2") == 0) {
            read_vec3(&cfg->axis2, "axis2", curline);
        }
        else if (strcasecmp(token, "angle1") == 0) {
            read_double(&cfg->angle1, "angle1", curline);
        }
        else if (strcasecmp(token, "angle2") == 0) {
            read_double(&cfg->angle2, "angle2", curline);
        }
        else if (strcasecmp(token, "rotmat1_x") == 0) {
            vec3 row;
            read_vec3(&row, "rotmat1_x", curline);
            cfg->rotmat1[0][0] = row.x;
            cfg->rotmat1[0][1] = row.y;
            cfg->rotmat1[0][2] = row.z;
            cfg->rotmat1_given = 1;
        }
        else if (strcasecmp(token, "rotmat1_y") == 0) {
            vec3 row;
            read_vec3(&row, "rotmat1_y", curline);
            cfg->rotmat1[1][0] = row.x;
            cfg->rotmat1[1][1] = row.y;
            cfg->rotmat1[1][2] = row.z;
            cfg->rotmat1_given = 1;
        }
        else if (strcasecmp(token, "rotmat1_z") == 0) {
            vec3 row;
            read_vec3(&row, "rotmat1_z", curline);
            cfg->rotmat1[2][0] = row.x;
            cfg->rotmat1[2][1] = row.y;
            cfg->rotmat1[2][2] = row.z;
            cfg->rotmat1_given = 1;
        }
        else if (strcasecmp(token, "cut_norm1") == 0) {
            read_vec3(&cfg->cut_norm1, "cut_norm1", curline);
        }
        else if (strcasecmp(token, "cut_norm2") == 0) {
            read_vec3(&cfg->cut_norm2, "cut_norm2", curline);
        }
        else if (strcasecmp(token, "plane_distance1") == 0) {
            read_double(&cfg->plane_dist1, "plane_distance1", curline);
            cfg->plane_dist1_given = 1;
        }
        else if (strcasecmp(token, "plane_distance2") == 0) {
            read_double(&cfg->plane_dist2, "plane_distance2", curline);
            cfg->plane_dist2_given = 1;
        }
        else if (strcasecmp(token, "plane_point1") == 0) {
            read_vec3(&cfg->plane_point1, "plane_point1", curline);
            cfg->plane_point1_given = 1;
        }
        else if (strcasecmp(token, "plane_point2") == 0) {
            read_vec3(&cfg->plane_point2, "plane_point2", curline);
            cfg->plane_point2_given = 1;
        }
        else if (strcasecmp(token, "cutdir1") == 0) {
            read_double(&cfg->cutdir1, "cutdir1", curline);
        }
        else if (strcasecmp(token, "cutdir2") == 0) {
            read_double(&cfg->cutdir2, "cutdir2", curline);
        }
        else if (strcasecmp(token, "pos_8d_init") == 0) {
            read_vec3(&cfg->pos_8d_init, "pos_8d_init", curline);
            cfg->has_wyckoff = 1;
        }
        else if (strcasecmp(token, "pos_4c_init") == 0) {
            read_vec3(&cfg->pos_4c_init, "pos_4c_init", curline);
            cfg->has_wyckoff = 1;
        }
        else if (strcasecmp(token, "pos_4c_init1") == 0) {
            read_vec3(&cfg->pos_4c_init1, "pos_4c_init1", curline);
            cfg->has_wyckoff = 1;
        }
        else {
            fprintf(stderr, "WARNING: Unknown parameter '%s' at line %d "
                    "(ignored)\n", token, curline);
        }

    } while (!feof(infile));

    fclose(infile);

    /* Validation */
    if (cfg->outfile[0] == '\0') {
        fprintf(stderr, "Error: 'outfile' must be specified\n");
        exit(2);
    }
    if (cfg->structure_name[0] == '\0') {
        fprintf(stderr, "Error: 'structure' must be specified\n");
        exit(2);
    }
    if (cfg->box_target[0] <= 0 || cfg->box_target[1] <= 0 ||
        cfg->box_target[2] <= 0) {
        fprintf(stderr, "Error: box_x, box_y, box_z must all be positive\n");
        exit(2);
    }

    return 0;
}
