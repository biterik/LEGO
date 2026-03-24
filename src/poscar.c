/* LEGO - poscar.c
 * VASP POSCAR file reader and crystal structure loader.
 */

#include "lego.h"
#include <ctype.h>
#include <libgen.h>
#include <unistd.h>

/* Default masses for common elements (AMU).
 * Used when POSCAR species are recognized. */
static double lookup_mass(const char *symbol)
{
    /* Common elements in materials science */
    if (strcasecmp(symbol, "H")  == 0) return 1.008;
    if (strcasecmp(symbol, "He") == 0) return 4.003;
    if (strcasecmp(symbol, "Li") == 0) return 6.941;
    if (strcasecmp(symbol, "Be") == 0) return 9.012;
    if (strcasecmp(symbol, "B")  == 0) return 10.81;
    if (strcasecmp(symbol, "C")  == 0) return 12.011;
    if (strcasecmp(symbol, "N")  == 0) return 14.007;
    if (strcasecmp(symbol, "O")  == 0) return 15.999;
    if (strcasecmp(symbol, "F")  == 0) return 18.998;
    if (strcasecmp(symbol, "Na") == 0) return 22.990;
    if (strcasecmp(symbol, "Mg") == 0) return 24.305;
    if (strcasecmp(symbol, "Al") == 0) return 26.982;
    if (strcasecmp(symbol, "Si") == 0) return 28.086;
    if (strcasecmp(symbol, "P")  == 0) return 30.974;
    if (strcasecmp(symbol, "S")  == 0) return 32.06;
    if (strcasecmp(symbol, "Cl") == 0) return 35.45;
    if (strcasecmp(symbol, "K")  == 0) return 39.098;
    if (strcasecmp(symbol, "Ca") == 0) return 40.078;
    if (strcasecmp(symbol, "Ti") == 0) return 47.867;
    if (strcasecmp(symbol, "V")  == 0) return 50.942;
    if (strcasecmp(symbol, "Cr") == 0) return 51.996;
    if (strcasecmp(symbol, "Mn") == 0) return 54.938;
    if (strcasecmp(symbol, "Fe") == 0) return 55.845;
    if (strcasecmp(symbol, "Co") == 0) return 58.933;
    if (strcasecmp(symbol, "Ni") == 0) return 58.693;
    if (strcasecmp(symbol, "Cu") == 0) return 63.546;
    if (strcasecmp(symbol, "Zn") == 0) return 65.38;
    if (strcasecmp(symbol, "Ga") == 0) return 69.723;
    if (strcasecmp(symbol, "Ge") == 0) return 72.63;
    if (strcasecmp(symbol, "Zr") == 0) return 91.224;
    if (strcasecmp(symbol, "Nb") == 0) return 92.906;
    if (strcasecmp(symbol, "Mo") == 0) return 95.95;
    if (strcasecmp(symbol, "Ag") == 0) return 107.868;
    if (strcasecmp(symbol, "Sn") == 0) return 118.71;
    if (strcasecmp(symbol, "Ta") == 0) return 180.948;
    if (strcasecmp(symbol, "W")  == 0) return 183.84;
    if (strcasecmp(symbol, "Pt") == 0) return 195.084;
    if (strcasecmp(symbol, "Au") == 0) return 196.967;
    if (strcasecmp(symbol, "Pb") == 0) return 207.2;

    /* Generic placeholder 'X' or unknown */
    if (strcasecmp(symbol, "X") == 0) return 1.0;

    fprintf(stderr, "WARNING: Unknown element '%s', using mass 1.0 AMU\n",
            symbol);
    return 1.0;
}

/* Check if a string starts with a letter (to detect element symbols line) */
static int starts_with_letter(const char *s)
{
    while (*s && isspace(*s)) s++;
    return isalpha(*s);
}

/* ========== POSCAR Reader ========== */

int poscar_read(const char *filename, Crystal *cryst)
{
    FILE *f;
    char line[1024];
    int i, atom_idx;

    memset(cryst, 0, sizeof(Crystal));

    f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "Can't open POSCAR file: %s\n", filename);
        return -1;
    }

    /* Line 1: Comment */
    if (!fgets(line, sizeof(line), f)) goto err;
    line[strcspn(line, "\n")] = '\0';
    strncpy(cryst->comment, line, sizeof(cryst->comment) - 1);

    /* Line 2: Scale factor */
    if (!fgets(line, sizeof(line), f)) goto err;
    cryst->scale = atof(line);
    if (cryst->scale <= 0.0) {
        fprintf(stderr, "Invalid POSCAR scale factor: %f\n", cryst->scale);
        goto err;
    }

    /* Lines 3-5: Lattice vectors */
    for (i = 0; i < 3; i++) {
        if (!fgets(line, sizeof(line), f)) goto err;
        if (sscanf(line, "%lf %lf %lf",
                   &cryst->lattice[i][0],
                   &cryst->lattice[i][1],
                   &cryst->lattice[i][2]) != 3)
            goto err;
        /* Apply scale factor */
        cryst->lattice[i][0] *= cryst->scale;
        cryst->lattice[i][1] *= cryst->scale;
        cryst->lattice[i][2] *= cryst->scale;
    }

    /* Line 6: Element names (VASP 5+) or atom counts (VASP 4) */
    if (!fgets(line, sizeof(line), f)) goto err;

    if (starts_with_letter(line)) {
        /* Parse element symbols */
        char *tok = strtok(line, " \t\n");
        cryst->ntypes = 0;
        while (tok && cryst->ntypes < MAX_TYPES) {
            strncpy(cryst->species[cryst->ntypes], tok, 7);
            cryst->species[cryst->ntypes][7] = '\0';
            cryst->masses[cryst->ntypes] = lookup_mass(tok);
            cryst->ntypes++;
            tok = strtok(NULL, " \t\n");
        }

        /* Next line: atom counts */
        if (!fgets(line, sizeof(line), f)) goto err;
    } else {
        /* No element names line (VASP 4 format) */
        /* We'll count types from the counts line */
        cryst->ntypes = 0;
    }

    /* Parse atom counts */
    {
        char *tok;
        char linecopy[1024];
        strncpy(linecopy, line, sizeof(linecopy));
        linecopy[sizeof(linecopy) - 1] = '\0';

        if (cryst->ntypes == 0) {
            /* Count how many numbers we have to determine ntypes */
            tok = strtok(linecopy, " \t\n");
            while (tok && cryst->ntypes < MAX_TYPES) {
                cryst->counts[cryst->ntypes] = atoi(tok);
                snprintf(cryst->species[cryst->ntypes], 8, "T%d",
                         cryst->ntypes + 1);
                cryst->masses[cryst->ntypes] = 1.0;
                cryst->ntypes++;
                tok = strtok(NULL, " \t\n");
            }
        } else {
            tok = strtok(linecopy, " \t\n");
            for (i = 0; i < cryst->ntypes && tok; i++) {
                cryst->counts[i] = atoi(tok);
                tok = strtok(NULL, " \t\n");
            }
        }
    }

    /* Compute total atom count */
    cryst->natoms = 0;
    for (i = 0; i < cryst->ntypes; i++)
        cryst->natoms += cryst->counts[i];

    if (cryst->natoms > MAX_BASIS) {
        fprintf(stderr, "POSCAR has %d atoms, exceeds MAX_BASIS=%d\n",
                cryst->natoms, MAX_BASIS);
        goto err;
    }
    if (cryst->natoms <= 0) {
        fprintf(stderr, "POSCAR has no atoms\n");
        goto err;
    }

    /* Next line: Selective dynamics (optional) or coordinate type */
    if (!fgets(line, sizeof(line), f)) goto err;

    /* Skip 'Selective dynamics' line if present */
    {
        char *p = line;
        while (*p && isspace(*p)) p++;
        if (*p == 'S' || *p == 's') {
            /* Selective dynamics line — read the next line for coord type */
            if (!fgets(line, sizeof(line), f)) goto err;
        }
    }

    /* Coordinate type: Direct or Cartesian */
    int is_cartesian = 0;
    {
        char *p = line;
        while (*p && isspace(*p)) p++;
        if (*p == 'C' || *p == 'c' || *p == 'K' || *p == 'k')
            is_cartesian = 1;
        /* else: Direct (fractional) */
    }

    /* Read atom coordinates */
    atom_idx = 0;
    for (i = 0; i < cryst->ntypes; i++) {
        int j;
        for (j = 0; j < cryst->counts[i]; j++) {
            if (!fgets(line, sizeof(line), f)) goto err;
            double c1, c2, c3;
            if (sscanf(line, "%lf %lf %lf", &c1, &c2, &c3) != 3) goto err;

            if (is_cartesian) {
                /* Convert Cartesian to fractional using lattice inverse */
                mat3 lat_inv;
                if (mat3_inverse(cryst->lattice, lat_inv) != 0) {
                    fprintf(stderr, "Singular lattice matrix in POSCAR\n");
                    goto err;
                }
                vec3 cart = {c1 * cryst->scale, c2 * cryst->scale,
                             c3 * cryst->scale};
                vec3 frac = mat3_mul_vec(lat_inv, cart);
                cryst->basis[atom_idx].x = frac.x;
                cryst->basis[atom_idx].y = frac.y;
                cryst->basis[atom_idx].z = frac.z;
            } else {
                cryst->basis[atom_idx].x = c1;
                cryst->basis[atom_idx].y = c2;
                cryst->basis[atom_idx].z = c3;
            }
            cryst->basis_type[atom_idx] = i;
            atom_idx++;
        }
    }

    fclose(f);
    printf("Read POSCAR: '%s' (%d atoms, %d types)\n",
           cryst->comment, cryst->natoms, cryst->ntypes);
    return 0;

err:
    fclose(f);
    fprintf(stderr, "Error reading POSCAR file: %s\n", filename);
    return -1;
}

/* ========== Crystal Structure Loader ========== */

/* Map built-in structure names to POSCAR filenames */
static const struct {
    const char *name;
    const char *poscar;
} builtin_structures[] = {
    {"sc",    "sc.poscar"},
    {"fcc",   "fcc.poscar"},
    {"bcc",   "bcc.poscar"},
    {"hex",   "hex.poscar"},
    {"hcp",   "hex.poscar"},     /* alias */
    {"l12",   "L12_Ni3Al.poscar"},
    {"b2",    "B2_NiAl.poscar"},
    {"do11",  "D011_Fe3C.poscar"},
    {"d011",  "D011_Fe3C.poscar"}, /* alias */
    {NULL, NULL}
};

/* Try to find a POSCAR file in various search paths */
static int find_poscar(const char *poscar_name, char *result, int maxlen,
                       const char *argv0)
{
    char path[1024];

    /* 1. Check LEGO_STRUCTURES environment variable */
    const char *env = getenv("LEGO_STRUCTURES");
    if (env) {
        snprintf(path, sizeof(path), "%s/%s", env, poscar_name);
        if (access(path, R_OK) == 0) {
            strncpy(result, path, maxlen - 1);
            result[maxlen - 1] = '\0';
            return 0;
        }
    }

    /* 2. Check compile-time STRUCTURES_DIR relative to CWD */
    snprintf(path, sizeof(path), "%s/%s", STRUCTURES_DIR, poscar_name);
    if (access(path, R_OK) == 0) {
        strncpy(result, path, maxlen - 1);
        result[maxlen - 1] = '\0';
        return 0;
    }

    /* 3. Check relative to executable */
    if (argv0) {
        char exe_copy[1024];
        strncpy(exe_copy, argv0, sizeof(exe_copy) - 1);
        exe_copy[sizeof(exe_copy) - 1] = '\0';
        char *dir = dirname(exe_copy);
        snprintf(path, sizeof(path), "%s/structures/%s", dir, poscar_name);
        if (access(path, R_OK) == 0) {
            strncpy(result, path, maxlen - 1);
            result[maxlen - 1] = '\0';
            return 0;
        }
    }

    return -1;
}

/* Global argv0 for find_poscar (set by main) */
static const char *g_argv0 = NULL;

void poscar_set_argv0(const char *argv0)
{
    g_argv0 = argv0;
}

int crystal_load(Config *cfg)
{
    /* Check if structure is "poscar" (user-supplied file) */
    if (strcasecmp(cfg->structure_name, "poscar") == 0) {
        if (cfg->structure_file[0] == '\0') {
            fprintf(stderr, "Error: 'structure_file' required when "
                    "structure=poscar\n");
            exit(2);
        }
        if (poscar_read(cfg->structure_file, &cfg->crystal) != 0)
            exit(2);
        return 0;
    }

    /* Look up built-in structure name */
    int i;
    const char *poscar_name = NULL;
    for (i = 0; builtin_structures[i].name; i++) {
        if (strcasecmp(cfg->structure_name,
                       builtin_structures[i].name) == 0) {
            poscar_name = builtin_structures[i].poscar;
            break;
        }
    }

    if (!poscar_name) {
        fprintf(stderr, "Unknown structure '%s'.\n"
                "Built-in: sc, fcc, bcc, hex, hcp, L12, B2, D011\n"
                "Or use: structure poscar\n"
                "        structure_file path/to/file.poscar\n",
                cfg->structure_name);
        exit(2);
    }

    /* Find the POSCAR file */
    char poscar_path[1024];
    if (find_poscar(poscar_name, poscar_path, sizeof(poscar_path),
                    g_argv0) != 0) {
        fprintf(stderr, "Can't find built-in POSCAR '%s'.\n"
                "Set LEGO_STRUCTURES environment variable or ensure "
                "'structures/' directory is accessible.\n", poscar_name);
        exit(2);
    }

    if (poscar_read(poscar_path, &cfg->crystal) != 0)
        exit(2);

    return 0;
}
