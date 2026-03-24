/* LEGO - lattice.c
 * Lattice operations: rotation matrix construction, M/M^{-1} computation,
 * automatic loop bounds, periodic box sizing.
 */

#include "lego.h"

/* ========== Internal helpers ========== */

/* Check if a value is close to an integer */
static int is_near_int(double v, double tol)
{
    return fabs(v - round(v)) < tol;
}

/* ========== Rotation Matrix Construction ========== */

/* Build the rotation matrix from user specification.
 * Priority: explicit rotmat1 > new_x/y/z directions > axis/angle. */
static void build_rotation(Config *cfg)
{
    mat3 R1, R2;

    if (cfg->rotmat1_given) {
        /* User supplied explicit rotation matrix */
        mat3_copy(R1, cfg->rotmat1);

        /* Check if determinant is valid */
        double d = mat3_det(R1);
        if (fabs(d) < 1.0e-10) {
            fprintf(stderr, "Error: supplied rotation matrix is singular "
                    "(det=%.6e)\n", d);
            exit(2);
        }
    }
    else if ((cfg->has_new_x && cfg->has_new_z) ||
             (cfg->has_new_x && cfg->has_new_y) ||
             (cfg->has_new_y && cfg->has_new_z)) {

        vec3 nx = cfg->new_x;
        vec3 ny = cfg->new_y;
        vec3 nz = cfg->new_z;

        if (cfg->has_new_x && cfg->has_new_z) {
            nx = vec3_normalize(nx);
            nz = vec3_normalize(nz);
            ny = vec3_cross(nz, nx);
            double ny_len = vec3_len(ny);
            if (ny_len < 1.0e-10) {
                fprintf(stderr, "Error: new_x and new_z are parallel\n");
                exit(2);
            }
            ny = vec3_scale(ny, 1.0 / ny_len);
            printf("Computed new_y: [%.6f %.6f %.6f], |y|=%.6f\n",
                   ny.x, ny.y, ny.z, vec3_len(ny));
        }
        else if (cfg->has_new_x && cfg->has_new_y) {
            nx = vec3_normalize(nx);
            ny = vec3_normalize(ny);
            nz = vec3_cross(nx, ny);
            double nz_len = vec3_len(nz);
            if (nz_len < 1.0e-10) {
                fprintf(stderr, "Error: new_x and new_y are parallel\n");
                exit(2);
            }
            nz = vec3_scale(nz, 1.0 / nz_len);
            printf("Computed new_z: [%.6f %.6f %.6f]\n",
                   nz.x, nz.y, nz.z);
        }
        else { /* has_new_y && has_new_z */
            ny = vec3_normalize(ny);
            nz = vec3_normalize(nz);
            nx = vec3_cross(ny, nz);
            double nx_len = vec3_len(nx);
            if (nx_len < 1.0e-10) {
                fprintf(stderr, "Error: new_y and new_z are parallel\n");
                exit(2);
            }
            nx = vec3_scale(nx, 1.0 / nx_len);
            printf("Computed new_x: [%.6f %.6f %.6f]\n",
                   nx.x, nx.y, nx.z);
        }

        /* Rotation matrix: rows = new coordinate axes */
        R1[0][0] = nx.x; R1[0][1] = nx.y; R1[0][2] = nx.z;
        R1[1][0] = ny.x; R1[1][1] = ny.y; R1[1][2] = ny.z;
        R1[2][0] = nz.x; R1[2][1] = nz.y; R1[2][2] = nz.z;
    }
    else {
        /* Rotation from axis1/angle1 */
        mat3_from_axis_angle(R1, cfg->axis1, cfg->angle1);
    }

    /* Validate R1 */
    double d1 = mat3_det(R1);
    if (fabs(fabs(d1) - 1.0) > 1.0e-6) {
        fprintf(stderr, "Error: rotation matrix 1 has det=%.8f "
                "(expected +/-1)\n", d1);
        exit(2);
    }

    mat3_print(R1, "Rotation matrix 1");

    /* Second rotation from axis2/angle2 */
    mat3_from_axis_angle(R2, cfg->axis2, cfg->angle2);

    double d2 = mat3_det(R2);
    if (fabs(fabs(d2) - 1.0) > 1.0e-6) {
        fprintf(stderr, "Error: rotation matrix 2 has det=%.8f "
                "(expected +/-1)\n", d2);
        exit(2);
    }

    mat3_print(R2, "Rotation matrix 2");

    /* Combined rotation: R_total = R2 * R1 */
    mat3_mul_mat(R2, R1, cfg->R_total);
    mat3_print(cfg->R_total, "Combined rotation matrix");

    /* Rotate shift vector */
    cfg->shiftvec_rotated = mat3_mul_vec(cfg->R_total, cfg->shiftvec);
}

/* ========== Lattice Matrix Construction ========== */

/* Build M = R_total * L (rotated lattice matrix).
 * L is the crystal's lattice matrix, optionally scaled by lattice_const. */
static void build_M(Config *cfg)
{
    mat3 L;
    mat3_copy(L, cfg->crystal.lattice);

    /* Apply lattice constant override */
    if (cfg->lattice_const > 0.0) {
        /* Determine current scale from the lattice.
         * For cubic cells, all vectors have same length = lattice parameter.
         * For orthorhombic/hexagonal, we scale uniformly unless
         * lattice_const_b/c are given. */
        double scale = cfg->lattice_const;

        /* Uniform scaling: scale entire lattice by lattice_const / |a1| */
        double a1_len = sqrt(L[0][0]*L[0][0] + L[0][1]*L[0][1] +
                             L[0][2]*L[0][2]);
        if (a1_len > 1.0e-15) {
            double factor = scale / a1_len;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    L[i][j] *= factor;
        }

        /* Non-uniform: override b axis if lattice_const_b given */
        if (cfg->lattice_const_b > 0.0) {
            double a2_len = sqrt(L[1][0]*L[1][0] + L[1][1]*L[1][1] +
                                 L[1][2]*L[1][2]);
            if (a2_len > 1.0e-15) {
                double factor = cfg->lattice_const_b / a2_len;
                /* Only scale this row's current length */
                L[1][0] *= factor;
                L[1][1] *= factor;
                L[1][2] *= factor;
            }
        }

        /* Non-uniform: override c axis if lattice_const_c given */
        if (cfg->lattice_const_c > 0.0) {
            double a3_len = sqrt(L[2][0]*L[2][0] + L[2][1]*L[2][1] +
                                 L[2][2]*L[2][2]);
            if (a3_len > 1.0e-15) {
                double factor = cfg->lattice_const_c / a3_len;
                L[2][0] *= factor;
                L[2][1] *= factor;
                L[2][2] *= factor;
            }
        }
    }

    /* M = R_total * L^T (we need columns of L as lattice vectors,
     * but L is stored as rows. M columns should be rotated lattice vectors.)
     *
     * If L rows are a1, a2, a3, then L^T has columns a1, a2, a3.
     * M = R * L^T means M columns = R*a1, R*a2, R*a3.
     *
     * Then for lattice point (i,j,k): pos = M * (i,j,k)^T
     *   = i * R*a1 + j * R*a2 + k * R*a3. Correct!
     */
    mat3 Lt;
    mat3_transpose(L, Lt);
    mat3_mul_mat(cfg->R_total, Lt, cfg->M);

    printf("Lattice vectors (after scaling):\n");
    printf("  a1 = [%.8f %.8f %.8f]\n", L[0][0], L[0][1], L[0][2]);
    printf("  a2 = [%.8f %.8f %.8f]\n", L[1][0], L[1][1], L[1][2]);
    printf("  a3 = [%.8f %.8f %.8f]\n\n", L[2][0], L[2][1], L[2][2]);

    mat3_print(cfg->M, "Rotated lattice matrix M (columns = rotated vectors)");
}

/* ========== Automatic Loop Bounds ========== */

static void compute_bounds(Config *cfg)
{
    /* Transform 8 box corners via M^{-1} to find (i,j,k) range */
    double bx = cfg->box_actual[0];
    double by = cfg->box_actual[1];
    double bz = cfg->box_actual[2];

    double imin = 1e30, jmin = 1e30, kmin = 1e30;
    double imax = -1e30, jmax = -1e30, kmax = -1e30;

    int ci, cj, ck;
    for (ci = 0; ci <= 1; ci++)
        for (cj = 0; cj <= 1; cj++)
            for (ck = 0; ck <= 1; ck++) {
                vec3 corner = {ci * bx, cj * by, ck * bz};
                vec3 frac = mat3_mul_vec(cfg->Minv, corner);

                if (frac.x < imin) imin = frac.x;
                if (frac.x > imax) imax = frac.x;
                if (frac.y < jmin) jmin = frac.y;
                if (frac.y > jmax) jmax = frac.y;
                if (frac.z < kmin) kmin = frac.z;
                if (frac.z > kmax) kmax = frac.z;
            }

    /* Add margin for basis offsets (basis coords in [0,1)) and safety */
    int margin = 2;
    cfg->ijk_min[0] = (int)floor(imin) - margin;
    cfg->ijk_min[1] = (int)floor(jmin) - margin;
    cfg->ijk_min[2] = (int)floor(kmin) - margin;
    cfg->ijk_max[0] = (int)ceil(imax) + margin;
    cfg->ijk_max[1] = (int)ceil(jmax) + margin;
    cfg->ijk_max[2] = (int)ceil(kmax) + margin;

    long ni = cfg->ijk_max[0] - cfg->ijk_min[0] + 1;
    long nj = cfg->ijk_max[1] - cfg->ijk_min[1] + 1;
    long nk = cfg->ijk_max[2] - cfg->ijk_min[2] + 1;

    printf("Loop bounds (automatic):\n");
    printf("  i: [%d, %d]  (%ld)\n", cfg->ijk_min[0], cfg->ijk_max[0], ni);
    printf("  j: [%d, %d]  (%ld)\n", cfg->ijk_min[1], cfg->ijk_max[1], nj);
    printf("  k: [%d, %d]  (%ld)\n", cfg->ijk_min[2], cfg->ijk_max[2], nk);
    printf("  Total iterations: %ld (x %d basis atoms = %ld candidates)\n\n",
           ni * nj * nk, cfg->crystal.natoms,
           ni * nj * nk * cfg->crystal.natoms);
}

/* ========== Min Repeat Distance ========== */

/* Find the minimum repeat distance along axis d (0=x, 1=y, 2=z).
 *
 * A lattice translation purely along axis d requires:
 *   M * (i,j,k)^T = (L, 0, 0)  for d=0 (etc.)
 *   => (i,j,k)^T = M^{-1} * (L, 0, 0) = L * column_d(M^{-1})
 *
 * For (i,j,k) to be integers, L * M^{-1}[r][d] must be integer for all r.
 * We search for the smallest positive L satisfying this.
 */
static double find_repeat_distance(const Config *cfg, int d)
{
    double col[3]; /* d-th column of M^{-1} */
    int r;

    for (r = 0; r < 3; r++)
        col[r] = cfg->Minv[r][d];

    /* Find non-zero entries */
    int n_nonzero = 0;
    double nonzero_vals[3];
    for (r = 0; r < 3; r++) {
        if (fabs(col[r]) > 1.0e-12) {
            nonzero_vals[n_nonzero++] = col[r];
        }
    }

    if (n_nonzero == 0) {
        fprintf(stderr, "WARNING: M^{-1} column %d is zero — "
                "no periodicity along this axis\n", d);
        return cfg->box_actual[d]; /* no periodicity */
    }

    /* Search for smallest L>0 such that L * col[r] is integer for all r */
    /* Start with the period of the first non-zero entry and try multiples */
    double base = 1.0 / fabs(nonzero_vals[0]);

    for (int n = 1; n <= MAX_REPEAT_SEARCH; n++) {
        double L = n * base;
        int all_int = 1;

        for (r = 0; r < 3; r++) {
            if (fabs(col[r]) > 1.0e-12) {
                if (!is_near_int(L * col[r], INT_TOL)) {
                    all_int = 0;
                    break;
                }
            }
        }

        if (all_int)
            return L;
    }

    /* Try from other non-zero entries as base */
    for (int start = 1; start < n_nonzero; start++) {
        base = 1.0 / fabs(nonzero_vals[start]);
        for (int n = 1; n <= MAX_REPEAT_SEARCH; n++) {
            double L = n * base;
            int all_int = 1;

            for (r = 0; r < 3; r++) {
                if (fabs(col[r]) > 1.0e-12) {
                    if (!is_near_int(L * col[r], INT_TOL)) {
                        all_int = 0;
                        break;
                    }
                }
            }

            if (all_int)
                return L;
        }
    }

    fprintf(stderr, "WARNING: Could not find repeat distance for axis %d "
            "within %d iterations\n", d, MAX_REPEAT_SEARCH);
    return cfg->box_actual[d];
}

static void compute_repeat_distances(Config *cfg)
{
    for (int d = 0; d < 3; d++) {
        cfg->repeat_dist[d] = find_repeat_distance(cfg, d);
    }

    printf("Min repeat distances:\n");
    printf("  x: %.16f\n", cfg->repeat_dist[0]);
    printf("  y: %.16f\n", cfg->repeat_dist[1]);
    printf("  z: %.16f\n\n", cfg->repeat_dist[2]);
}

/* ========== Periodic Box Sizing ========== */

static void adjust_box_periodic(Config *cfg)
{
    for (int d = 0; d < 3; d++) {
        if (cfg->box_mode == MODE_PERIODIC) {
            int n = (int)floor(cfg->box_target[d] / cfg->repeat_dist[d]);
            if (n < 1) n = 1;
            cfg->box_actual[d] = n * cfg->repeat_dist[d] - cfg->epsilon;

            printf("  Axis %c: target=%.4f, repeat=%.8f, n=%d, "
                   "actual=%.16f\n",
                   'x' + d, cfg->box_target[d], cfg->repeat_dist[d],
                   n, cfg->box_actual[d]);
        } else {
            cfg->box_actual[d] = cfg->box_target[d];
        }
    }
    printf("\n");
}

/* ========== Cut Plane Setup ========== */

static void setup_cut_planes(Config *cfg)
{
    /* Normalize cut plane normals */
    double n1len = vec3_len(cfg->cut_norm1);
    double n2len = vec3_len(cfg->cut_norm2);

    if (n1len > 1.0e-15)
        cfg->cut_norm1 = vec3_normalize(cfg->cut_norm1);
    if (n2len > 1.0e-15)
        cfg->cut_norm2 = vec3_normalize(cfg->cut_norm2);

    /* Compute plane distances from points if not given directly */
    if (!cfg->plane_dist1_given && cfg->plane_point1_given) {
        cfg->plane_dist1 = vec3_dot(cfg->plane_point1, cfg->cut_norm1);
    }
    if (!cfg->plane_dist2_given && cfg->plane_point2_given) {
        cfg->plane_dist2 = vec3_dot(cfg->plane_point2, cfg->cut_norm2);
    }
}

/* ========== Main Setup Function ========== */

int lattice_setup(Config *cfg)
{
    printf("=== Lattice Setup ===\n\n");

    /* Apply mass overrides */
    for (int t = 0; t < cfg->crystal.ntypes; t++) {
        if (cfg->mass_override[t] >= 0.0) {
            printf("Mass override type %d: %.6f -> %.6f AMU\n",
                   t, cfg->crystal.masses[t], cfg->mass_override[t]);
            cfg->crystal.masses[t] = cfg->mass_override[t];
        }
    }

    /* Build rotation matrix */
    build_rotation(cfg);

    /* Build M = R * L^T */
    build_M(cfg);

    /* Compute M^{-1} */
    if (mat3_inverse(cfg->M, cfg->Minv) != 0) {
        fprintf(stderr, "Error: rotated lattice matrix M is singular\n");
        exit(2);
    }
    mat3_print(cfg->Minv, "M^{-1}");

    /* Compute repeat distances (needs M^{-1}) */
    compute_repeat_distances(cfg);

    /* Adjust box dimensions for periodicity */
    printf("Box sizing (%s mode):\n",
           cfg->box_mode == MODE_PERIODIC ? "periodic" : "exact");

    /* For periodic mode, we need repeat distances first, then adjust box.
     * For exact mode, box_actual = box_target. */
    if (cfg->box_mode == MODE_EXACT) {
        for (int d = 0; d < 3; d++)
            cfg->box_actual[d] = cfg->box_target[d];
        printf("  Using exact box dimensions as specified.\n\n");
    } else {
        adjust_box_periodic(cfg);
    }

    printf("Final box dimensions: %.16f  %.16f  %.16f\n\n",
           cfg->box_actual[0], cfg->box_actual[1], cfg->box_actual[2]);

    /* Compute automatic loop bounds (needs box_actual and M^{-1}) */
    compute_bounds(cfg);

    /* Set up cut planes */
    setup_cut_planes(cfg);

    printf("Crystal structure: %s\n", cfg->crystal.comment);
    printf("  %d types, %d atoms per cell\n",
           cfg->crystal.ntypes, cfg->crystal.natoms);
    for (int t = 0; t < cfg->crystal.ntypes; t++) {
        printf("  Type %d (%s): %d atoms, mass %.6f AMU\n",
               t, cfg->crystal.species[t], cfg->crystal.counts[t],
               cfg->crystal.masses[t]);
    }
    printf("\n");

    if (cfg->cutdir1 != 0.0) {
        printf("Cut plane 1: norm=[%.4f %.4f %.4f], dist=%.4f, dir=%.1f\n",
               cfg->cut_norm1.x, cfg->cut_norm1.y, cfg->cut_norm1.z,
               cfg->plane_dist1, cfg->cutdir1);
    }
    if (cfg->cutdir2 != 0.0) {
        printf("Cut plane 2: norm=[%.4f %.4f %.4f], dist=%.4f, dir=%.1f\n",
               cfg->cut_norm2.x, cfg->cut_norm2.y, cfg->cut_norm2.z,
               cfg->plane_dist2, cfg->cutdir2);
    }

    printf("=== End Lattice Setup ===\n\n");
    return 0;
}
