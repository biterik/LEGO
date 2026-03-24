/* LEGO - math3d.c
 * Vector and matrix math for 3D operations.
 * Based on mk_config_calc.c by Erik Bitzek.
 */

#include "lego.h"

/* ========== Vector Operations ========== */

double vec3_len(vec3 v)
{
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

vec3 vec3_add(vec3 a, vec3 b)
{
    vec3 r = {a.x + b.x, a.y + b.y, a.z + b.z};
    return r;
}

vec3 vec3_sub(vec3 a, vec3 b)
{
    vec3 r = {a.x - b.x, a.y - b.y, a.z - b.z};
    return r;
}

vec3 vec3_scale(vec3 v, double s)
{
    vec3 r = {v.x * s, v.y * s, v.z * s};
    return r;
}

vec3 vec3_normalize(vec3 v)
{
    double len = vec3_len(v);
    if (len < 1.0e-15) {
        vec3 r = {0, 0, 0};
        return r;
    }
    return vec3_scale(v, 1.0 / len);
}

double vec3_dot(vec3 a, vec3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

vec3 vec3_cross(vec3 a, vec3 b)
{
    vec3 r;
    r.x = a.y * b.z - a.z * b.y;
    r.y = a.z * b.x - a.x * b.z;
    r.z = a.x * b.y - a.y * b.x;
    return r;
}

/* ========== Matrix Operations ========== */

void mat3_identity(mat3 m)
{
    int i, j;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            m[i][j] = (i == j) ? 1.0 : 0.0;
}

void mat3_zero(mat3 m)
{
    int i, j;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            m[i][j] = 0.0;
}

void mat3_copy(mat3 dst, const mat3 src)
{
    int i, j;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            dst[i][j] = src[i][j];
}

/* Determinant using Sarrus rule */
double mat3_det(const mat3 m)
{
    return m[0][0] * m[1][1] * m[2][2]
         + m[1][0] * m[2][1] * m[0][2]
         + m[2][0] * m[0][1] * m[1][2]
         - m[0][2] * m[1][1] * m[2][0]
         - m[1][2] * m[2][1] * m[0][0]
         - m[2][2] * m[0][1] * m[1][0];
}

/* Matrix inverse using cofactor/adjugate method.
 * Returns 0 on success, -1 if singular. */
int mat3_inverse(const mat3 m, mat3 inv)
{
    double d = mat3_det(m);
    if (fabs(d) < 1.0e-15)
        return -1;

    double inv_det = 1.0 / d;

    /* Cofactor matrix, transposed (adjugate) */
    inv[0][0] =  (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det;
    inv[0][1] = -(m[0][1] * m[2][2] - m[0][2] * m[2][1]) * inv_det;
    inv[0][2] =  (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det;

    inv[1][0] = -(m[1][0] * m[2][2] - m[1][2] * m[2][0]) * inv_det;
    inv[1][1] =  (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det;
    inv[1][2] = -(m[0][0] * m[1][2] - m[0][2] * m[1][0]) * inv_det;

    inv[2][0] =  (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det;
    inv[2][1] = -(m[0][0] * m[2][1] - m[0][1] * m[2][0]) * inv_det;
    inv[2][2] =  (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det;

    return 0;
}

/* Matrix-vector product: result = m * v */
vec3 mat3_mul_vec(const mat3 m, vec3 v)
{
    vec3 r;
    r.x = m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z;
    r.y = m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z;
    r.z = m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z;
    return r;
}

/* Matrix-matrix product: result = a * b
 * result must not alias a or b */
void mat3_mul_mat(const mat3 a, const mat3 b, mat3 result)
{
    int i, j, k;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++) {
            result[i][j] = 0.0;
            for (k = 0; k < 3; k++)
                result[i][j] += a[i][k] * b[k][j];
        }
}

/* Rotation matrix from axis and angle (Rodrigues' formula).
 * Axis is normalized internally. Angle in radians. */
void mat3_from_axis_angle(mat3 m, vec3 axis, double angle)
{
    double len = vec3_len(axis);
    double ax = axis.x / len;
    double ay = axis.y / len;
    double az = axis.z / len;
    double co = cos(angle);
    double si = sin(angle);
    double t = 1.0 - co;

    m[0][0] = t * ax * ax + co;
    m[0][1] = t * ax * ay - az * si;
    m[0][2] = t * ax * az + ay * si;

    m[1][0] = t * ax * ay + az * si;
    m[1][1] = t * ay * ay + co;
    m[1][2] = t * ay * az - ax * si;

    m[2][0] = t * ax * az - ay * si;
    m[2][1] = t * ay * az + ax * si;
    m[2][2] = t * az * az + co;
}

void mat3_transpose(const mat3 m, mat3 result)
{
    int i, j;
    /* Handle aliasing by using temp */
    mat3 tmp;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            tmp[i][j] = m[j][i];
    mat3_copy(result, tmp);
}

void mat3_print(const mat3 m, const char *label)
{
    int i;
    printf("%s:\n", label);
    for (i = 0; i < 3; i++)
        printf("  %12.8f %12.8f %12.8f\n", m[i][0], m[i][1], m[i][2]);
    printf("\n");
}
