/* Math routines for the rotation of vektors   *
 * in mk_conf.c, simple and stupid realisation
 * Matrix[Zeile][Spalte]
 */

#include <math.h>
#include "mk_config_globals.h" 

/* absolute value of a vector */
double vabs(vektor a)
{
  double res;
  res = sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
  return(res);
}

/* Determinant of 3x3 Matrix */
double det(double m[3][3])
{
  double erg;
  erg =m[0][0]*m[1][1]*m[2][2]
    +  m[1][0]*m[2][1]*m[0][2]
    +  m[2][0]*m[0][1]*m[1][2]
    -  m[0][2]*m[1][1]*m[2][0]
    -  m[1][2]*m[2][1]*m[0][0]
    -  m[2][2]*m[0][1]*m[1][0];
  return(erg);
}

/* create rotation matrix from a given axis and angle */

void creatematrix(double matrix[3][3], vektor achse, double winkel)
{
  double betrag = vabs(achse);
  double a = achse.x / betrag;
  double b = achse.y / betrag;
  double c = achse.z / betrag;
  double co = cos(winkel);
  double si = sin(winkel);
  /* correct ? or transposed ?*/
  matrix[0][0] = (1-co)*a*a + co;
  matrix[1][0] = (1-co)*a*b + c*si;
  matrix[2][0] = (1-co)*a*c - b*si;

  matrix[0][1] = (1-co)*a*b - c*si;
  matrix[1][1] = (1-co)*b*b + co;
  matrix[2][1] = (1-co)*b*c + a*si;

  matrix[0][2] = (1-co)*a*c + b*si;
  matrix[1][2] = (1-co)*b*c - a*si;
  matrix[2][2] = (1-co)*c*c + co;

}

/* Multiplikation Matrix - Vektor for 3D */

void mvmult(double matrix[3][3],vektor v, vektor *res)
{
  double c[3],a[3];
  int i,j;
  
  a[0]=v.x; 
  a[1]=v.y; 
  a[2]=v.z;
  
 /*  printf("a: %f %f %f\n",a[0],a[1],a[2]); */
/*   printf("matrix:\n"); */
/*   for(i=0;i<3;i++) */
/*     printf(" %f %f %f\n",matrix[i][0],matrix[i][1],matrix[i][2]); */
/*   printf("\n"); */

  for (i=0;i<3;i++)
    c[i]=0.0;

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      c[i] += a[j]*matrix[i][j];
  
  res->x = c[0];
  res->y = c[1];
  res->z = c[2];

 }


/* Scalar Multiplikation Vektor - Vektor for 3D */
double smult(vektor a, vektor b)
{
  double res;
  res = (a.x*b.x) + (a.y*b.y) +(a.z*b.z);
  return(res);
}

/* Vector product, quick and dirty */
void cross(vektor a, vektor b, vektor *res)
{
  double x,y,z;
  
  x = a.y*b.z - a.z*b.y;
  y = a.z*b.x - a.x*b.z;
  z = a.x*b.y - a.y*b.x;
  
  res->x = x;
  res->y = y;
  res->z = z;

}
