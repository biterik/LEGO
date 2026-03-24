/* mk_config_globals.h
   Global variables for mk_config.c */

#include <stdio.h>

typedef double real;
typedef struct {double x; double y; double z; } vektor;
typedef char str255[255];

#define INIT(data) =data
#define nullvektor  { 0.0, 0.0, 0.0 }

str255 struktur;

str255 outfile;
FILE *out;

int lower, upper;
double e;

int start,type;
int nbasisvects;
int basistype[17];

int nr;

vektor scell_min;
vektor scell_max;    /* super cell dimension */

double mass,m,m1,a,b,c;

double box_x,box_y,box_z;

vektor shiftvec;
vektor n1,n2;
vektor point1,point2;
double d1;
double d2;
double cutdir1,cutdir2;

vektor i8d;
vektor i4c;
vektor i4c1;
vektor nx;
vektor ny;
vektor nz;
vektor m1_x,m1_y,m1_z;

double winkel1;
double winkel2;
vektor axe1;
vektor axe2; 
double drehmatrix1[3][3];
double drehmatrix2[3][3];
