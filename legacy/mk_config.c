/* LEGO (Lattice Elemental Geometry Operations)
 * Program to create arbitary rotated and cut crystal structures in IMD format 
 * Erik Bitzek Jan 2000 2011
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "mk_config_globals.h"

#define EPS 1.0e-4         /* epsilon for determination of supercell */
#define MAX_SCELL 100000000000000.0    /* max. supercell dimension               */

/* Function Prototypes */
void getparams(char *inputfilename);
void mvmult(double matrix[3][3],vektor v, vektor *res);
void creatematrix(double matrix[3][3], vektor achse, double winkel);
double vabs(vektor a);
double smult(vektor a, vektor b);
void cross(vektor a,vektor b, vektor *res);
double det(double matrix[3][3]);
int main (int argc, char *argv[]);

void write_atom(vektor p, int typ,double mass);

int main(int argc, char *argv[])
{
  int i,j,k;
  vektor p,tp;
  
  double n1abs,n2abs;
  double nxabs,nyabs,nzabs;

  
  
  int hasbasis =0;

  double to_ps=0.000103650;     /* amu -> eV*psec^2/angstrom^2 */
  
  /* (rotated) translation vectors */
  vektor a1,a2,a3;
  vektor a1_rot INIT(nullvektor);
  vektor a2_rot INIT(nullvektor);
  vektor a3_rot INIT(nullvektor);
  vektor a1_rot2 INIT(nullvektor);
  vektor a2_rot2 INIT(nullvektor);
  vektor a3_rot2 INIT(nullvektor);


/* (rotated) basis vectors */
  vektor b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16;
  vektor b1_rot INIT(nullvektor);
  vektor b2_rot INIT(nullvektor);
  vektor b3_rot INIT(nullvektor);
  vektor b4_rot INIT(nullvektor);
  vektor b5_rot INIT(nullvektor);
  vektor b6_rot INIT(nullvektor);
  vektor b7_rot INIT(nullvektor);
  vektor b8_rot INIT(nullvektor);
  vektor b9_rot INIT(nullvektor);
  vektor b10_rot INIT(nullvektor);
  vektor b11_rot INIT(nullvektor);
  vektor b12_rot INIT(nullvektor);
  vektor b13_rot INIT(nullvektor);
  vektor b14_rot INIT(nullvektor);
  vektor b15_rot INIT(nullvektor);
  vektor b16_rot INIT(nullvektor);
  vektor b1_rot2 INIT(nullvektor);
  vektor b2_rot2 INIT(nullvektor);
  vektor b3_rot2 INIT(nullvektor);
  vektor b4_rot2 INIT(nullvektor);
  vektor b5_rot2 INIT(nullvektor);
  vektor b6_rot2 INIT(nullvektor);
  vektor b7_rot2 INIT(nullvektor);
  vektor b8_rot2 INIT(nullvektor);
  vektor b9_rot2 INIT(nullvektor);
  vektor b10_rot2 INIT(nullvektor);
  vektor b11_rot2 INIT(nullvektor);
  vektor b12_rot2 INIT(nullvektor);
  vektor b13_rot2 INIT(nullvektor);
  vektor b14_rot2 INIT(nullvektor);
  vektor b15_rot2 INIT(nullvektor);
  vektor b16_rot2 INIT(nullvektor);

  
   /* (rotated) standard vectors */
  vektor e1,e2,e3;
  vektor e1_rot INIT(nullvektor);
  vektor e2_rot INIT(nullvektor);
  vektor e3_rot INIT(nullvektor);
  vektor e1_rot2 INIT(nullvektor);
  vektor e2_rot2 INIT(nullvektor);
  vektor e3_rot2 INIT(nullvektor);
  
  vektor sv_rot INIT(nullvektor);
  vektor sv_rot2 INIT(nullvektor);

  
  e1.x = 1.0;
  e1.y = 0;
  e1.z = 0;

  e2.x = 0.0;
  e2.y = 1.0;
  e2.z = 0;
      
  e3.x = 0.0;
  e3.y = 0.0;
  e3.z = 1.0;

  
  /* default values */
  shiftvec.x = 0.0;
  shiftvec.y = 0.0;
  shiftvec.z = 0.0;
  type=0;
  c=0.0;
  
  nx.x=0.0;
  nx.y=0.0;
  nx.z=0.0;
  ny.x=0.0;
  ny.y=0.0;
  ny.z=0.0;
  nz.x=0.0;
  nz.y=0.0;
  nz.z=0.0;

  axe1.x =1.0;
  axe1.y =0.0;
  axe1.z =0.0;
  axe2.x =1.0;
  axe2.y =0.0;
  axe2.z =0.0;


  drehmatrix1[0][0] = 0.0;
  drehmatrix1[0][1] = 0.0;
  drehmatrix1[0][2] = 0.0;
  drehmatrix1[1][0] = 0.0;
  drehmatrix1[1][1] = 0.0;
  drehmatrix1[1][2] = 0.0;
  drehmatrix1[2][0] = 0.0;
  drehmatrix1[2][1] = 0.0;
  drehmatrix1[2][2] = 0.0;

  
  d1=-999999.9999;
  d2=-999999.9999;;
   /* get global parameters */
  getparams(argv[1]);
  
  /* preparations for output */
  m *= to_ps;
  nr = start;
  n1abs = vabs(n1);
  n2abs = vabs(n2);
  n1.x *= 1.0/n1abs;
  n1.y *= 1.0/n1abs;
  n1.z *= 1.0/n1abs;
  n2.x *= 1.0/n2abs;
  n2.y *= 1.0/n2abs;
  n2.z *= 1.0/n2abs;

 /* info for hessian normal form */
  if(d1==-999999.9999)
  {
      d1=smult(point1,n1);
  }
  if(d2==-999999.9999)
  {
      d2=smult(point2,n2);
  }
  
  
  /* Informations about input values */
  printf("Crystalstructure: %s\n",struktur);
  printf("lower: %d upper: %d\n",lower,upper);
  printf("Epsilon: %f\n",e);
  printf("outfile: %s\n",outfile);
  printf("Atomnr. starting value: %d\n",start);
  printf("mass (rescaled): %f\n",m);
  printf("a: %f\n",a);
  printf("x -> [ %f %f %f ]\n",nx.x,nx.y,nx.z);
  printf("y -> [ %f %f %f ]\n",ny.x,ny.y,ny.z);
  printf("z -> [ %f %f %f ]\n",nz.x,nz.y,nz.z);
  printf("box_x: %f box_y: %f box_z: %f\n",box_x,box_y,box_z);
  printf("cut_norm1: %f %f %f\n",n1.x,n1.y,n1.z);
  printf("distance cutplane - 0: %f  side to cut away: %f\n",d1,cutdir1);
  printf("cut_norm2 %f %f %f\n",n2.x,n2.y,n2.z);
  printf("distance cutplane - 0: %f  side to cut away: %f\n",d2,cutdir2);
  printf("axis1: %f %f %f\n",axe1.x,axe1.y,axe1.z);
  printf("axis2: %f %f %f\n",axe2.x,axe2.y,axe2.z);
  printf("angle1: %f\n",winkel1);
  printf("angle2: %f\n",winkel2);fflush(stdout);
  
  /* open outfile */
  out = fopen(outfile,"w");
  if (NULL==out){
    printf("Can't open outputfile %s\n",outfile);
    exit(2);
  }

  if (strcasecmp(struktur,"sc")==0)
    {
      /* primitve translations of sc */
        // see Ashcroft
      a1.x = a;       
      a1.y = 0.0;       
      a1.z = 0.0;

      a2.x = 0.0;
      a2.y = a;
      a2.z = 0.0;
      
      a3.x = 0.0;
      a3.y = 0.0;
      a3.z = a;
    }
  
  else if (strcasecmp(struktur,"fcc")==0)
      // diamond = 2*fcc, w. shiftvector a1/4 +a2/4+a3/4
    {
      /* primitve translations of fcc */
        // see Ashcroft
      a1.x = 0.5*a;       
      a1.y = 0.5*a;       
      a1.z = 0.0;

      a2.x = 0.0;
      a2.y = 0.5*a;
      a2.z = 0.5*a;
      
      a3.x = 0.5*a;
      a3.y = 0.0;
      a3.z = 0.5*a;
    }
  else if (strcasecmp(struktur,"bcc")==0)
    {
    /* primitve translations of bcc */
      a1.x = 0.5*a;       
      a1.y = 0.5*a;       
      a1.z =-0.5*a;

      a2.x =-0.5*a;
      a2.y = 0.5*a;
      a2.z = 0.5*a;
      
      a3.x = 0.5*a;
      a3.y =-0.5*a;
      a3.z = 0.5*a;
    }
  else if (strcasecmp(struktur,"hex")==0)
      // hcp = 2*hex, w. shiftvector a1/3 +a2/3+a3/2
  {
      if (c<=0)
      {
          printf("Hey, %s needs more lattice parameters, but c=%f\n",struktur,c);
          exit(2);
      }

      /* primitve translations of hcp */
      a1.x = 1.0*a;       
      a1.y = 0.0;
      a1.z = 0.0;

      a2.x = 0.5*a;
      a2.y = sqrt(3.0)*a/2.0;
      a2.z = 0.0;
      
      a3.x = 0.0;
      a3.y = 0.0;
      a3.z = c;
    }
  else if (strcasecmp(struktur,"L12")==0)
      // e.g. Ni3Al, now we use for the first time basis vectors ;-)
    {
      /* primitive vectors */
      // see http://cst-www.nrl.navy.mil/lattice/struk/l1_2.html
      a1.x = 1.0*a;       
      a1.y = 0.0;       
      a1.z = 0.0;

      a2.x = 0.0;
      a2.y = 1.0*a;
      a2.z = 0.0;
      
      a3.x = 0.0;
      a3.y = 0.0;
      a3.z = 1.0*a;

      m1 *=to_ps;
      printf("structure with 2 atoms:\n 0: mass %f amu =%f 1: mass %f amu=%f \n",m,m/to_ps,m1,m1/to_ps);
      
      hasbasis =1;
      nbasisvects=4;

      // Al
      b1.x = 0.0;       
      b1.y = 0.0;       
      b1.z = 0.0;
      basistype[1]=0; // same numbering as for basisvector

      // Ni

      b2.x = 0.0;
      b2.y = 0.5*a;
      b2.z = 0.5*a;
      basistype[2]=1;

      b3.x = 0.5*a;
      b3.y = 0.0;
      b3.z = 0.5*a;
      basistype[3]=1;
      
      b4.x = 0.5*a;
      b4.y = 0.5*a;
      b4.z = 0.0;
      basistype[4]=1;

    }
  else if (strcasecmp(struktur,"B2")==0)
      // e.g. NiAl
      // adapted from L12 and bcc structure
    {
      /* primitive vectors */
      // see http://cst-www.nrl.navy.mil/lattice/struk/l1_2.html
      a1.x = 1.0*a;       
      a1.y = 0.0;       
      a1.z = 0.0;

      a2.x = 0.0;
      a2.y = 1.0*a;
      a2.z = 0.0;
      
      a3.x = 0.0;
      a3.y = 0.0;
      a3.z = 1.0*a;

      m1 *=to_ps;
      printf("structure with 2 atoms:\n 0: mass %f amu =%f 1: mass %f amu=%f \n",m,m/to_ps,m1,m1/to_ps);

      hasbasis =2;
      nbasisvects=2;

      // Al
      b1.x = 0.0;       
      b1.y = 0.0;       
      b1.z = 0.0;
      basistype[1]=0; // same numbering as for basisvector

      // Ni

      b2.x = 0.5*a;
      b2.y = 0.5*a;
      b2.z = 0.5*a;
      basistype[2]=1;


    /* primitve translations of bcc 
      a1.x = 0.5*a;
      a1.y = 0.5*a;
      a1.z =-0.5*a;

      a2.x =-0.5*a;
      a2.y = 0.5*a;
      a2.z = 0.5*a;

      a3.x = 0.5*a;
      a3.y =-0.5*a;
      a3.z = 0.5*a;
	*/

    }
  else if (strcasecmp(struktur,"DO11")==0)
      // e.g. Fe3C
      // implemented by Jo Moeller, Sun Jan 27 11:40:02 CET 2013
      // 
    {
      /* primitive vectors */
      // see http://cst-www.nrl.navy.mil/lattice/struk/l1_2.html
      a1.x = 1.0*a;       
      a1.y = 0.0;       
      a1.z = 0.0;

      a2.x = 0.0;
      a2.y = 1.0*b;
      a2.z = 0.0;
      
      a3.x = 0.0;
      a3.y = 0.0;
      a3.z = 1.0*c;

      m1 *=to_ps;
      printf("structure with 2 atoms:\n 0: mass %f amu =%f 1: mass %f amu=%f \n",m,m/to_ps,m1,m1/to_ps);
      printf("initial Wyckoff positions:\n 8d(0): %f, %f, %f \n 4c(0): %f, %f, %f \n 4c(1): %f, %f, %f \n",i8d.x,i8d.y,i8d.z,i4c.x,i4c.y,i4c.z,i4c1.x,i4c1.y,i4c1.z);
      
      hasbasis=3;
      nbasisvects=16;

      // the 8 general base vectors of atom type 0, e.g., FeI:
      b1.x = i8d.x*a;       
      b1.y = i8d.y*b;       
      b1.z = i8d.z*c;
      basistype[1]=0;

      b2.x = (1-i8d.x+0.5)*a;
      b2.y = (1-i8d.y)*b;
      b2.z = (i8d.z+0.5)*c;
      basistype[2]=0;

      b3.x = (1-i8d.x)*a;
      b3.y = (i8d.y+0.5)*b;
      b3.z = (1-i8d.z)*c;
      basistype[3]=0;
           
      b4.x = (i8d.x+0.5)*a;
      b4.y = (1-i8d.y+0.5)*b;
      b4.z = (1-i8d.z+0.5)*c;
      basistype[4]=0;

      b5.x = (1-i8d.x)*a;
      b5.y = (1-i8d.y)*b;
      b5.z = (1-i8d.z)*c;
      basistype[5]=0;
			            
      b6.x = (i8d.x+0.5)*a;
      b6.y = (i8d.y)*b;
      b6.z = (1-i8d.z+0.5)*c;
      basistype[6]=0;
				            
      b7.x = i8d.x*a;
      b7.y = (1-i8d.y+0.5)*b;
      b7.z = i8d.z*c;
      basistype[7]=0;

      b8.x = (-i8d.x+0.5)*a;
      b8.y = (i8d.y+0.5)*b;
      b8.z = (i8d.z+0.5)*c;
      basistype[8]=0;

      // the 4 special base vectors of atom type 0, e.g., FeII:
      b9.x = i4c.x*a;
      b9.y = i4c.y*b;
      b9.z = i4c.z*c;
      basistype[9]=0;

      b10.x = (-i4c.x+0.5)*a;
      b10.y = (1-i4c.y)*b;
      b10.z = (i4c.z+0.5)*c;
      basistype[10]=0;
      
      b11.x = (1-i4c.x)*a;       
      b11.y = (1-i4c.y)*b;       
      b11.z = (1-i4c.z)*c;
      basistype[11]=0;

      b12.x = (i4c.x+0.5)*a;
      b12.y = (i4c.y)*b;
      b12.z = (-i4c.z+0.5)*c;
      basistype[12]=0;

      // the 4 special base vectors of atom type 1, e.g., C:
      b13.x = i4c1.x*a;
      b13.y = i4c1.y*b;
      b13.z = i4c1.z*c;
      basistype[13]=1;

      b14.x = (-i4c1.x+0.5)*a;
      b14.y = (1-i4c1.y)*b;
      b14.z = (i4c1.z+0.5)*c;
      basistype[14]=1;
      
      b15.x = (1-i4c1.x)*a;       
      b15.y = (1-i4c1.y)*b;       
      b15.z = (1-i4c1.z)*c;
      basistype[15]=1;

      b16.x = (i4c1.x+0.5)*a;
      b16.y = (i4c1.y)*b;
      b16.z = (-i4c1.z+0.5)*c;
      basistype[16]=1;

    }
  else 
    {
      printf("Unknown Crystalstructure : %s\n",struktur);
      exit(2);
    }

 

  if(det(drehmatrix1)==0.0)
    {
       
        /* creation of transformation matrix */
        nxabs = vabs(nx);
        nyabs = vabs(ny); nzabs = vabs(nz); 
        if ((nxabs!=0.0) && (nzabs!=0.0)){
            nx.x *= 1.0/nxabs;
            nx.y *= 1.0/nxabs;
            nx.z *= 1.0/nxabs;
            nz.x *= 1.0/nzabs;
            nz.y *= 1.0/nzabs;
            nz.z *= 1.0/nzabs; 
            
            cross(nz,nx,&ny);
            printf("y -> [ %f %f %f ], |y|= %f\n",ny.x,ny.y,ny.z,vabs(ny));
            
            drehmatrix1[0][0] = nx.x;
            drehmatrix1[0][1] = nx.y;
            drehmatrix1[0][2] = nx.z;
            drehmatrix1[1][0] = ny.x;
            drehmatrix1[1][1] = ny.y;
            drehmatrix1[1][2] = ny.z;
            drehmatrix1[2][0] = nz.x;
            drehmatrix1[2][1] = nz.y;
            drehmatrix1[2][2] = nz.z;
        }
        else if((nxabs!=0.0) && (nyabs!=0.0)){
            nx.x *= 1.0/nxabs;
            nx.y *= 1.0/nxabs;
            nx.z *= 1.0/nxabs; 
            ny.x *= 1.0/nyabs;
            ny.y *= 1.0/nyabs;
            ny.z *= 1.0/nyabs;
    
            cross(nx,ny,&nz);
            printf("z -> [ %f %f %f ]\n",nz.x,nz.y,nz.z);
            

            drehmatrix1[0][0] = nx.x;
            drehmatrix1[0][1] = nx.y;
            drehmatrix1[0][2] = nx.z;
            drehmatrix1[1][0] = ny.x;
            drehmatrix1[1][1] = ny.y;
            drehmatrix1[1][2] = ny.z;
            drehmatrix1[2][0] = nz.x;
            drehmatrix1[2][1] = nz.y;
            drehmatrix1[2][2] = nz.z;

        }
        else if ((nyabs!=0.0) && (nzabs!=0.0)){
            ny.x *= 1.0/nyabs;
            ny.y *= 1.0/nyabs;
            ny.z *= 1.0/nyabs;
            nz.x *= 1.0/nzabs;
	    nz.z *= 1.0/nzabs;
            
            cross(ny,nz,&nx);
            printf("x -> [ %f %f %f ]\n",nx.x,nx.y,nx.z);
            
            drehmatrix1[0][0] = nx.x;
            drehmatrix1[0][1] = nx.y;
            drehmatrix1[0][2] = nx.z;
            drehmatrix1[1][0] = ny.x;
            drehmatrix1[1][1] = ny.y;
            drehmatrix1[1][2] = ny.z;
            drehmatrix1[2][0] = nz.x;
            drehmatrix1[2][1] = nz.y;
            drehmatrix1[2][2] = nz.z;
            
        }
        else
        {
            /* rotate the primitive Translation Vectors */
            creatematrix(drehmatrix1,axe1,winkel1);
            
        }
    }
   printf("Rotation matrix 1:\n");
   for (i=0;i<3;i++)
     printf(" %f %f %f\n",drehmatrix1[i][0],drehmatrix1[i][1],drehmatrix1[i][2]);
   printf("\n");

  if(fabs(det(drehmatrix1)-1.0) > 1.0e-6)
    {
      printf("Det(drehmatrix1)=%f\n",det(drehmatrix1));
      exit(2);
    }
  creatematrix(drehmatrix2,axe2,winkel2);
  printf("Rotation matrix 2:\n"); 
  for (i=0;i<3;i++) 
    printf(" %f %f %f\n",drehmatrix2[i][0],drehmatrix2[i][1],drehmatrix2[i][2]); 
  printf("\n"); 
  if(fabs(det(drehmatrix2)-1.0) > 1.0e-6)
    {
      printf("Det(drehmatrix2)=%f\n",det(drehmatrix2));
      exit(2);
    }
  mvmult(drehmatrix1,a1,&a1_rot);
  mvmult(drehmatrix1,a2,&a2_rot);
  mvmult(drehmatrix1,a3,&a3_rot);

  /* L12 */
  if (hasbasis==1)
    {
      mvmult(drehmatrix1,b1,&b1_rot);
      mvmult(drehmatrix1,b2,&b2_rot);
      mvmult(drehmatrix1,b3,&b3_rot);
      mvmult(drehmatrix1,b4,&b4_rot);
    }
  
  /* B2: */
  if (hasbasis==2) 
    {
      mvmult(drehmatrix1,b1,&b1_rot);
      mvmult(drehmatrix1,b2,&b2_rot);
    }

  /* D011: */
  if (hasbasis==3)
    {
      mvmult(drehmatrix1,b1,&b1_rot);
      mvmult(drehmatrix1,b2,&b2_rot);
      mvmult(drehmatrix1,b3,&b3_rot);
      mvmult(drehmatrix1,b4,&b4_rot);
      mvmult(drehmatrix1,b5,&b5_rot);
      mvmult(drehmatrix1,b6,&b6_rot);
      mvmult(drehmatrix1,b7,&b7_rot);
      mvmult(drehmatrix1,b8,&b8_rot);
      mvmult(drehmatrix1,b9,&b9_rot);
      mvmult(drehmatrix1,b10,&b10_rot);
      mvmult(drehmatrix1,b11,&b11_rot);
      mvmult(drehmatrix1,b12,&b12_rot);
      mvmult(drehmatrix1,b13,&b13_rot);
      mvmult(drehmatrix1,b14,&b14_rot);
      mvmult(drehmatrix1,b15,&b15_rot);
      mvmult(drehmatrix1,b16,&b16_rot);
    }
/*   printf("1. Rotation:\n"); */
/*   printf(" %f %f %f\n\n",a1_rot.x,a1_rot.y,a1_rot.z); */
/*   printf(" %f %f %f\n\n",a2_rot.x,a2_rot.y,a2_rot.z); */
/*   printf(" %f %f %f\n\n",a3_rot.x,a3_rot.y,a3_rot.z); */

  mvmult(drehmatrix2,a1_rot,&a1_rot2);
  mvmult(drehmatrix2,a2_rot,&a2_rot2);
  mvmult(drehmatrix2,a3_rot,&a3_rot2);

  /* L12: */
  if (hasbasis==1)
    {
      mvmult(drehmatrix2,b1_rot,&b1_rot2);
      mvmult(drehmatrix2,b2_rot,&b2_rot2);
      mvmult(drehmatrix2,b3_rot,&b3_rot2);
      mvmult(drehmatrix2,b4_rot,&b4_rot2);
    }
  /* B2: */
  if (hasbasis==2)
    {
      mvmult(drehmatrix2,b1_rot,&b1_rot2);
      mvmult(drehmatrix2,b2_rot,&b2_rot2);
    }

  /* D011: */
  if (hasbasis==3)
    {
      mvmult(drehmatrix2,b1_rot,&b1_rot2);
      mvmult(drehmatrix2,b2_rot,&b2_rot2);
      mvmult(drehmatrix2,b3_rot,&b3_rot2);
      mvmult(drehmatrix2,b4_rot,&b4_rot2);
      mvmult(drehmatrix2,b5_rot,&b5_rot2);
      mvmult(drehmatrix2,b6_rot,&b6_rot2);
      mvmult(drehmatrix2,b7_rot,&b7_rot2);
      mvmult(drehmatrix2,b8_rot,&b8_rot2);
      mvmult(drehmatrix2,b9_rot,&b9_rot2);
      mvmult(drehmatrix2,b10_rot,&b10_rot2);
      mvmult(drehmatrix2,b11_rot,&b11_rot2);
      mvmult(drehmatrix2,b12_rot,&b12_rot2);
      mvmult(drehmatrix2,b13_rot,&b13_rot2);
      mvmult(drehmatrix2,b14_rot,&b14_rot2);
      mvmult(drehmatrix2,b15_rot,&b15_rot2);
      mvmult(drehmatrix2,b16_rot,&b16_rot2);
    }
  printf("unit vectors after 2. Rotation:\n");
  printf(" %.8f %.8f %.8f\n",a1_rot2.x,a1_rot2.y,a1_rot2.z);
  printf(" %.8f %.8f %.8f\n",a2_rot2.x,a2_rot2.y,a2_rot2.z);
  printf(" %.8f %.8f %.8f\n\n",a3_rot2.x,a3_rot2.y,a3_rot2.z);
  
  /* L12: */
  if (hasbasis==1)
    {
      printf("basis vectors after 2. Rotation:\n");
      printf("type %d: %.8f %.8f %.8f\n",  basistype[1],b1_rot2.x,b1_rot2.y,b1_rot2.z);
      printf("type %d: %.8f %.8f %.8f\n",  basistype[2],b2_rot2.x,b2_rot2.y,b2_rot2.z);
      printf("type %d: %.8f %.8f %.8f\n",  basistype[3],b3_rot2.x,b3_rot2.y,b3_rot2.z);
      printf("type %d: %.8f %.8f %.8f\n\n",  basistype[4],b4_rot2.x,b4_rot2.y,b4_rot2.z);
    }
  /* B2: */
  if (hasbasis==2)
    {
      printf("basis vectors after 2. Rotation:\n");
      printf("type %d: %.8f %.8f %.8f\n",  basistype[1],b1_rot2.x,b1_rot2.y,b1_rot2.z);
      printf("type %d: %.8f %.8f %.8f\n\n",  basistype[2],b2_rot2.x,b2_rot2.y,b2_rot2.z);
  
    }
  /* D011: */
  if (hasbasis==3)
    {
      printf("basis vectors after 2. Rotation:\n");
      printf("type %d: %.8f %.8f %.8f\n",  basistype[1],b1_rot2.x,b1_rot2.y,b1_rot2.z);
      printf("type %d: %.8f %.8f %.8f\n",  basistype[2],b2_rot2.x,b2_rot2.y,b2_rot2.z);
      printf("type %d: %.8f %.8f %.8f\n",  basistype[3],b3_rot2.x,b3_rot2.y,b3_rot2.z);
      printf("type %d: %.8f %.8f %.8f\n",  basistype[4],b4_rot2.x,b4_rot2.y,b4_rot2.z);
      printf("type %d: %.8f %.8f %.8f\n",  basistype[5],b5_rot2.x,b5_rot2.y,b5_rot2.z);
      printf("type %d: %.8f %.8f %.8f\n",  basistype[6],b6_rot2.x,b6_rot2.y,b6_rot2.z);
      printf("type %d: %.8f %.8f %.8f\n",  basistype[7],b7_rot2.x,b7_rot2.y,b7_rot2.z);
      printf("type %d: %.8f %.8f %.8f\n",  basistype[8],b8_rot2.x,b8_rot2.y,b8_rot2.z);
      printf("type %d: %.8f %.8f %.8f\n",  basistype[9],b9_rot2.x,b9_rot2.y,b9_rot2.z);
      printf("type %d: %.8f %.8f %.8f\n",  basistype[10],b10_rot2.x,b10_rot2.y,b10_rot2.z);
      printf("type %d: %.8f %.8f %.8f\n",  basistype[11],b11_rot2.x,b11_rot2.y,b11_rot2.z);
      printf("type %d: %.8f %.8f %.8f\n",  basistype[12],b12_rot2.x,b12_rot2.y,b12_rot2.z);
      printf("type %d: %.8f %.8f %.8f\n",  basistype[13],b13_rot2.x,b13_rot2.y,b13_rot2.z);
      printf("type %d: %.8f %.8f %.8f\n",  basistype[14],b14_rot2.x,b14_rot2.y,b14_rot2.z);
      printf("type %d: %.8f %.8f %.8f\n",  basistype[15],b15_rot2.x,b15_rot2.y,b15_rot2.z);
      printf("type %d: %.8f %.8f %.8f\n\n",basistype[16],b16_rot2.x,b16_rot2.y,b16_rot2.z);
  
    }
  mvmult(drehmatrix1,e1,&e1_rot);
  mvmult(drehmatrix1,e2,&e2_rot);
  mvmult(drehmatrix1,e3,&e3_rot);

  mvmult(drehmatrix2,e1_rot,&e1_rot2);
  mvmult(drehmatrix2,e2_rot,&e2_rot2);
  mvmult(drehmatrix2,e3_rot,&e3_rot2);

  // shift vector for multi-base structures
  /* shiftvec.x *= a; */
  /* shiftvec.y *= a; */
  /* shiftvec.z *= a; */
  mvmult(drehmatrix1,shiftvec,&sv_rot);
  mvmult(drehmatrix2,sv_rot,&sv_rot2);
  
  printf("Rotation matrix :\n");
   // this matrix is used by GBcharacterizer and in matlab
  // this is also how it will be read in 
  printf(" %.12f %.12f %.12f\n",e1_rot2.x,e1_rot2.y,e1_rot2.z);
  printf(" %.12f %.12f %.12f\n",e2_rot2.x,e2_rot2.y,e2_rot2.z);
  printf(" %.12f %.12f %.12f\n\n",e3_rot2.x,e3_rot2.y,e3_rot2.z);
  
  /* printf(" %.12f %.12f %.12f\n",e1_rot2.x,e2_rot2.x,e3_rot2.x); */
  /* printf(" %.12f %.12f %.12f\n",e1_rot2.y,e2_rot2.y,e3_rot2.y); */
  /* printf(" %.12f %.12f %.12f\n\n",e1_rot2.z,e2_rot2.z,e3_rot2.z); */
  
  
  scell_min.x =  MAX_SCELL;
  scell_min.y =  MAX_SCELL;
  scell_min.z =  MAX_SCELL;

  scell_max.x = 0.0;
  scell_max.y = 0.0;
  scell_max.z = 0.0;

 
  /* calculate atom positions */
 
  for(i=lower;i<upper;i++)
    for(j=lower;j<upper;j++)
      for(k=lower;k<upper;k++){
	p.x = i*a1_rot2.x +  j*a2_rot2.x + k*a3_rot2.x;
	p.y = i*a1_rot2.y +  j*a2_rot2.y + k*a3_rot2.y;
	p.z = i*a1_rot2.z +  j*a2_rot2.z + k*a3_rot2.z;

        /* bcc, fcc , sc, hex: */
        if (hasbasis==0) {
          p.x += sv_rot2.x;
          p.y += sv_rot2.y;
          p.z += sv_rot2.z;
          write_atom(p,type,m);
        }

        /* L12: */
	if(hasbasis==1) {
            type = basistype[1];
            tp.x = p.x + b1_rot2.x;
            tp.y = p.y + b1_rot2.y;
            tp.z = p.z + b1_rot2.z;
            write_atom(tp,type,m);

	    type = basistype[2];
	    if(type==0) mass=m; else mass=m1;
	    tp.x = p.x + b2_rot2.x;
	    tp.y = p.y + b2_rot2.y;
	    tp.z = p.z + b2_rot2.z;
	    write_atom(tp,type,mass);

	    type = basistype[3];
	    if(type==0) mass=m; else mass=m1;
	    tp.x = p.x + b3_rot2.x;
	    tp.y = p.y + b3_rot2.y;
	    tp.z = p.z + b3_rot2.z;
	    write_atom(tp,type,mass);

	    type = basistype[4];
	    if(type==0) mass=m; else mass=m1;
	    tp.x = p.x + b4_rot2.x;
	    tp.y = p.y + b4_rot2.y;
	    tp.z = p.z + b4_rot2.z;
	    write_atom(tp,type,mass);
	}

	/* B2: */
        if(hasbasis==2) {
            type = basistype[1];
            tp.x = p.x + b1_rot2.x;
            tp.y = p.y + b1_rot2.y;
            tp.z = p.z + b1_rot2.z;
            write_atom(tp,type,m);

            type = basistype[2];
            if(type==0) mass=m; else mass=m1;
            tp.x = p.x + b2_rot2.x;
            tp.y = p.y + b2_rot2.y;
            tp.z = p.z + b2_rot2.z;
            write_atom(tp,type,mass);
        }

        /* D011: */
	if(hasbasis==3) {
            type = basistype[1];
            tp.x = p.x + b1_rot2.x;
            tp.y = p.y + b1_rot2.y;
            tp.z = p.z + b1_rot2.z;
            write_atom(tp,type,m);

	    type = basistype[2];
	    if(type==0) mass=m; else mass=m1;
	    tp.x = p.x + b2_rot2.x;
	    tp.y = p.y + b2_rot2.y;
	    tp.z = p.z + b2_rot2.z;
	    write_atom(tp,type,mass);

	    type = basistype[3];
	    if(type==0) mass=m; else mass=m1;
	    tp.x = p.x + b3_rot2.x;
	    tp.y = p.y + b3_rot2.y;
	    tp.z = p.z + b3_rot2.z;
	    write_atom(tp,type,mass);

	    type = basistype[4];
	    if(type==0) mass=m; else mass=m1;
	    tp.x = p.x + b4_rot2.x;
	    tp.y = p.y + b4_rot2.y;
	    tp.z = p.z + b4_rot2.z;
	    write_atom(tp,type,mass);

	    type = basistype[5];
	    if(type==0) mass=m; else mass=m1;
	    tp.x = p.x + b5_rot2.x;
	    tp.y = p.y + b5_rot2.y;
	    tp.z = p.z + b5_rot2.z;
	    write_atom(tp,type,mass);

	    type = basistype[6];
	    if(type==0) mass=m; else mass=m1;
	    tp.x = p.x + b6_rot2.x;
	    tp.y = p.y + b6_rot2.y;
	    tp.z = p.z + b6_rot2.z;
	    write_atom(tp,type,mass);

	    type = basistype[7];
	    if(type==0) mass=m; else mass=m1;
	    tp.x = p.x + b7_rot2.x;
	    tp.y = p.y + b7_rot2.y;
	    tp.z = p.z + b7_rot2.z;
	    write_atom(tp,type,mass);

	    type = basistype[8];
	    if(type==0) mass=m; else mass=m1;
	    tp.x = p.x + b8_rot2.x;
	    tp.y = p.y + b8_rot2.y;
	    tp.z = p.z + b8_rot2.z;
	    write_atom(tp,type,mass);

	    type = basistype[9];
	    if(type==0) mass=m; else mass=m1;
	    tp.x = p.x + b9_rot2.x;
	    tp.y = p.y + b9_rot2.y;
	    tp.z = p.z + b9_rot2.z;
	    write_atom(tp,type,mass);

	    type = basistype[10];
	    if(type==0) mass=m; else mass=m1;
	    tp.x = p.x + b10_rot2.x;
	    tp.y = p.y + b10_rot2.y;
	    tp.z = p.z + b10_rot2.z;
	    write_atom(tp,type,mass);

	    type = basistype[11];
	    if(type==0) mass=m; else mass=m1;
	    tp.x = p.x + b11_rot2.x;
	    tp.y = p.y + b11_rot2.y;
	    tp.z = p.z + b11_rot2.z;
	    write_atom(tp,type,mass);

	    type = basistype[12];
	    if(type==0) mass=m; else mass=m1;
	    tp.x = p.x + b12_rot2.x;
	    tp.y = p.y + b12_rot2.y;
	    tp.z = p.z + b12_rot2.z;
	    write_atom(tp,type,mass);

	    type = basistype[13];
	    if(type==0) mass=m; else mass=m1;
	    tp.x = p.x + b13_rot2.x;
	    tp.y = p.y + b13_rot2.y;
	    tp.z = p.z + b13_rot2.z;
	    write_atom(tp,type,mass);

	    type = basistype[14];
	    if(type==0) mass=m; else mass=m1;
	    tp.x = p.x + b14_rot2.x;
	    tp.y = p.y + b14_rot2.y;
	    tp.z = p.z + b14_rot2.z;
	    write_atom(tp,type,mass);

	    type = basistype[15];
	    if(type==0) mass=m; else mass=m1;
	    tp.x = p.x + b15_rot2.x;
	    tp.y = p.y + b15_rot2.y;
	    tp.z = p.z + b15_rot2.z;
	    write_atom(tp,type,mass);

	    type = basistype[16];
	    if(type==0) mass=m; else mass=m1;
	    tp.x = p.x + b16_rot2.x;
	    tp.y = p.y + b16_rot2.y;
	    tp.z = p.z + b16_rot2.z;
	    write_atom(tp,type,mass);
	  }
	

	
      }
  printf("Min. Super Cell dim: %.16f %.16f %.16f\n",scell_min.x, scell_min.y, scell_min.z);
  printf("Max. Super Cell dim: %.16f %.16f %.16f\n",scell_max.x, scell_max.y, scell_max.z);
  printf("Next atom nr.: %d\n",nr);
  fclose(out);
  return(0);

}



void write_atom(vektor p, int typ,double masse)
{
  double pn1,pn2;
  /* print if... */
  /* ...atom is in box */
  if( p.x+e >= 0.0 && p.x+e < box_x &&         /* rounding errors... */
      p.y+e >= 0.0 && p.y+e < box_y &&
      p.z+e >= 0.0 && p.z+e < box_z ){
    
    pn1 = smult(p,n1);
    pn2 = smult(p,n2);
    
    /* 	  printf(" %f %f %f \n",p.x, p.y,p.z); */
    /* 	  printf(" %f \n",pn1); */
    /* 	  printf(" %f \n",pn2); */
    /* 	  fflush(stdout); */
    
    /* and atom is not cut away */
    if(cutdir1 < 0.0) /* cut away if on the same side as 0 */
      if(pn1 > d1){ 
	if(cutdir2 < 0.0)
	  if(pn2 > d2){ 
	    fprintf(out,"%d %d %.16f %.16f %.16f %.16f \n",
		    nr++, typ, masse, p.x, p.y, p.z);
	  }
	if(cutdir2 > 0.0)
	  //		if(pn2 <= d2){
	  if(pn2 < d2){
	    fprintf(out,"%d %d %.16f %.16f %.16f %.16f \n",
		    nr++, typ, masse, p.x, p.y, p.z);
	  }
	if(cutdir2 == 0.0){
	  fprintf(out,"%d %d %.16f %.16f %.16f %.16f \n",
		  nr++, typ, masse, p.x, p.y, p.z);
	}
      }
    if(cutdir1 > 0.0) /* cut away if on the other side as 0 */
      //	    if(pn1 <= d1){
      if(pn1 < d1){
	if(cutdir2 < 0.0)
	  if(pn2 > d2){
	    fprintf(out,"%d %d %.16f %.16f %.16f %.16f \n",
		    nr++, typ, masse, p.x, p.y, p.z);
	  }
	if(cutdir2 > 0.0)
	  //	if(pn2 <= d2){
	  if(pn2 < d2){
	    fprintf(out,"%d %d %.16f %.16f %.16f %.16f \n",
		    nr++, typ, masse, p.x, p.y, p.z);
	  }
	if(cutdir2 == 0.0){
	  fprintf(out,"%d %d %.16f %.16f %.16f %.16f \n",
		  nr++, typ, masse, p.x, p.y, p.z);
	}
      }
    if(cutdir1 == 0.0){ /* don't cut away anything */
      if(cutdir2 < 0.0)
	if(pn2 > d2){ 
	  fprintf(out,"%d %d %.16f %.16f %.16f %.16f \n",
			nr++, typ, masse, p.x, p.y, p.z);
	}
      if(cutdir2 > 0.0)
	//      if(pn2 <= d2){
                if(pn2 < d2){ 
		  fprintf(out,"%d %d %.16f %.16f %.16f %.16f \n",
			  nr++, typ, masse, p.x, p.y, p.z);
		}
      if(cutdir2 == 0.0){
	fprintf(out,"%d %d %.16f %.16f %.16f %.16f \n",
		nr++, typ, masse, p.x, p.y, p.z);
	/*  printf(" %d %d %d\n",i,j,k);  */
	
	/* find (smallest) cellsize for periodic boundaries* 
	 * only meaningful when nothing is cut away        */
	if (!( (p.x ==0.0) && (p.y ==0.0) && (p.z == 0.0) )){
	  p.x = fabs(p.x);
	  p.y = fabs(p.y);
	  p.z = fabs(p.z);
	  
	  if((p.y < EPS) && (p.z < EPS)){
	    /*   printf("cellvektor x: %.16f %.16f %.16f \n", p.x, p.y, p.z); */
	    if (p.x < scell_min.x)
	      scell_min.x = p.x;
	    if (p.x > scell_max.x)
	      scell_max.x = p.x;
	  }
	  if((p.x < EPS) && (p.z < EPS)){
	    /* 		  printf("cellvektor y: %.16f %.16f %.16f \n", p.x, p.y, p.z);  */
	    if (p.y < scell_min.y)
	      scell_min.y = p.y;
	    if (p.y > scell_max.y)
	      scell_max.y = p.y;
	  }
	  if((p.x < EPS) && (p.y < EPS)){
	    /* 		  printf("cellvektor z: %.16f %.16f %.16f \n", p.x, p.y, p.z);  */
	    if (p.z < scell_min.z)
	      scell_min.z = p.z;
	    if (p.z > scell_max.z)
	      scell_max.z = p.z;
	  }
	}
      }
    }
  }
}
