/* LEGO (Lattice Elemental Geometry Operations)
 * Program to create arbitary rotated and cut FCC structures in IMD format 
 * Erik Bitzek Jan 2000
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

int main(int argc, char *argv[])
{
  int i,j,k;
  vektor p;
  int nr;
  double pn1,pn2,n1abs,n2abs;
  double nxabs,nyabs,nzabs;

  vektor scell_min;
  vektor scell_max;    /* super cell dimension */
  

  double to_ps=0.000103650;     /* amu -> eV*psec^2/angstrom^2 */
  
  /* (rotated) translation vectors */
  vektor a1,a2,a3;
  vektor a1_rot INIT(nullvektor);
  vektor a2_rot INIT(nullvektor);
  vektor a3_rot INIT(nullvektor);
  vektor a1_rot2 INIT(nullvektor);
  vektor a2_rot2 INIT(nullvektor);
  vektor a3_rot2 INIT(nullvektor);

  
   /* (rotated) standard vectors */
  vektor e1,e2,e3;
  vektor e1_rot INIT(nullvektor);
  vektor e2_rot INIT(nullvektor);
  vektor e3_rot INIT(nullvektor);
  vektor e1_rot2 INIT(nullvektor);
  vektor e2_rot2 INIT(nullvektor);
  vektor e3_rot2 INIT(nullvektor);
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
  if (strcasecmp(struktur,"fcc")==0)
    {
      /* primitve translations of fcc */
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
  else 
    {
      printf("Unknown Crystalstructure : %s\n",struktur);
      exit(2);
    }


  if(det(drehmatrix1)==0.0)
    {
       
        /* creation of transformation matrix */
        nxabs = vabs(nx);
        nyabs = vabs(ny);
        nzabs = vabs(nz);
        
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
            nz.y *= 1.0/nzabs;
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
  
/*   printf("1. Rotation:\n"); */
/*   printf(" %f %f %f\n\n",a1_rot.x,a1_rot.y,a1_rot.z); */
/*   printf(" %f %f %f\n\n",a2_rot.x,a2_rot.y,a2_rot.z); */
/*   printf(" %f %f %f\n\n",a3_rot.x,a3_rot.y,a3_rot.z); */

  mvmult(drehmatrix2,a1_rot,&a1_rot2);
  mvmult(drehmatrix2,a2_rot,&a2_rot2);
  mvmult(drehmatrix2,a3_rot,&a3_rot2);

  printf("unit vectors after 2. Rotation:\n");
  printf(" %.8f %.8f %.8f\n",a1_rot2.x,a1_rot2.y,a1_rot2.z);
  printf(" %.8f %.8f %.8f\n",a2_rot2.x,a2_rot2.y,a2_rot2.z);
  printf(" %.8f %.8f %.8f\n\n",a3_rot2.x,a3_rot2.y,a3_rot2.z);
  
  mvmult(drehmatrix1,e1,&e1_rot);
  mvmult(drehmatrix1,e2,&e2_rot);
  mvmult(drehmatrix1,e3,&e3_rot);

  mvmult(drehmatrix2,e1_rot,&e1_rot2);
  mvmult(drehmatrix2,e2_rot,&e2_rot2);
  mvmult(drehmatrix2,e3_rot,&e3_rot2);

  printf("Rotation matrix :\n");
  /* printf(" %.12f %.12f %.12f\n",e1_rot2.x,e1_rot2.y,e1_rot2.z); */
  /* printf(" %.12f %.12f %.12f\n",e2_rot2.x,e2_rot2.y,e2_rot2.z); */
  /* printf(" %.12f %.12f %.12f\n\n",e3_rot2.x,e3_rot2.y,e3_rot2.z); */
  // this matrix is used by GBcharacterizer and in matlab
  // this is also how it will be read in 
  printf(" %.12f %.12f %.12f\n",e1_rot2.x,e2_rot2.x,e3_rot2.x);
  printf(" %.12f %.12f %.12f\n",e1_rot2.y,e2_rot2.y,e3_rot2.y);
  printf(" %.12f %.12f %.12f\n\n",e1_rot2.z,e2_rot2.z,e3_rot2.z);
  
  
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
		  fprintf(out,"%d 0 %.16f %.16f %.16f %.16f \n",
			  nr++, m, p.x, p.y, p.z);
		}
	      if(cutdir2 > 0.0)
		if(pn2 < d2){
		  fprintf(out,"%d 0 %.16f %.16f %.16f %.16f \n",
			  nr++, m, p.x, p.y, p.z);
		}
	      if(cutdir2 == 0.0){
		fprintf(out,"%d 0 %.16f %.16f %.16f %.16f \n",
			nr++, m, p.x, p.y, p.z);
	      }
	    }
	  if(cutdir1 > 0.0) /* cut away if on the other side as 0 */
	    if(pn1 < d1){
	      if(cutdir2 < 0.0)
		if(pn2 > d2){
		  fprintf(out,"%d 0 %.16f %.16f %.16f %.16f \n",
			  nr++, m, p.x, p.y, p.z);
		}
	      if(cutdir2 > 0.0)
		if(pn2 < d2){
		  fprintf(out,"%d 0 %.16f %.16f %.16f %.16f \n",
			  nr++, m, p.x, p.y, p.z);
		}
	      if(cutdir2 == 0.0){
		fprintf(out,"%d 0 %.16f %.16f %.16f %.16f \n",
			nr++, m, p.x, p.y, p.z);
	      }
	    }
	  if(cutdir1 == 0.0){ /* don't cut away anything */
	    if(cutdir2 < 0.0)
	      if(pn2 > d2){ 
		fprintf(out,"%d 0 %.16f %.16f %.16f %.16f \n",
			nr++, m, p.x, p.y, p.z);
	      }
	    if(cutdir2 > 0.0)
	      if(pn2 < d2){ 
		fprintf(out,"%d 0 %.16f %.16f %.16f %.16f \n",
			nr++, m, p.x, p.y, p.z);
	      }
	    if(cutdir2 == 0.0){
	      fprintf(out,"%d 0 %.16f %.16f %.16f %.16f \n",
		      nr++, m, p.x, p.y, p.z);
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
  printf("Min. Super Cell dim: %.16f %.16f %.16f\n",scell_min.x, scell_min.y, scell_min.z);
  printf("Max. Super Cell dim: %.16f %.16f %.16f\n",scell_max.x, scell_max.y, scell_max.z);
  printf("Next atom nr.: %d\n",nr);
  fclose(out);
  return(0);

}
