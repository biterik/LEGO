/* Paramerter File handling for mk_config.c */

/* stolen from imd_param.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mk_config_globals.h" 


/******************************************************************************
*
* error -- Complain and abort
*
******************************************************************************/

void error(char *mesg)
{
  fprintf(stderr,"Error: %s\n",mesg);
  exit(2);
}


typedef enum ParamType {
  PARAM_STR, PARAM_STRPTR,
  PARAM_INT, PARAM_INT_COPY,
  PARAM_INTEGER, PARAM_INTEGER_COPY,
  PARAM_REAL, PARAM_REAL_COPY
} PARAMTYPE;

int curline; /* number of current line */

/*****************************************************************************
*
* Parameter aus Zeile auslesen / get parameter from line
*
*****************************************************************************/

/* Parameter:
   param_name ... Parametername (fuer Fehlermeldungen)
   param ........ Adresse der Variable fuer den Parameter
   ptype ........ Parametertyp
                  folgende Werte sind zulaessig:
                  PARAM_STR : String, deklariert als char[]
                  PARAM_STRPTR : String, deklariert als Zeiger auf char*
                  PARAM_INT : Integer-Wert(e)
                  PARAM_INT_COPY : Integer-Wert(e), kopierend
                  PARAM_REAL : Real-Wert(e)
                  PARAM_REAL_COPY : Real-Wert(e), kopierend
                  
   pnum_min ..... Minimale Anzahl der einzulesenden Werte
                  (Falls weniger Werte gelesen werden koennen als verlangt,
                  wird ein Fehler gemeldet).
   pnum_max ..... Maximale Anzahl der einzulesenden Werte
                  (Die nicht kopierenden Routinen lesen hoechstens
                  pnum_max Werte aus der uebergebenen Zeile ein,
                  weitere Werte werden ignoriert. Falls weniger als
                  pnum_max Werte vorhanden sind, wird das Lesen
                  abgebrochen, es wird kein Fehler gemeldet,
                  wenn mindestens pnum_min Werte abgesaettigt wurden.
                  Die kopierenden Routinen melden ebenfalls keinen
                  Fehler, wenn mindestens pnum_min Werte abgesaettigt
                  wurden. Falls weniger als pnum_max Werte vorhanden sind,
                  werden die restlichen Werte mit Kopien des zuletzt gelesenen
                  Werts aufgefuellt.

  Resultat:
  nichtkopierende Routinen: Die Anzahl der gelesenen Werte wird zurueckgegeben.
  kopierende Routinen: Die Anzahl der tatsaechlich gelesenen Werte wird
                       zurueckgegeben. Resultat = pnum_max - Anzahl der Kopien
*/

int getparam(char *param_name, void *param, PARAMTYPE ptype, 
             int pnum_min, int pnum_max)
{
  static char errmsg[256];
  char *str;
  int i;
  int numread;

  numread = 0;
  if (ptype == PARAM_STR) {
    str = strtok(NULL," \t\n");
    if (str == NULL) {
      fprintf(stderr,"parameter for %s missing in line %u\n",
              param_name,curline);
      error("string expected\n");
    }
    else strncpy((char *)param,str,pnum_max);
    numread++;
  }
  else if (ptype == PARAM_STRPTR) {
    str = strtok(NULL," \t\n");
    if (str == NULL) {
      fprintf(stderr,"parameter for %s missing in line %u\n",
              param_name,curline);
      error("string expected\n");
    }
    else *((char**)param) = strdup(str);
    numread++;
  }
  else if (ptype == PARAM_INT) {
    for (i=0; i<pnum_min; i++) {
      str = strtok(NULL," \t\n");
      if (str == NULL) {
        fprintf(stderr,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg,"integer vector of dim %u expected\n",
                (unsigned)pnum_min);
        error(errmsg);
      }
      else ((int*)param)[i] = atoi(str);
      numread++;
    }
    for (i=pnum_min; i<pnum_max; i++) {
      if ((str = strtok(NULL," \t\n")) != NULL) {
        ((int*)param)[i] = atoi(str);
        numread++;
      }
      else break;
    }
  }
  else if (ptype == PARAM_INT_COPY) {
    int ival = 0;
    for (i=0; i<pnum_max; i++) {
      str = strtok(NULL," \t\n");
      if (str != NULL) {
        ival = atoi(str);
        numread++; /* return number of parameters actually read */
      }
      else if (i<pnum_min) { 
        fprintf(stderr,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg,"integer vector of dim %u expected\n",
                (unsigned)pnum_min);
        error(errmsg);
      };
      ((int*)param)[i] = ival;
    }
  }
  else if (ptype == PARAM_INTEGER) {
    for (i=0; i<pnum_min; i++) {
      str = strtok(NULL," \t\n");
      if (str == NULL) {
        fprintf(stderr,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg,"integer vector of dim %u expected\n",
                (unsigned)pnum_min);
        error(errmsg);
      }
      else ((int*)param)[i] = atoi(str);
      numread++;
    }
    for (i=pnum_min; i<pnum_max; i++) {
      if ((str = strtok(NULL," \t\n")) != NULL) {
        ((int*)param)[i] = atoi(str);
        numread++;
      }
      else break;
    }
  }
  else if (ptype == PARAM_INTEGER_COPY) {
    int ival = 0;
    for (i=0; i<pnum_max; i++) {
      str = strtok(NULL," \t\n");
      if (str != NULL) {
        ival = atoi(str);
        numread++; /* return number of parameters actually read */
      }
      else if (i<pnum_min) {
        fprintf(stderr,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg,"integer vector of dim %u expected\n",
                (unsigned)pnum_min);
        error(errmsg);
      };
      ((int*)param)[i] = (int)ival;
    }
  }
  else if (ptype == PARAM_REAL) {
    for (i=0; i<pnum_min; i++) {
      str = strtok(NULL," \t\n");
      if (str == NULL) {
        fprintf(stderr,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg,"real vector of dim %u expected\n",
                (unsigned)pnum_min);
        error(errmsg);
      }
      else ((real*)param)[i] = atof(str);
      numread++;
    }
    for (i=pnum_min; i<pnum_max; i++) {
      if ((str = strtok(NULL," \t\n")) != NULL) {
        ((real*)param)[i] = atof(str);
        numread++;
      }
      else break;
    }
  }
  else if (ptype == PARAM_REAL_COPY) {
    real rval = 0;
    for (i=0; i<pnum_max; i++) {
      str = strtok(NULL," \t\n");
      if (str != NULL) {
        rval = atof(str);
        numread++; /* return number of parameters actually read */
      }
      else if (i<pnum_min) {
        fprintf(stderr,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg,"real vector of dim %u expected\n",
                (unsigned)pnum_min);
        error(errmsg);
      };
      ((real*)param)[i] = rval;
    }
  }
  return numread;
} /* getparam */






void getparams(char *inputfilename)
{
  FILE *infile;
  char buffer[1024];
  char *token;
  char *res;
  curline = 0;
  
  /* open file */
  infile = fopen(inputfilename,"r");
  if (NULL==infile){
    printf("Can't open inputfile %s\n",inputfilename);
    exit(2);
  }

  /* get parameter */
  do {
    res=fgets(buffer,1024,infile);
    if (NULL == res) { break; }; /* probably EOF reached */
    curline++;
    token = strtok(buffer," \t\n");
    if (NULL == token) continue; /* skip blank lines */
    if (token[0]=='#') continue; /* skip comments */

    if (strcasecmp(token,"outfile")==0) {
      /* output filename */
      getparam("outfile",outfile,PARAM_STR,1,255);
    }
    else if (strcasecmp(token,"structure")==0) {
      /* output filename */
      getparam("structure",struktur,PARAM_STR,1,255);
    }
     else if (strcasecmp(token,"lower")==0) {
      /* start value for linear superpos. of vectors */
      getparam("lower",&lower,PARAM_INT,1,1);
    }
     else if (strcasecmp(token,"upper")==0) {
      /* end value for linear superpos. of vectors */
      getparam("upper",&upper,PARAM_INT,1,1);
    }
     else if (strcasecmp(token,"epsilon")==0) {
      /* epsilon to try to avoid rounding errors */
      getparam("epsilon",&e,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"startnr")==0) {
      /* start value for atom numbering */
      getparam("startnr",&start,PARAM_INT,1,1);
    }
     else if (strcasecmp(token,"atomtype")==0) {
      /* start value for atom numbering */
      getparam("atomtype",&type,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"mass")==0) {
      /* Atom mass in u */
      getparam("mass",&m,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"mass1")==0) {
      /* Atom mass in u */
      getparam("mass1",&m1,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"lattice_const")==0) {
      /* Lattice Constant in Angstrom */
      getparam("lattice_const",&a,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"lattice_const_b")==0) {
      /* 2nd Lattice Constant in Angstrom */
      getparam("lattice_consti_b",&b,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"lattice_const_c")==0) {
      /* 3rd Lattice Constant in Angstrom */
      getparam("lattice_const_c",&c,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"shiftvec")==0) {
      /* shiftvector for construction of multibase struct */
      getparam("shiftvec",&shiftvec,PARAM_REAL,3,3);
    }
   else if (strcasecmp(token,"pos_8d_init")==0) {
      /* initial values of general Wyckoff 8d sites */
      getparam("pos_8d_init",&i8d,PARAM_REAL,3,3);
    }
    else if (strcasecmp(token,"pos_4c_init")==0) {
      /* initial values of special Wyckoff 4c sites  */
      getparam("pos_4c_init",&i4c,PARAM_REAL,3,3);
    }
    else if (strcasecmp(token,"pos_4c_init1")==0) {
      /* initial values of special Wyckoff 4c sites for atom type 1 */
      getparam("pos_4c_init",&i4c1,PARAM_REAL,3,3);
    }	
    else if (strcasecmp(token,"box_x")==0) {
      /* box dimension 'x' in Angstrom */
      getparam("box_x",&box_x,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"box_y")==0) {
      /* box dimension 'y' */
      getparam("box_y",&box_y,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"box_z")==0) {
      /* box dimension 'z' */
      getparam("box_z",&box_z,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"cut_norm1")==0) {
      /*  normal of cut plane */
      getparam("cut_norm1",&n1,PARAM_REAL,3,3);
    }
    else if (strcasecmp(token,"plane_distance1")==0) {
      /* cut plane distance to 0 */
      getparam("plane_distance1",&d1,PARAM_REAL,1,1);
    }
     else if (strcasecmp(token,"cutdir1")==0) {
      /* which side to cut away */
      getparam("cutdir1",&cutdir1,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"cut_norm2")==0) {
      /*  normal of cut plane */
      getparam("cut_norm2",&n2,PARAM_REAL,3,3);
    }
    else if (strcasecmp(token,"plane_distance2")==0) {
      /* cut plane distance to 0 */
      getparam("plane_distance2",&d2,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"plane_point2")==0) {
      /* cut plane distance to 0 */
      getparam("plane_point2",&point2,PARAM_REAL,3,3);
    }
    else if (strcasecmp(token,"plane_point1")==0) {
      /* cut plane distance to 0 */
      getparam("plane_point1",&point1,PARAM_REAL,3,3);
    }
     else if (strcasecmp(token,"cutdir2")==0) {
      /* which side to cut away */
      getparam("cutdir2",&cutdir2,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"new_x")==0) {
      /* new axis should correspond to */
      getparam("new_x",&nx,PARAM_REAL,3,3);
    } 
    else if (strcasecmp(token,"new_y")==0) {
      /* new axis should correspond to */
      getparam("new_y",&ny,PARAM_REAL,3,3);
    } 
    else if (strcasecmp(token,"new_z")==0) {
      /* new axis should correspond to */
      getparam("new_z",&nz,PARAM_REAL,3,3);
    } 
    else if (strcasecmp(token,"axis1")==0) {
      /* vector of rotation axis1 */
      getparam("axis1",&axe1,PARAM_REAL,3,3);
    } 
    else if (strcasecmp(token,"axis2")==0) {
      /* vector of rotation axis2 */
      getparam("axis2",&axe2,PARAM_REAL,3,3);
    }
     else if (strcasecmp(token,"rotmat1_x")==0) {
      /* vector of rotation axis2 */
      getparam("rotmat1_x",&m1_x,PARAM_REAL,3,3);
      drehmatrix1[0][0]=m1_x.x;
      drehmatrix1[1][0]=m1_x.y;
      drehmatrix1[2][0]=m1_x.z;
     }
      else if (strcasecmp(token,"rotmat1_y")==0) {
      /* vector of rotation axis2 */
      getparam("rotmat1_2",&m1_y,PARAM_REAL,3,3);
      drehmatrix1[0][1]=m1_y.x;
      drehmatrix1[1][1]=m1_y.y;
      drehmatrix1[2][1]=m1_y.z;
    }
       else if (strcasecmp(token,"rotmat1_z")==0) {
      /* vector of rotation axis2 */
      getparam("rotmat1_3",&m1_z,PARAM_REAL,3,3);
      drehmatrix1[0][2]=m1_z.x;
      drehmatrix1[1][2]=m1_z.y;
      drehmatrix1[2][2]=m1_z.z;
    }
    else if (strcasecmp(token,"angle1")==0) {
      /* rotation angle 1*/
      getparam("angle1",&winkel1,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"angle2")==0) {
      /* rotation angle 1*/
      getparam("angle2",&winkel2,PARAM_REAL,1,1);
    }
    else {
      fprintf(stderr,"**** WARNING: Unknown TAG %s ignored ****\n",token);
    }
  } while (!feof(infile));
  fclose(infile);

} /* getparamfile */

