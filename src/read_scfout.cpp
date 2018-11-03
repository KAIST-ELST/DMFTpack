/**********************************************************************
  read_scfout.c:

     read_scfout.c is a subroutine to read a binary file,
     filename.scfout.

  Log of read_scfout.c:

     2/July/2003  Released by T.Ozaki

***********************************************************************/
//(modified by Jae-Hoon Sim, 20180816)

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <string.h>
#include <fstream>
//#include <cstring>
#include "read_scfout.h"

#define MAX_LINE_SIZE 256
#define fp_bsize         1048576     /* buffer size for setvbuf */

static void Input( FILE *fp );

int read_scfout(std::string scfout_file)
{
    FILE *fp;
    char buf[fp_bsize];          /* setvbuf */

    if ((fp = fopen(scfout_file.c_str(),"r")) != NULL) {

#ifdef xt3
        setvbuf(fp,   buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

        printf("\nRead the scfout file (%s)\n",scfout_file.c_str());
        fflush(stdout);
        Input(fp);
        fclose(fp);
        return 1;
    }
    else {
        printf("Failure of reading the scfout file (%s).\n",scfout_file.c_str());
        fflush(stdout);
        return 0;
    }

}


void Input( FILE *fp )
{
    static int Gc_AN,ct_AN,h_AN,i,j,can,Gh_AN;
    static int wan1,wan2,TNO1,TNO2,spin,Rn,num_lines;
    static int k,q_AN,Gq_AN;
    static int i_vec[20],*p_vec;
    static double d_vec[20];
    static char makeinp[100];
    static char strg[MAX_LINE_SIZE];
    FILE *fp_makeinp;
    char buf[fp_bsize];          /* setvbuf */

    /****************************************************
       atomnum
       spinP_switch
    ****************************************************/

    fread(i_vec,sizeof(int),6,fp);
    atomnum      = i_vec[0];
    SpinP_switch = i_vec[1];
    Catomnum =     i_vec[2];
    Latomnum =     i_vec[3];
    Ratomnum =     i_vec[4];
    TCpyCell =     i_vec[5];

    /****************************************************
      allocation of arrays:

      double atv[TCpyCell+1][4];
    ****************************************************/

    atv = (double**)malloc(sizeof(double*)*(TCpyCell+1));
    for (Rn=0; Rn<=TCpyCell; Rn++) {
        atv[Rn] = (double*)malloc(sizeof(double)*4);
    }

    /****************************************************
                  read atv[TCpyCell+1][4];
    ****************************************************/

    for (Rn=0; Rn<=TCpyCell; Rn++) {
        fread(atv[Rn],sizeof(double),4,fp);
    }

    /****************************************************
      allocation of arrays:

      int atv_ijk[TCpyCell+1][4];
    ****************************************************/

    atv_ijk = (int**)malloc(sizeof(int*)*(TCpyCell+1));
    for (Rn=0; Rn<=TCpyCell; Rn++) {
        atv_ijk[Rn] = (int*)malloc(sizeof(int)*4);
    }

    /****************************************************
              read atv_ijk[TCpyCell+1][4];
    ****************************************************/

    for (Rn=0; Rn<=TCpyCell; Rn++) {
        fread(atv_ijk[Rn],sizeof(int),4,fp);
    }

    /****************************************************
      allocation of arrays:

      int Total_NumOrbs[atomnum+1];
      int FNAN[atomnum+1];
    ****************************************************/

    Total_NumOrbs = (int*)malloc(sizeof(int)*(atomnum+1));
    FNAN = (int*)malloc(sizeof(int)*(atomnum+1));

    /****************************************************
           the number of orbitals in each atom
    ****************************************************/

    p_vec = (int*)malloc(sizeof(int)*atomnum);
    fread(p_vec,sizeof(int),atomnum,fp);
    Total_NumOrbs[0] = 1;
    for (ct_AN=1; ct_AN<=atomnum; ct_AN++) {
        Total_NumOrbs[ct_AN] = p_vec[ct_AN-1];
    }
    free(p_vec);

    /****************************************************
     FNAN[]:
     the number of first nearest neighbouring atoms
    ****************************************************/

    p_vec = (int*)malloc(sizeof(int)*atomnum);
    fread(p_vec,sizeof(int),atomnum,fp);
    FNAN[0] = 0;
    for (ct_AN=1; ct_AN<=atomnum; ct_AN++) {
        FNAN[ct_AN] = p_vec[ct_AN-1];
    }
    free(p_vec);

    /****************************************************
      allocation of arrays:

      int natn[atomnum+1][FNAN[ct_AN]+1];
      int ncn[atomnum+1][FNAN[ct_AN]+1];
    ****************************************************/

    natn = (int**)malloc(sizeof(int*)*(atomnum+1));
    for (ct_AN=0; ct_AN<=atomnum; ct_AN++) {
        natn[ct_AN] = (int*)malloc(sizeof(int)*(FNAN[ct_AN]+1));
    }

    ncn = (int**)malloc(sizeof(int*)*(atomnum+1));
    for (ct_AN=0; ct_AN<=atomnum; ct_AN++) {
        ncn[ct_AN] = (int*)malloc(sizeof(int)*(FNAN[ct_AN]+1));
    }

    /****************************************************
      natn[][]:
      grobal index of neighboring atoms of an atom ct_AN
     ****************************************************/

    for (ct_AN=1; ct_AN<=atomnum; ct_AN++) {
        fread(natn[ct_AN],sizeof(int),FNAN[ct_AN]+1,fp);
    }

    /****************************************************
      ncn[][]:
      grobal index for cell of neighboring atoms
      of an atom ct_AN
    ****************************************************/

    for (ct_AN=1; ct_AN<=atomnum; ct_AN++) {
        fread(ncn[ct_AN],sizeof(int),FNAN[ct_AN]+1,fp);
    }

    /****************************************************
      tv[4][4]:
      unit cell vectors in Bohr
    ****************************************************/

    fread(tv[1],sizeof(double),4,fp);
    fread(tv[2],sizeof(double),4,fp);
    fread(tv[3],sizeof(double),4,fp);

    /****************************************************
      rtv[4][4]:
      unit cell vectors in Bohr
    ****************************************************/

    fread(rtv[1],sizeof(double),4,fp);
    fread(rtv[2],sizeof(double),4,fp);
    fread(rtv[3],sizeof(double),4,fp);

    /****************************************************
      Gxyz[][1-3]:
      atomic coordinates in Bohr
    ****************************************************/

    Gxyz = (double**)malloc(sizeof(double*)*(atomnum+1));
    for (i=0; i<(atomnum+1); i++) {
        Gxyz[i] = (double*)malloc(sizeof(double)*60);
    }

    for (ct_AN=1; ct_AN<=atomnum; ct_AN++) {
        fread(Gxyz[ct_AN],sizeof(double),4,fp);
    }

    /****************************************************
      allocation of arrays:

      Kohn-Sham Hamiltonian

       dooble Hks[SpinP_switch+1]
                 [atomnum+1]
                 [FNAN[ct_AN]+1]
                 [Total_NumOrbs[ct_AN]]
                 [Total_NumOrbs[h_AN]];

      Overlap matrix

       dooble OLP[atomnum+1]
                 [FNAN[ct_AN]+1]
                 [Total_NumOrbs[ct_AN]]
                 [Total_NumOrbs[h_AN]];

      Overlap matrix with position operator x, y, z

       dooble OLPpox,y,z
                   [atomnum+1]
                   [FNAN[ct_AN]+1]
                   [Total_NumOrbs[ct_AN]]
                   [Total_NumOrbs[h_AN]];

      Density matrix

       dooble DM[SpinP_switch+1]
                [atomnum+1]
                [FNAN[ct_AN]+1]
                [Total_NumOrbs[ct_AN]]
                [Total_NumOrbs[h_AN]];
    ****************************************************/

    Hks = (double*****)malloc(sizeof(double****)*(SpinP_switch+1));
    for (spin=0; spin<=SpinP_switch; spin++) {

        Hks[spin] = (double****)malloc(sizeof(double***)*(atomnum+1));
        for (ct_AN=0; ct_AN<=atomnum; ct_AN++) {
            TNO1 = Total_NumOrbs[ct_AN];
            Hks[spin][ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN]+1));
            for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++) {
                Hks[spin][ct_AN][h_AN] = (double**)malloc(sizeof(double*)*TNO1);

                if (ct_AN==0) {
                    TNO2 = 1;
                }
                else {
                    Gh_AN = natn[ct_AN][h_AN];
                    TNO2 = Total_NumOrbs[Gh_AN];
                }
                for (i=0; i<TNO1; i++) {
                    Hks[spin][ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*TNO2);
                }
            }
        }
    }

    iHks = (double*****)malloc(sizeof(double****)*3);
    for (spin=0; spin<3; spin++) {

        iHks[spin] = (double****)malloc(sizeof(double***)*(atomnum+1));
        for (ct_AN=0; ct_AN<=atomnum; ct_AN++) {
            TNO1 = Total_NumOrbs[ct_AN];
            iHks[spin][ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN]+1));
            for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++) {
                iHks[spin][ct_AN][h_AN] = (double**)malloc(sizeof(double*)*TNO1);

                if (ct_AN==0) {
                    TNO2 = 1;
                }
                else {
                    Gh_AN = natn[ct_AN][h_AN];
                    TNO2 = Total_NumOrbs[Gh_AN];
                }
                for (i=0; i<TNO1; i++) {
                    iHks[spin][ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*TNO2);
                    for (j=0; j<TNO2; j++) iHks[spin][ct_AN][h_AN][i][j] = 0.0;
                }
            }
        }
    }

    OLP = (double****)malloc(sizeof(double***)*(atomnum+1));
    for (ct_AN=0; ct_AN<=atomnum; ct_AN++) {
        TNO1 = Total_NumOrbs[ct_AN];
        OLP[ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN]+1));
        for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++) {
            OLP[ct_AN][h_AN] = (double**)malloc(sizeof(double*)*TNO1);

            if (ct_AN==0) {
                TNO2 = 1;
            }
            else {
                Gh_AN = natn[ct_AN][h_AN];
                TNO2 = Total_NumOrbs[Gh_AN];
            }
            for (i=0; i<TNO1; i++) {
                OLP[ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*TNO2);
            }
        }
    }


    OLPpox = (double****)malloc(sizeof(double***)*(atomnum+1));
    for (ct_AN=0; ct_AN<=atomnum; ct_AN++) {
        TNO1 = Total_NumOrbs[ct_AN];
        OLPpox[ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN]+1));
        for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++) {
            OLPpox[ct_AN][h_AN] = (double**)malloc(sizeof(double*)*TNO1);

            if (ct_AN==0) {
                TNO2 = 1;
            }
            else {
                Gh_AN = natn[ct_AN][h_AN];
                TNO2 = Total_NumOrbs[Gh_AN];
            }
            for (i=0; i<TNO1; i++) {
                OLPpox[ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*TNO2);
            }
        }
    }

    OLPpoy = (double****)malloc(sizeof(double***)*(atomnum+1));
    for (ct_AN=0; ct_AN<=atomnum; ct_AN++) {
        TNO1 = Total_NumOrbs[ct_AN];
        OLPpoy[ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN]+1));
        for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++) {
            OLPpoy[ct_AN][h_AN] = (double**)malloc(sizeof(double*)*TNO1);

            if (ct_AN==0) {
                TNO2 = 1;
            }
            else {
                Gh_AN = natn[ct_AN][h_AN];
                TNO2 = Total_NumOrbs[Gh_AN];
            }
            for (i=0; i<TNO1; i++) {
                OLPpoy[ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*TNO2);
            }
        }
    }

    OLPpoz = (double****)malloc(sizeof(double***)*(atomnum+1));
    for (ct_AN=0; ct_AN<=atomnum; ct_AN++) {
        TNO1 = Total_NumOrbs[ct_AN];
        OLPpoz[ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN]+1));
        for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++) {
            OLPpoz[ct_AN][h_AN] = (double**)malloc(sizeof(double*)*TNO1);

            if (ct_AN==0) {
                TNO2 = 1;
            }
            else {
                Gh_AN = natn[ct_AN][h_AN];
                TNO2 = Total_NumOrbs[Gh_AN];
            }
            for (i=0; i<TNO1; i++) {
                OLPpoz[ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*TNO2);
            }
        }
    }

    DM = (double*****)malloc(sizeof(double****)*(SpinP_switch+1));
    for (spin=0; spin<=SpinP_switch; spin++) {

        DM[spin] = (double****)malloc(sizeof(double***)*(atomnum+1));
        for (ct_AN=0; ct_AN<=atomnum; ct_AN++) {
            TNO1 = Total_NumOrbs[ct_AN];
            DM[spin][ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN]+1));
            for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++) {
                DM[spin][ct_AN][h_AN] = (double**)malloc(sizeof(double*)*TNO1);

                if (ct_AN==0) {
                    TNO2 = 1;
                }
                else {
                    Gh_AN = natn[ct_AN][h_AN];
                    TNO2 = Total_NumOrbs[Gh_AN];
                }
                for (i=0; i<TNO1; i++) {
                    DM[spin][ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*TNO2);
                }
            }
        }
    }

    /****************************************************
                      Hamiltonian matrix
    ****************************************************/

    for (spin=0; spin<=SpinP_switch; spin++) {
        for (ct_AN=1; ct_AN<=atomnum; ct_AN++) {
            TNO1 = Total_NumOrbs[ct_AN];
            for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++) {
                Gh_AN = natn[ct_AN][h_AN];
                TNO2 = Total_NumOrbs[Gh_AN];
                for (i=0; i<TNO1; i++) {
                    fread(Hks[spin][ct_AN][h_AN][i],sizeof(double),TNO2,fp);
                }
            }
        }
    }

    /****************************************************
    iHks:
    imaginary Kohn-Sham matrix elements of basis orbitals
    for alpha-alpha, beta-beta, and alpha-beta spin matrices
    of which contributions come from spin-orbit coupling
    and Hubbard U effective potential.
    ****************************************************/

    if (SpinP_switch==3) {
        for (spin=0; spin<3; spin++) {
            for (ct_AN=1; ct_AN<=atomnum; ct_AN++) {
                TNO1 = Total_NumOrbs[ct_AN];
                for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++) {
                    Gh_AN = natn[ct_AN][h_AN];
                    TNO2 = Total_NumOrbs[Gh_AN];
                    for (i=0; i<TNO1; i++) {
                        fread(iHks[spin][ct_AN][h_AN][i],sizeof(double),TNO2,fp);
                    }
                }
            }
        }
    }

    /****************************************************
                       Overlap matrix
    ****************************************************/

    for (ct_AN=1; ct_AN<=atomnum; ct_AN++) {
        TNO1 = Total_NumOrbs[ct_AN];
        for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++) {
            Gh_AN = natn[ct_AN][h_AN];
            TNO2 = Total_NumOrbs[Gh_AN];
            for (i=0; i<TNO1; i++) {
                fread(OLP[ct_AN][h_AN][i],sizeof(double),TNO2,fp);
            }
        }
    }

    /****************************************************
            Overlap matrix with position operator x
    ****************************************************/

    for (ct_AN=1; ct_AN<=atomnum; ct_AN++) {
        TNO1 = Total_NumOrbs[ct_AN];
        for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++) {
            Gh_AN = natn[ct_AN][h_AN];
            TNO2 = Total_NumOrbs[Gh_AN];
            for (i=0; i<TNO1; i++) {
                fread(OLPpox[ct_AN][h_AN][i],sizeof(double),TNO2,fp);
            }
        }
    }

    /****************************************************
            Overlap matrix with position operator y
    ****************************************************/

    for (ct_AN=1; ct_AN<=atomnum; ct_AN++) {
        TNO1 = Total_NumOrbs[ct_AN];
        for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++) {
            Gh_AN = natn[ct_AN][h_AN];
            TNO2 = Total_NumOrbs[Gh_AN];
            for (i=0; i<TNO1; i++) {
                fread(OLPpoy[ct_AN][h_AN][i],sizeof(double),TNO2,fp);
            }
        }
    }

    /****************************************************
            Overlap matrix with position operator z
    ****************************************************/

    for (ct_AN=1; ct_AN<=atomnum; ct_AN++) {
        TNO1 = Total_NumOrbs[ct_AN];
        for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++) {
            Gh_AN = natn[ct_AN][h_AN];
            TNO2 = Total_NumOrbs[Gh_AN];
            for (i=0; i<TNO1; i++) {
                fread(OLPpoz[ct_AN][h_AN][i],sizeof(double),TNO2,fp);
            }
        }
    }

    /****************************************************
                      Density matrix
    ****************************************************/

    for (spin=0; spin<=SpinP_switch; spin++) {
        for (ct_AN=1; ct_AN<=atomnum; ct_AN++) {
            TNO1 = Total_NumOrbs[ct_AN];
            for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++) {
                Gh_AN = natn[ct_AN][h_AN];
                TNO2 = Total_NumOrbs[Gh_AN];
                for (i=0; i<TNO1; i++) {
                    fread(DM[spin][ct_AN][h_AN][i],sizeof(double),TNO2,fp);
                }
            }
        }
    }

    /****************************************************
        Solver
    ****************************************************/

    fread(i_vec,sizeof(int),1,fp);
    Solver = i_vec[0];

    /****************************************************
        ChemP
        Temp
    ****************************************************/

    fread(d_vec,sizeof(double),10,fp);
    ChemP  = d_vec[0];
    E_Temp = d_vec[1];
    dipole_moment_core[1] = d_vec[2];
    dipole_moment_core[2] = d_vec[3];
    dipole_moment_core[3] = d_vec[4];
    dipole_moment_background[1] = d_vec[5];
    dipole_moment_background[2] = d_vec[6];
    dipole_moment_background[3] = d_vec[7];
    Valence_Electrons = d_vec[8];
    Total_SpinS = d_vec[9];

    /****************************************************
        input file
    ****************************************************/

    fread(i_vec, sizeof(int), 1, fp);
    num_lines = i_vec[0];

    sprintf(makeinp,"temporal_12345.input");

    if ((fp_makeinp = fopen(makeinp,"w")) != NULL) {

#ifdef xt3
        setvbuf(fp_makeinp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

        for (i=1; i<=num_lines; i++) {
            fread(strg, sizeof(char), MAX_LINE_SIZE, fp);
            fprintf(fp_makeinp,"%s",strg);
        }

        fclose(fp_makeinp);
    }
    else {
        printf("error in making temporal_12345.input\n");
    }

}





void analysis_example(std::string  scfout_file)
{

    static int ct_AN,h_AN,Gh_AN,i,j,TNO1,TNO2;
    static int spin,Rn;

    static double *a;
//    static FILE *fp;

//    double Ebond[30],Es,Ep;

    int file_exist=read_scfout(scfout_file);

    /*
    for (spin=0; spin<=SpinP_switch; spin++){
      for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
        TNO1 = Total_NumOrbs[ct_AN];
        for (i=0; i<TNO1; i++) Ebond[i] = 0.0;
        for (i=0; i<TNO1; i++){
    for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
      Gh_AN = natn[ct_AN][h_AN];
      Rn = ncn[ct_AN][h_AN];
      TNO2 = Total_NumOrbs[Gh_AN];
            for (j=0; j<TNO2; j++){
              Ebond[i] += DM[spin][ct_AN][h_AN][i][j]*Hks[spin][ct_AN][h_AN][i][j];
      }
    }
          printf("spin=%2d ct_AN=%2d i=%2d Ebond=%15.12f\n",spin,ct_AN,i,Ebond[i]);
        }

        Es = Ebond[0] + Ebond[1] + Ebond[2] + Ebond[3] + Ebond[5] + Ebond[6] + Ebond[9] + Ebond[10];
        Ep = Ebond[4] + Ebond[7] + Ebond[8] + Ebond[11] + + Ebond[12];

        printf("%15.12f %15.12f\n",4.0*Es,4.0*Ep);

      }
    }
    */

    /**********************************************************************
     Example 1:

     Print the numbers of atoms

       atomnum:   the number of atoms in the total system
       Catomnum:  the number of atoms in the central region
       Latomnum:  the number of atoms in the left lead
       Ratomnum:  the number of atoms in the right lead

       grobal index of atom runs
       Catomnum -> Catomnum + Latomnum -> Catomnum + Latomnum + Ratomnum
    ***********************************************************************/

    /*
        printf("atomnum=%i\n", atomnum);
        printf("Catomnum=%i\n",Catomnum);
        printf("Latomnum=%i\n",Latomnum);
        printf("Ratomnum=%i\n",Ratomnum);
    */

    /**********************************************************************
     Example 2:

     Print Kohn-Sham Hamiltonian

        Hks[spin][ct_AN][h_AN][i][j]

        spin:  spin=0, up
               spin=1, down

        ct_AN: global index of atoms
        h_AN   local index of neighbouring atoms for the atom ct_AN
        i:     orbital index in the atom ct_AN
        j:     orbital index in the atom h_AN

     NOTE:

        For instance, if the basis specification of the atom ct_AN is s2p2,
        then the obital index runs in order of
                      s, s', px, py, pz, px', py', pz'

        Transformation of the local index h_AN to the grobal index Gh_AN
        is made as

                         Gh_AN = natn[ct_AN][h_AN];

        Also, the cell index is given by

                         Rn = ncn[ct_AN][h_AN];

        Each component l, m, or n (Rn = l*a + m*b + n*c) are given by

                         l = atv_ijk[Rn][1];
                         m = atv_ijk[Rn][2];
                         n = atv_ijk[Rn][3];
    ***********************************************************************/

    /*
        for (spin=0; spin<=SpinP_switch; spin++) {
            printf("\n\nKohn-Sham Hamiltonian spin=%i\n",spin);
            for (ct_AN=1; ct_AN<=atomnum; ct_AN++) {
                TNO1 = Total_NumOrbs[ct_AN];
                for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++) {
                    Gh_AN = natn[ct_AN][h_AN];
                    Rn = ncn[ct_AN][h_AN];
                    TNO2 = Total_NumOrbs[Gh_AN];
                    printf("glbal index=%i  local index=%i (grobal=%i, Rn=%i)\n",
                           ct_AN,h_AN,Gh_AN,Rn);
                    for (i=0; i<TNO1; i++) {
                        for (j=0; j<TNO2; j++) {
                            printf("%10.7f ",Hks[spin][ct_AN][h_AN][i][j]);
                        }
                        printf("\n");
                    }
                }
            }
        }
    */

    /**********************************************************************
     Example 3:

     Print overlap matrix

        OLP[ct_AN][h_AN][i][j]

        ct_AN: global index of atoms
        h_AN   local index of neighbouring atoms for the atom ct_AN
        i:     orbital index in the atom ct_AN
        j:     orbital index in the atom h_AN

     NOTE:

        For instance, if the basis specification of the atom ct_AN is s2p2,
        then the obital index runs in order of
                      s, s', px, py, pz, px', py', pz'

        Transformation of the local index h_AN to the grobal index Gh_AN
        is made as

                         Gh_AN = natn[ct_AN][h_AN];

        Also, the cell index is given by

                         Rn = ncn[ct_AN][h_AN];

        Each component l, m, or n (Rn = l*a + m*b + n*c) are given by

                         l = atv_ijk[Rn][1];
                         m = atv_ijk[Rn][2];
                         n = atv_ijk[Rn][3];
    ***********************************************************************/

    /*
        printf("\n\nOverlap matrix\n");
        for (ct_AN=1; ct_AN<=atomnum; ct_AN++) {
            TNO1 = Total_NumOrbs[ct_AN];
            for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++) {
                Gh_AN = natn[ct_AN][h_AN];
                Rn = ncn[ct_AN][h_AN];
                TNO2 = Total_NumOrbs[Gh_AN];
                printf("glbal index=%i  local index=%i (grobal=%i, Rn=%i)\n",
                       ct_AN,h_AN,Gh_AN,Rn);
                for (i=0; i<TNO1; i++) {
                    for (j=0; j<TNO2; j++) {
                        printf("%10.7f ",OLP[ct_AN][h_AN][i][j]);
                    }
                    printf("\n");
                }
            }
        }
    */


    /**********************************************************************
     Example 4:

     Print overlap matrix with position operator x

        OLPpox[ct_AN][h_AN][i][j]

        ct_AN: global index of atoms
        h_AN   local index of neighbouring atoms for the atom ct_AN
        i:     orbital index in the atom ct_AN
        j:     orbital index in the atom h_AN

     NOTE:

        For instance, if the basis specification of the atom ct_AN is s2p2,
        then the obital index runs in order of
                      s, s', px, py, pz, px', py', pz'

        Transformation of the local index h_AN to the grobal index Gh_AN
        is made as

                         Gh_AN = natn[ct_AN][h_AN];

        Also, the cell index is given by

                         Rn = ncn[ct_AN][h_AN];

        Each component l, m, or n (Rn = l*a + m*b + n*c) are given by

                         l = atv_ijk[Rn][1];
                         m = atv_ijk[Rn][2];
                         n = atv_ijk[Rn][3];
    ***********************************************************************/

    /*
        printf("\n\nOverlap matrix with position operator x\n");
        for (ct_AN=1; ct_AN<=atomnum; ct_AN++) {
            TNO1 = Total_NumOrbs[ct_AN];
            for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++) {
                Gh_AN = natn[ct_AN][h_AN];
                Rn = ncn[ct_AN][h_AN];
                TNO2 = Total_NumOrbs[Gh_AN];
                printf("glbal index=%i  local index=%i (grobal=%i, Rn=%i)\n",
                       ct_AN,h_AN,Gh_AN,Rn);
                for (i=0; i<TNO1; i++) {
                    for (j=0; j<TNO2; j++) {
                        printf("%10.7f ",OLPpox[ct_AN][h_AN][i][j]);
                    }
                    printf("\n");
                }
            }
        }
    */

    /**********************************************************************
     Example 5:

     Print density matrix

        DM[spin][ct_AN][h_AN][i][j]

        spin:  spin=0, up
               spin=1, down

        ct_AN: global index of atoms
        h_AN   local index of neighbouring atoms for the atom ct_AN
        i:     orbital index in the atom ct_AN
        j:     orbital index in the atom h_AN

     NOTE:

        For instance, if the basis specification of the atom ct_AN is s2p2,
        then the obital index runs in order of
                      s, s', px, py, pz, px', py', pz'

        Transformation of the local index h_AN to the grobal index Gh_AN
        is made as

                         Gh_AN = natn[ct_AN][h_AN];

        Also, the cell index is given by

                         Rn = ncn[ct_AN][h_AN];

        Each component l, m, or n (Rn = l*a + m*b + n*c) are given by

                         l = atv_ijk[Rn][1];
                         m = atv_ijk[Rn][2];
                         n = atv_ijk[Rn][3];
    ***********************************************************************/

    /*
        for (spin=0; spin<=SpinP_switch; spin++) {
            printf("\n\nDensity matrix spin=%i\n",spin);
            for (ct_AN=1; ct_AN<=atomnum; ct_AN++) {
                TNO1 = Total_NumOrbs[ct_AN];
                for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++) {
                    Gh_AN = natn[ct_AN][h_AN];
                    Rn = ncn[ct_AN][h_AN];
                    TNO2 = Total_NumOrbs[Gh_AN];

                    printf("glbal index=%i  local index=%i (grobal=%i, Rn=%i)\n",
                           ct_AN,h_AN,Gh_AN,Rn);

                    for (i=0; i<TNO1; i++) {
                        for (j=0; j<TNO2; j++) {
                            printf("%10.7f ",DM[spin][ct_AN][h_AN][i][j]);
                        }
                        printf("\n");
                    }
                }
            }
        }
    */


    /**********************************************************************
    Write Hk and Overlap matrix S  (modified by Jae-Hoon Sim, 20150611)
    DM =    Only for collinear calculation.
    **********************************************************************/

    if(file_exist==1) {
        double tol_local=1e-20;
        int l,m,n;
        FILE *fp_dmjh;
        FILE *fp_jh;//SJH Hk.HWR
        FILE *fp_jh2; //OverlapMatrix.HWR
        double H_up_down_r, H_up_down_i, H_spin_r, H_spin_i;
        double DM_up_down_r, DM_up_down_i, DM_spin_r, DM_spin_i;
        fp_jh = fopen("Hk.HWR","w");
        fp_dmjh = fopen("DM.HWR","w");
        for (ct_AN=1; ct_AN<=atomnum; ct_AN++) {
            TNO1 = Total_NumOrbs[ct_AN];
            fprintf(fp_jh, "Norbital of Atom = %d\n",  (2)*TNO1 );
        }
// (Hks[spin][ct_AN][h_AN][i][j] with Rn=ncn[c][h])  = <ci|H|hjR> =  <hjR | H | ci>*
// state1={ct_AN, TNO1, i}   , state2={h_AN, Gh_AN ,j, Rn=(l,m,n)}
        for (spin=0; spin<=SpinP_switch; spin++) {
            if(spin == 0 || spin ==1) {
                for (ct_AN=1; ct_AN<=atomnum; ct_AN++) {
                    TNO1 = Total_NumOrbs[ct_AN];
                    for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++) {
                        Gh_AN = natn[ct_AN][h_AN];
                        TNO2 = Total_NumOrbs[Gh_AN];
                        Rn = ncn[ct_AN][h_AN];
                        l = atv_ijk[Rn][1];
                        m = atv_ijk[Rn][2];
                        n = atv_ijk[Rn][3];
                        for (i=0; i<TNO1; i++) {
                            for (j=0; j<TNO2; j++) {
                                H_spin_r =  Hks[spin][ct_AN][h_AN][i][j]  ;
                                H_spin_i = iHks[spin][ct_AN][h_AN][i][j] ;
                                DM_spin_r = DM[spin][ct_AN][h_AN][i][j]  ;
                                DM_spin_i = 0.0;
                                /*Hamiltonian, colinear part*/
                                if( (std::pow(H_spin_r,2) + std::pow(H_spin_i,2)) > tol_local)
                                    fprintf(fp_jh, "%+d %+d %+d %d %d %d %d        %e   %e \n",
                                            l,   m,   n,
                                            Gh_AN, ct_AN,
                                            ((2)*j+spin)+1,   ((2)*i+spin)+1,
                                            (H_spin_r*27.2113845), -(H_spin_i*27.2113845)          ); /*Hartree To eV*/
                                /*DM*/
                                if( (std::pow(DM_spin_r,2) + std::pow(DM_spin_i,2)) > tol_local)
                                    fprintf(fp_dmjh, "%+d %+d %+d %d %d %d %d        %e   %e \n",
                                            l,   m,   n,
                                            Gh_AN, ct_AN,
                                            ((2)*j+spin)+1,   ((2)*i+spin)+1,
                                            (DM_spin_r),(DM_spin_i)          );
                                if(SpinP_switch==0) {
                                    if( (std::pow(H_spin_r,2) + std::pow(H_spin_i,2)) > tol_local)
                                        fprintf(fp_jh, "%+d %+d %+d %d %d %d %d        %e   %e \n",
                                                l,   m,   n,
                                                Gh_AN, ct_AN,
                                                ((2)*j+spin)+2,   ((2)*i+spin)+2,
                                                (H_spin_r*27.2113845), -(H_spin_i*27.2113845)          ); /*Hartree To eV*/
                                    if( (std::pow(DM_spin_r,2) + std::pow(DM_spin_i,2)) > tol_local)
                                        fprintf(fp_dmjh, "%+d %+d %+d %d %d %d %d        %e   %e \n",
                                                l,   m,   n,
                                                Gh_AN, ct_AN,
                                                ((2)*j+spin)+2,   ((2)*i+spin)+2,
                                                (DM_spin_r),(DM_spin_i)          );
                                }
                            }
                        }
                    }
                }
            }/*collinear part*/
        }/*SpinP_switch*/
        if(SpinP_switch ==3 ) {
// (Hks[2 or 3 ][ct_AN][h_AN][i][j] with Rn=ncn[c][h])  = <ci_up|H|hjR_down> =  <hjR_down | H | ci_up>*
            for (ct_AN=1; ct_AN<=atomnum; ct_AN++) {
                TNO1 = Total_NumOrbs[ct_AN];
                for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++) {
                    Gh_AN = natn[ct_AN][h_AN];
                    TNO2 = Total_NumOrbs[Gh_AN];
                    Rn = ncn[ct_AN][h_AN];
                    l = atv_ijk[Rn][1];
                    m = atv_ijk[Rn][2];
                    n = atv_ijk[Rn][3];
                    for (i=0; i<TNO1; i++) {
                        for (j=0; j<TNO2; j++) {
                            H_up_down_r = Hks[2][ct_AN][h_AN][i][j] ;
                            H_up_down_i = Hks[3][ct_AN][h_AN][i][j] + iHks[2][ct_AN][h_AN][i][j];
                            /*
                                                    DM_up_down_r = Hks[2][ct_AN][h_AN][i][j] ;
                                                    DM_up_down_i = Hks[3][ct_AN][h_AN][i][j] + iHks[2][ct_AN][h_AN][i][j];
                            */
                            if( (std::pow(H_up_down_r,2) + std::pow(H_up_down_i,2)) > tol_local)
                                fprintf(fp_jh, "%+d %+d %+d %d %d %d %d        %e   %e \n",
                                        l,   m,   n,
                                        Gh_AN, ct_AN,
                                        ((2)*j+1)+1,   ((2)*i+0)+1,
                                        ( H_up_down_r*27.2113845), -(H_up_down_i*27.2113845)          ); /*Hartree To eV*/
                            if( (std::pow(H_up_down_r,2) + std::pow(H_up_down_i,2)) > tol_local)
                                fprintf(fp_jh, "%+d %+d %+d %d %d %d %d        %e   %e \n",
                                        -l,   -m,   -n,
                                        ct_AN,Gh_AN,
                                        ((2)*i+0)+1,   ((2)*j+1)+1,
                                        ( H_up_down_r*27.2113845),   (H_up_down_i*27.2113845)          ); /*Hartree To eV*/
                        }
                    }
                }
            }
        }
        fclose(fp_jh);

        fp_jh2 = fopen("OverlapMatrix.HWR","w");
        for (ct_AN=1; ct_AN<=atomnum; ct_AN++) {
            TNO1 = Total_NumOrbs[ct_AN];
            for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++) {
                Gh_AN = natn[ct_AN][h_AN];
                TNO2 = Total_NumOrbs[Gh_AN];
                Rn = ncn[ct_AN][h_AN];
                l = atv_ijk[Rn][1];
                m = atv_ijk[Rn][2];
                n = atv_ijk[Rn][3];
                for (i=0; i<TNO1; i++) {
                    for (j=0; j<TNO2; j++) {
                        if( (std::pow(OLP[ct_AN][h_AN][i][j],2) ) > tol_local)
                            fprintf(fp_jh2, " %+d %+d %+d %d %d %d %d      %e    %e  \n",
                                    l,m,n,         Gh_AN, ct_AN,            j+1, i+1,            OLP[ct_AN][h_AN][i][j], 0.       );
                    }
                }
            }
        }
        fclose(fp_jh2);


        //OLPpox[A1][A2][i][j] =  <A1,i| x-r_A1 | A2, j R> in Bhor unit. Here r_A1 = position of the atom A1.
        FILE *fp_jh3; //OverlapMatrix_x.HWR
        FILE *fp_jh4; //OverlapMatrix_y.HWR
        FILE *fp_jh5; //OverlapMatrix_z.HWR
        fp_jh3 = fopen("OverlapMatrix_x.HWR","w");
        fp_jh4 = fopen("OverlapMatrix_y.HWR","w");
        fp_jh5 = fopen("OverlapMatrix_z.HWR","w");
        for (ct_AN=1; ct_AN<=atomnum; ct_AN++) {
            TNO1 = Total_NumOrbs[ct_AN];
            for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++) {
                Gh_AN = natn[ct_AN][h_AN];
                TNO2 = Total_NumOrbs[Gh_AN];
                Rn = ncn[ct_AN][h_AN];
                l = atv_ijk[Rn][1];
                m = atv_ijk[Rn][2];
                n = atv_ijk[Rn][3];
                for (i=0; i<TNO1; i++) {
                    for (j=0; j<TNO2; j++) {

                        double x_expec_inANG = ( OLPpox[ct_AN][h_AN][i][j] +  ( Gxyz[ct_AN][1])*OLP[ct_AN][h_AN][i][j]  )*0.529177249;
                        double y_expec_inANG = ( OLPpoy[ct_AN][h_AN][i][j] +  ( Gxyz[ct_AN][2])*OLP[ct_AN][h_AN][i][j]  )*0.529177249;
                        double z_expec_inANG = ( OLPpoz[ct_AN][h_AN][i][j] +  ( Gxyz[ct_AN][3])*OLP[ct_AN][h_AN][i][j]  )*0.529177249;
                        //x_expec_inANG = <i|x|jR> = <jR|x|i>*
                        if( (std::pow(x_expec_inANG,2) ) > tol_local)
                            fprintf(fp_jh3, " %+d %+d %+d %d %d %d %d      %e    %e  \n",
                                    l,m,n,         Gh_AN,ct_AN,            j+1, i+1,            x_expec_inANG, 0.       );
                        if( (std::pow(y_expec_inANG,2) ) > tol_local)
                            fprintf(fp_jh4, " %+d %+d %+d %d %d %d %d      %e    %e  \n",
                                    l,m,n,         Gh_AN,ct_AN,            j+1, i+1,            y_expec_inANG, 0.       );
                        if( (std::pow(z_expec_inANG,2) ) > tol_local)
                            fprintf(fp_jh5, " %+d %+d %+d %d %d %d %d      %e    %e  \n",
                                    l,m,n,         Gh_AN,ct_AN,            j+1, i+1,            z_expec_inANG, 0.       );
                    }
                }
            }
        }
        fclose(fp_jh3);
        fclose(fp_jh4);
        fclose(fp_jh5);
    }
}
