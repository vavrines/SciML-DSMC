//******************************************************************************
//Title:
//        Direct Simulation Monte Carlo (DSMC) program for two-dimensional flows
//Version:
//        Version2 2.3.05
//History:
//        Written by Bijan Gosaheyshi, 7/2021.
//******************************************************************************

#define _XOPEN_SOURCE /* ALWAYS BEFORE the include statement */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <tgmath.h>
#include <float.h>
#include <time.h>
#include "dsmcBasicInput.h"
//--------------Function declarations-----------------------//
void dsmcINITIALIZATION();
void initAftRead();
void OUTPUT();
void TRAckParticle();
void TRAckParticle_Groups();
void INIT();
void INDEX();
void Move();
void Collision();
void SELECT(int N, int i, int j, int *L, int *M, long double *vr);
void ELASTIC(int M1, int M2, long double VR);
void Read();
void WriteOut();
void SAMPI();
void ContinueSample();
void Sample();
void ENTER();
void RVELC(long double *U, long double *V, long double VMP);
void Remove(int *N);
void REFLECT2D(int n, int k, long double X, long double Y, int m, int dir);
void init_by_array64(unsigned long long init_key[],
                     unsigned long long key_length);
double genrand64_real2(void);

long double **AllocateArray2D(int dim1, int dim2);
long double ***AllocateArray3D(int dim1, int dim2, int dim3);
int **AllocateArray2DINTEGER(int dim1, int dim2);
long double GAM(long double X);
long double RF(int X);
long double ERF(long double S);

//--------------------------------------------------------
#ifndef M_PI
#define M_PI (long double)3.14159265358979323846264338327950288
#endif
#define SPI (long double)sqrtl(M_PI)

#define BOLTZ (long double)1.380622E-23

double drand48(void);
//////---Global Variable declaration--------------------
clock_t startTime, endTime;

struct DSMC1
{
    int NM, NPR, SELT, MNM, MNC, NSMP, IFSurf, NX, NY, RandSeed;                                                  //9
    long unsigned int MOVT, NCOL, REPTN;                                                                          //3
    long double FND, FNUM, CharLength, landa, deltax, deltay, pp1, CW, CH, FW, FH, CD, VMP, CVR, C_Soundi, Machi; //16
    long double DTM, DDTM, SurfDistL, SurfHeightL;                                                                //4
    long double SPM0, SPM1, SPM2, SPM3, SPM4, SPM5;                                                               //6
    long double TIME, TIMI, SEPT;                                                                                 //3
    double time_taken;                                                                                            //1
};
//
struct DSMC2
{
    int *IB, *IR;
    long double *CB, *sp, *TSURF, *BME, *BMR;
    int **IC, **CtIJ, **IJtC, **ICol;
    long double **CT, **CC, **PP, **PV, **CS, **DCol;
    long double ***CG, ***CCG;
};
struct DSMC1 DST1;
struct DSMC2 DST2;

/*--CCG is for collisions
CCG[0] is the maximum value of (relative speed)*(coll. cross-section)
CCG[1] is the remainder when the collision selection number is rounded
IC[0-1][M]: information on the particles in cell M 
IC[0]: (start address -1) of the particle numbers in array IR
IC[1]: the number of particles in the cell M
*/
// Iperiod[0-1] X-Y directions =0Yes, 1 No
// SPM0: collision cross section eq(1-35)
// SPM1: the reference temperature
// SPM2: the viscosity-temperature power law
// SPM3: the reciprocal of the VSS scattering parameter
// SPM4: the reduced mass
// SPM5: the Gamma function of (5/2 - viscosity-temperature power law)
/*
*--TIME time
*--NM is the number of molecules
*--NPR the number of output/restart file update cycles
*--NCOL is the total number of collisions
*--MOVT the total number of molecular moves
*--SELT the total number of pair selections
*--SEPT the sum of collision pair separations
*/
/*
*--SP(N) information on gas
*----N=0 the reference cross-section (diameter in the data)
*----N=1 the reference temperature
*----N=2 the viscosity-temperature power law
*----N=3 the reciprocal of the VSS scattering parameter
*----N=4 the molecular mass
*/
/*
*--CG(N,X,Y) is the geometry related information on cell X Y
N=0 the minimum x coordinate
N=1 the maximum x coordinate
N=2 the cell width X
N=3 the minimum y coordinate
N=4 the maximum y coordinate
N=5 the cell width Y
N=6 the cell center X
N=7 the cell center Y
*--CCG(N,X,Y) is for collisions in cell X Y
*----N=0 is the maximum value of (relative speed)*(coll. cross-section)
*----N=1 is the remainder when the selection number is rounded
CS
N=0 number sum
N=1,2,3 sum of u,v,w
N=4,5,6 sum of u*u,v*v,w*w
*/
//IB=1 Stream; IB=2 Surface
//long double ***feq;    //long double feq[NX][NY][NPOP];
//CtIJ[i][j] convert column grid to NCELL
//IJtC[0-1][NCELL] convert NCELL to column i,j grid
/*DCol, ICol: information about Collision cells
DCol double numbers, ICol integer numbers
DCol[0]:separation distances;
ICol[0]:Number of collisions;
*/
///////////////////////////////////////////////////

#ifndef _DSMCINCLUDE_H

#define _DSMCINCLUDE_H

#endif /* _DSMCINCLUDE_H */
       ///////////////////////////////////////////////////