//******************************************************************************
//Title:
//        Direct Simulation Monte Carlo (DSMC) program for two-dimensional flows
//Version:
//        Version2 2.3.05
//History:
//        Written by Bijan Gosaheyshi, 7/2021.
//******************************************************************************

#include "dsmcINCLUDE.h"

////////////////////////////////////////////////////////////////////
void dsmcINITIALIZATION()
{
    int MM, NCELL, Nj, Ni;
    long double A, REM, SC;
    //srand48(time(NULL));
    srand48(time(0));
    //Make zero variables
    startTime = clock(); //start counting time
    DST1.NSMP = 0;
    DST1.IFSurf = 0;
    if ((SurfDist >= 0) && (SurfDist <= 1.0))
    {
        DST1.IFSurf = 1;
    } //if there is additional surface inside
    ////
    DST2.sp = malloc(5 * sizeof(long double));
    DST2.CB = malloc(2 * sizeof(long double));
    DST2.TSURF = malloc((4 + DST1.IFSurf) * sizeof(long double));
    DST2.IB = malloc((4 + DST1.IFSurf) * sizeof(int));
    DST2.BME = malloc(4 * sizeof(long double));
    DST2.BMR = malloc(4 * sizeof(long double));
    DST2.CB[0] = cb1;
    DST2.CB[1] = cb2;
    DST2.IB[0] = ib1;
    DST2.IB[1] = ib2;
    DST2.IB[2] = ib3;
    DST2.IB[3] = ib4;
    if (DST1.IFSurf == 1)
    {
        DST2.IB[4] = ib5;
        DST2.TSURF[4] = Tsurf5;
    }
    DST1.SurfHeightL = SurfHeight * DST2.CB[1];
    DST1.SurfDistL = SurfDist * DST2.CB[0];
    DST2.TSURF[0] = Tsurf1;
    DST2.TSURF[1] = Tsurf2;
    DST2.TSURF[2] = Tsurf3;
    DST2.TSURF[3] = Tsurf4;
    DST2.BME[0] = 0.;
    DST2.BME[1] = 0.;
    DST2.BME[2] = 0.;
    DST2.BME[3] = 0.;
    DST2.BMR[0] = 0.;
    DST2.BMR[1] = 0.;
    DST2.BMR[2] = 0.;
    DST2.BMR[3] = 0.;

    DST1.CharLength = DST2.CB[0];

    if (charL == 0)
    {
        DST1.CharLength = DST2.CB[0];
    }
    else if (charL == 1)
    {
        DST1.CharLength = DST2.CB[1];
    }

    //if (DST2.CB[1]<DST1.CharLength) DST1.CharLength=DST2.CB[1];//select smaller length as characteristic length
    //
    DST2.sp[0] = sp1D;
    DST2.sp[1] = sp2T;
    DST2.sp[2] = sp3VT;
    DST2.sp[3] = sp4ReVSS;
    DST2.sp[4] = sp5M;
    /*sp[] information on gas: 
0 the reference cross-section (diameter in the data),
1 the reference temperature, 2 the viscosity-temperature power law,
3 the reciprocal of the VSS scattering parameter, 4 the molecular mass*/

    DST1.SPM0 = M_PI * DST2.sp[0] * DST2.sp[0];                      //collision cross section (1-35)
    DST1.SPM1 = DST2.sp[1];                                          //the reference temperature
    DST1.SPM2 = DST2.sp[2];                                          //the viscosity-temperature power law
    DST1.SPM3 = DST2.sp[3];                                          //the reciprocal of the VSS scattering parameter
    DST1.SPM4 = DST2.sp[4] * DST2.sp[4] / (DST2.sp[4] + DST2.sp[4]); //the reduced mass is defined in eqn (2.7)
    DST1.SPM5 = GAM(2.5 - DST1.SPM2);                                //the Gamma function of (5/2 - viscosity-temperature power law)
    //calcualte least required number of divisions based on the Mean_Free_Path(landa) or knundsen number
    DST1.landa = DST1.CharLength * Kn; //Mean free path based on simulation Knudsen number

    DST1.deltax = DST1.landa / 5.;
    DST1.NX = DST2.CB[0] / DST1.deltax + 0.99;

    DST1.deltay = DST1.landa / 5.;
    DST1.NY = DST2.CB[1] / DST1.deltay + 0.99;

    printf("Before, NX= %d NY= %d\n", DST1.NX, DST1.NY);

    if (DST1.NX < MinDevisionX)
    {
        DST1.NX = MinDevisionX;
    }
    if (DST1.NY < MinDevisionY)
    {
        DST1.NY = MinDevisionY;
    }
    printf("After, NX= %d NY= %d\n", DST1.NX, DST1.NY);

    DST1.deltax = DST2.CB[0] / ((float)DST1.NX);
    DST1.deltay = DST2.CB[1] / ((float)DST1.NY);

    DST1.MNC = DST1.NX * DST1.NY;
    DST1.MNM = DST1.MNC * ppcel * 1.4;
    DST1.CW = DST2.CB[0] / DST1.NX;
    DST1.CH = DST2.CB[1] / DST1.NY;
    DST1.FW = DST2.CB[0];
    DST1.FH = DST2.CB[1];
    DST1.CD = sqrtl(DST1.CW * DST1.CW + DST1.CH * DST1.CH);
    //////Allocatoin of DSMC Arrays///
    DST2.CG = AllocateArray3D(8, DST1.NX, DST1.NY);
    DST2.CCG = AllocateArray3D(2, DST1.NX, DST1.NY);
    DST2.CT = AllocateArray2D(DST1.NX, DST1.NY);
    DST2.CC = AllocateArray2D(DST1.NX, DST1.NY);

    DST2.CS = AllocateArray2D(19, DST1.MNC);
    DST2.PP = AllocateArray2D(3, DST1.MNM);
    DST2.PV = AllocateArray2D(3, DST1.MNM);
    DST2.IC = AllocateArray2DINTEGER(2, DST1.MNC);
    DST2.IJtC = AllocateArray2DINTEGER(DST1.NX, DST1.NY);
    DST2.CtIJ = AllocateArray2DINTEGER(2, DST1.MNC);
    DST2.IR = malloc(DST1.MNM * sizeof(int));
    DST2.DCol = AllocateArray2D(1, DST1.MNC);
    DST2.ICol = AllocateArray2DINTEGER(1, DST1.MNC);
    //////Allocatoin of GRID Arrays///
    for (int j = 0; j < DST1.NY; j++)
    {
        for (int i = 0; i < DST1.NX; i++)
        {
            NCELL = j * DST1.NX + (i); //number of cells from 0:NX.NY-1
            DST2.IJtC[i][j] = NCELL;
            DST2.CtIJ[0][NCELL] = i;
            DST2.CtIJ[1][NCELL] = j;
        }
    }

    //////Allocatoin of DSMC Arrays///
    for (int i = 0; i < DST1.MNC; i++)
    {
        DST2.IC[0][i] = 0;
        DST2.IC[1][i] = 0;
        DST2.DCol[0][i] = 0.0; //separation collision distance
        DST2.ICol[0][i] = 0;   //number of accepted collisions
        for (int j = 0; j <= 18; j++)
        {
            DST2.CS[j][i] = 0.;
        }
    }
    for (int i = 0; i < DST1.MNM; i++)
    {
        DST2.IR[i] = 0;
        for (int k = 0; k < 3; k++)
        {
            DST2.PP[k][i] = -9000000;
            //
            DST2.PV[k][i] = -9000000;
        }
    }
    ///zero making

    DST1.FND = 1 / (sqrtl(2) * M_PI * DST2.sp[0] * DST2.sp[0] * DST1.CharLength * Kn);
    DST1.FNUM = DST1.FND * (DST2.CB[0] / DST1.NX) * (DST2.CB[1] / DST1.NY) / (ppcel);
    printf("FND= %Le DST1.FNUM= %Le\n", DST1.FND, DST1.FNUM);

    //pp1=(FND*(cb[0]/NX)*(cb[1]/NY))/FNUM;
    //printf("ppcel= %Lf pp1= %Lf\n",ppcel,pp1);
    printf("FND= %Le FNUM= %Le %Le %Le %Le %Le %Le\n", DST1.FND, DST1.FNUM, DST2.sp[0], DST2.sp[1], DST2.sp[2], DST2.sp[3], DST2.sp[4]);
    printf("CB= %Le, %Le, %f, %Le\n", DST2.CB[1], DST2.CB[0], M_PI, DST1.CharLength);
    DST1.DTM = 0.19 * (DST1.CharLength / DST1.NX) / sqrtl(2. * BOLTZ * FTMP / DST2.sp[4]); //time by VMP
    DST1.DDTM = 1e6;
    A = sqrtl((VFX + WallvelX) * (VFX + WallvelX) + (VFY + WallvelY) * (VFY + WallvelY));
    DST1.C_Soundi = sqrtl(5. / 3. * BOLTZ * FTMP / DST2.sp[4]); //initial sound speed:sqrtl(gama k T / m):322.6752
    DST1.Machi = A / DST1.C_Soundi;                             //initial mach number
    if (A > 0)
        DST1.DDTM = 0.3 * (DST1.CharLength / DST1.NX) / A; //time by max vel

    printf("Time Values, by VMP = %Le, by Wallvel= %Le\n", DST1.DTM, DST1.DDTM);
    if (DST1.DDTM < DST1.DTM)
    {
        DST1.DTM = DST1.DDTM; // select minimum DTM
    }

    /// initialize zero//
    DST1.TIME = 0.;
    DST1.TIMI = 0.;
    DST1.NM = 0;
    DST1.NPR = 0;
    DST1.REPTN = 0;
    DST1.NCOL = 0;
    DST1.MOVT = 0.;
    DST1.SELT = 0.;
    DST1.SEPT = 0.;
    //Initial Cell  evaluation

    for (int i = 0; i < DST1.NX; i++)
    {
        for (int j = 0; j < DST1.NY; j++)
        {
            //x
            DST2.CG[0][i][j] = (i)*DST1.CW;
            DST2.CG[1][i][j] = (i + 1.) * DST1.CW;
            DST2.CG[2][i][j] = DST2.CG[1][i][j] - DST2.CG[0][i][j];
            //y
            DST2.CG[3][i][j] = (j)*DST1.CH;
            DST2.CG[4][i][j] = (j + 1.) * DST1.CH;
            DST2.CG[5][i][j] = DST2.CG[4][i][j] - DST2.CG[3][i][j];
            DST2.CG[6][i][j] = 0.5 * (DST2.CG[0][i][j] + DST2.CG[1][i][j]);
            DST2.CG[7][i][j] = 0.5 * (DST2.CG[3][i][j] + DST2.CG[4][i][j]);
            DST2.CT[i][j] = FTMP;
            DST2.CC[i][j] = DST2.CG[2][i][j] * DST2.CG[5][i][j];
            // printf("cg position (%d,%d): %Le\n",i,j,CG[5][i][j]);
            DST2.CCG[0][i][j] = DST1.SPM0 * 300. * sqrtl(FTMP / 300.);
            DST2.CCG[1][i][j] = RF(Irandom);
        }
    }
    DST1.NM = -1;
    REM = 0;
    DST1.VMP = sqrtl(2. * BOLTZ * FTMP / DST2.sp[4]);
    ////Group particle definition////////////////////
    int detG;
    long double t_group;
#define groups 5
    double temp_group[groups] = {16, 160, 480, 960, 1920};
    double vmp_group[groups];
    for (int k = 0; k < groups; k++)
    {
        vmp_group[k] = sqrtl(2. * BOLTZ * temp_group[k] / DST2.sp[4]);
        t_group = 0.19 * (DST1.CharLength / DST1.NX) / vmp_group[k];
        printf("Group %d temp(k): %lf  vmp(m/s): %lf\n", k, temp_group[k], vmp_group[k]);
        if (DST1.VMP < vmp_group[k])
            DST1.VMP = vmp_group[k];
        if (t_group < DST1.DTM)
            DST1.DTM = t_group;
    }
#define div_group (int)(ppcel / groups)
    printf("Group length %d from tot particles %Lf\n", div_group, ppcel);
    printf("dt %Lf dt_group %Lf\n", DST1.DTM, t_group);

    ////////////////////Group particle definition///
    for (int j = 0; j < DST1.NY; j++)
    {
        for (int i = 0; i < DST1.NX; i++)
        {
            A = (DST1.FND * DST2.CC[i][j]) / (DST1.FNUM) + REM;
            MM = A;
            REM = (A - MM);
            //printf("A,MM,REM, %Lf %d %Lf\n",A,MM,REM);
            if (i == DST1.NX - 1 && j == DST1.NY - 1)
                MM = round(A);
            if (MM > 0)
            {
                for (int k = 0; k < MM; k++)
                {
                    if (DST1.NM > DST1.MNM)
                    {
                        printf("ExceedParticle %d > %d:: increase MNM\n", DST1.NM, DST1.MNM);
                        getchar();
                    }
                    if (DST1.NM <= DST1.MNM)
                    {
                        DST1.NM = DST1.NM + 1;
                        DST2.PP[0][DST1.NM] = DST2.CG[0][i][j] + RF(Irandom) * DST2.CG[2][i][j];
                        DST2.PP[1][DST1.NM] = DST2.CG[3][i][j] + RF(Irandom) * DST2.CG[5][i][j];
                        NCELL = j * DST1.NX + (i); //number of cells from 0:NX.NY-1
                        DST2.PP[2][DST1.NM] = NCELL;
                        Nj = (int)(NCELL / DST1.NX); //number of row from     0:NY-1
                        Ni = NCELL - (Nj * DST1.NX); //number of column from  0:NX-1

                        //determine group of particle
                        detG = (int)(DST1.NM / div_group);
                        DST1.VMP = vmp_group[detG];
                        //printf("Particle %d @group %d VMP %Lf\n", DST1.NM, detG, DST1.VMP);
                        //determine group of particle

                        for (int kk = 0; kk < 3; kk++)
                        {
                            RVELC(&DST2.PV[kk][DST1.NM], &A, DST1.VMP);
                        }
                        DST2.PV[0][DST1.NM] = DST2.PV[0][DST1.NM] + VFX;
                        DST2.PV[1][DST1.NM] = DST2.PV[1][DST1.NM] + VFY;
                    }
                }
            }
        }
    }
    DST1.NM = DST1.NM + 1;

    //so it means that:
    //now particle number varies from 0 till DST1.NM-1

    printf(" %d Particles are generated out of %d MNM, for %Lf ppcel in %d Cells\n", DST1.NM, DST1.MNM, ppcel, DST1.MNC);

    //calculate the number of particles that enter at each time step
    //4 borders: Xbtm/upp Ybtm/upp IB0123
    for (int k = 0; k < 4; k++)
    {
        if (DST2.IB[k] == 1)
        {
            printf("\n Side %d stream set up\n", k);
            DST1.VMP = sqrtl(2. * BOLTZ * FTMP / DST2.sp[4]);
            if (k == 0)
                SC = VFX / DST1.VMP;
            if (k == 1)
                SC = -VFX / DST1.VMP;
            if (k == 2)
                SC = VFY / DST1.VMP;
            if (k == 3)
                SC = -VFY / DST1.VMP;
            if (fabsl(SC) < 10.1)
            {
                A = (exp(-SC * SC) + SPI * SC * (1. + ERF(SC))) / (2. * SPI);
            }
            if (SC > 10.)
            {
                A = SC;
            }
            if (SC < -10.)
            {
                A = 0.;
            } //A is the non-dimensional flux of eqn (4.22)

            if ((k == 0) || (k == 1))
            {
                DST2.BME[k] = DST1.FND * A * DST1.VMP * DST1.DTM * DST1.FH / DST1.FNUM;
            }
            else
            {
                DST2.BME[k] = DST1.FND * A * DST1.VMP * DST1.DTM * DST1.FW / DST1.FNUM;
            }
            printf("entering mols %Lf\n", DST2.BME[k]);
        }
    }
    printf("\nEnd DSMC initialize\n \n");
    if (DST1.NPR == 0)
    {
        printf("\nWrite the initialization state 0\n \n");
        INDEX();
        Sample();
        OUTPUT();
    }

    FILE *fp;
    char buf[65];
    sprintf(buf, "../Output/DescrpFlow.dat");
    fp = fopen(buf, "w");
    long double p0;
    p0 = DST1.FND * BOLTZ * FTMP;

    fprintf(fp, "VARIABLES = \"FND\" \"DTM\" \"NM\" \"FNUM\" \"P0\" \"T0\" \"NX\" \"NY\" \"CBX\" \"CBY\"\n");
    fprintf(fp, "%Le %Le %d %Le %Lf %Lf %d %d %Le %Le\n", DST1.FND, DST1.DTM, DST1.NM, DST1.FNUM, p0, FTMP, DST1.NX, DST1.NY, DST2.CB[0], DST2.CB[1]);
    fclose(fp);
}
////////////////////////////////////////////////////////////////////

//function to allocate 3D array
long double ***AllocateArray3D(int dim1, int dim2, int dim3)
{

    int i, j, k;

    long double ***p = (long double ***)malloc(dim1 * sizeof(long double **));

    for (i = 0; i < dim1; i++)
    {
        p[i] = (long double **)malloc(dim2 * sizeof(long double *));
        for (j = 0; j < dim2; j++)
            p[i][j] = (long double *)malloc(dim3 * sizeof(long double));
    }

    return p;
}

//function to allocate 2D array

long double **AllocateArray2D(int dim1, int dim2)
{

    int i, j;

    long double **p = (long double **)malloc(dim1 * sizeof(long double *));

    for (i = 0; i < dim1; i++)
        p[i] = (long double *)malloc(dim2 * sizeof(long double));

    return p;
}

///////
//function to allocate 2D array integer

int **AllocateArray2DINTEGER(int dim1, int dim2)
{

    int i, j;

    int **p = (int **)malloc(dim1 * sizeof(int *));

    for (i = 0; i < dim1; i++)
        p[i] = (int *)malloc(dim2 * sizeof(int));

    return p;
}

///////

void INIT()
{
    printf("\n Initialize all the variables, for a new simulation\n");
}
///////
void Read()
{
    int i, j;
    printf("Read the information of an existing simulation\n");
    FILE *fptr;
    if ((fptr = fopen("./dsmcBinary.bin", "rb+")) == NULL)
    {
        printf("Error! opening file");
        exit(1);
    }
    ///DST1
    fread(&DST1, sizeof(struct DSMC1), 1, fptr);
    ///DST2
    fread(DST2.IB, (4 + DST1.IFSurf) * sizeof(int), 1, fptr);
    fread(DST2.IR, (DST1.MNM) * sizeof(int), 1, fptr);
    fread(DST2.CB, 2 * sizeof(long double), 1, fptr);
    fread(DST2.sp, 5 * sizeof(long double), 1, fptr);
    fread(DST2.TSURF, (4 + DST1.IFSurf) * sizeof(long double), 1, fptr);
    fread(DST2.BME, (4) * sizeof(long double), 1, fptr);
    fread(DST2.BMR, (4) * sizeof(long double), 1, fptr);
    for (i = 0; i <= 1; i++)
    {
        fread(DST2.IC[i], DST1.MNC * sizeof(int), 1, fptr);
        fread(DST2.CtIJ[i], DST1.MNC * sizeof(int), 1, fptr);
    }
    for (i = 0; i < DST1.NX; i++)
    {
        fread(DST2.IJtC[i], DST1.NY * sizeof(int), 1, fptr);
        fread(DST2.CT[i], DST1.NY * sizeof(long double), 1, fptr);
        fread(DST2.CC[i], DST1.NY * sizeof(long double), 1, fptr);
    }
    fread(DST2.ICol[0], DST1.MNC * sizeof(int), 1, fptr);
    fread(DST2.DCol[0], DST1.MNC * sizeof(long double), 1, fptr);
    for (i = 0; i <= 2; i++)
    {
        fread(DST2.PP[i], DST1.MNM * sizeof(long double), 1, fptr);
        fread(DST2.PV[i], DST1.MNM * sizeof(long double), 1, fptr);
    }
    for (i = 0; i < 19; i++)
    {
        fread(DST2.CS[i], DST1.MNC * sizeof(long double), 1, fptr);
    }

    for (i = 0; i < 8; i++)
    {
        for (j = 0; j < DST1.NX; j++)
        {
            fread(DST2.CG[i][j], DST1.NY * sizeof(long double), 1, fptr);
        }
    }
    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < DST1.NX; j++)
        {
            fread(DST2.CCG[i][j], DST1.NY * sizeof(long double), 1, fptr);
        }
    }
    //DST2 Finished
    fclose(fptr);
    //check whether required to call initAftRead function:
    if (Init_AftRead == 0)
        initAftRead();

    printf("\t %Le\n", DST2.PP[0][DST1.NM - 1]);
    //printf("\t %Lf %Lf %Lf\n", DST2.PP[0][DST1.NM-1],DST2.PP[1][DST1.NM-1],DST2.PP[2][DST1.NM-1]);
    printf("\t FND: %Le\t MNC: %d\t CT[1][2]: %Le, IC[0][MNC-1]: %d IR[NM-1]: %d\n", DST1.FND, DST1.MNC, DST2.CT[1][2], DST2.IC[0][DST1.MNC - 1], DST2.IR[DST1.NM - 1]);
}
///////
void WriteOut()
{
    int i, j;
    //printf("\n WriteOut the information of an existing simulation\n");
    FILE *fptr;
    if ((fptr = fopen("./dsmcBinary.bin", "wb+")) == NULL)
    {
        printf("Error! opening file");
        // Program exits if the file pointer returns NULL.
        exit(1);
    }
    ///DST1
    fwrite(&DST1, sizeof(struct DSMC1), 1, fptr);
    ///DST2
    fwrite(DST2.IB, (4 + DST1.IFSurf) * sizeof(int), 1, fptr);
    fwrite(DST2.IR, (DST1.MNM) * sizeof(int), 1, fptr);
    fwrite(DST2.CB, 2 * sizeof(long double), 1, fptr);
    fwrite(DST2.sp, 5 * sizeof(long double), 1, fptr);
    fwrite(DST2.TSURF, (4 + DST1.IFSurf) * sizeof(long double), 1, fptr);
    fwrite(DST2.BME, (4) * sizeof(long double), 1, fptr);
    fwrite(DST2.BMR, (4) * sizeof(long double), 1, fptr);
    for (i = 0; i <= 1; i++)
    {
        fwrite(DST2.IC[i], DST1.MNC * sizeof(int), 1, fptr);
        fwrite(DST2.CtIJ[i], DST1.MNC * sizeof(int), 1, fptr);
    }
    for (i = 0; i < DST1.NX; i++)
    {
        fwrite(DST2.IJtC[i], DST1.NY * sizeof(int), 1, fptr);
        fwrite(DST2.CT[i], DST1.NY * sizeof(long double), 1, fptr);
        fwrite(DST2.CC[i], DST1.NY * sizeof(long double), 1, fptr);
    }
    fwrite(DST2.ICol[0], DST1.MNC * sizeof(int), 1, fptr);
    fwrite(DST2.DCol[0], DST1.MNC * sizeof(long double), 1, fptr);
    for (i = 0; i <= 2; i++)
    {
        fwrite(DST2.PP[i], DST1.MNM * sizeof(long double), 1, fptr);
        fwrite(DST2.PV[i], DST1.MNM * sizeof(long double), 1, fptr);
    }
    for (i = 0; i < 19; i++)
    {
        fwrite(DST2.CS[i], DST1.MNC * sizeof(long double), 1, fptr);
    }

    for (i = 0; i < 8; i++)
    {
        for (j = 0; j < DST1.NX; j++)
        {
            fwrite(DST2.CG[i][j], DST1.NY * sizeof(long double), 1, fptr);
        }
    }
    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < DST1.NX; j++)
        {
            fwrite(DST2.CCG[i][j], DST1.NY * sizeof(long double), 1, fptr);
        }
    }
    fclose(fptr);
}
///////
void SAMPI()
{
    int L, N, k, j, M;
    long double A;
    if (DST1.NPR <= NPS)
        printf("\n Initialize samples bc (NPR:%d<NPS:%d)\n", DST1.NPR, NPS);

    DST1.NSMP = 0; //sample number increment
    DST1.NCOL = 0;
    DST1.MOVT = 0;
    DST1.TIMI = DST1.TIME;

    for (N = 0; N < DST1.MNC; N++)
    { //for each cell N
        DST2.ICol[0][N] = 0;
        DST2.DCol[0][N] = 0.0;
        for (k = 0; k <= 18; k++)
        {
            DST2.CS[k][N] = 0;
            //k=0 sample number sum
            //k=1 sample sum of u
            //k=2 sample sum of v
            //k=3 sample sum of w
            //k=4 sample sum of uu
            //k=5 sample sum of vv
            //k=6 sample sum of ww
            //k=7 sample sum of uv
            //k=8 sample sum of uw
            //k=9 sample sum of vw

            //k=10 sample sum of uuu
            //k=11 sample sum of uuv
            //k=12 sample sum of uuw

            //k=13 sample sum of vvu
            //k=14 sample sum of vvv
            //k=15 sample sum of vvw

            //k=16 sample sum of wwu
            //k=17 sample sum of wwv
            //k=18 sample sum of www

        } //k loop
    }     //MNC loop
}
///////
void ContinueSample()
{
    printf("\n Store the current sampling (NSMP=%d), and continue\n", DST1.NSMP);
    printf("\n Please Confirm, otherwise stop\n");
}

///////
void Sample()
{
    int L, N, k, j, M;
    long double A;
    //printf("\n Sample the flow for each cell\n");

    DST1.NSMP = DST1.NSMP + 1; //sample number increment

    for (N = 0; N < DST1.MNC; N++)
    {                      //for each cell N
        L = DST2.IC[1][N]; // the number of particles at each cell L
        if (L > 0)
        { //if there are any particle
            for (j = 1; j < L + 1; j++)
            {                          //for each particle j
                k = DST2.IC[0][N] + j; //get the number of particle
                k = k - 1;
                M = DST2.IR[k];                                                //get the number of particle based on cross referencig
                DST2.CS[0][N] = DST2.CS[0][N] + 1;                             //sample number sum
                DST2.CS[1][N] = DST2.CS[1][N] + DST2.PV[0][M];                 //sample sum of u
                DST2.CS[2][N] = DST2.CS[2][N] + DST2.PV[1][M];                 //sample sum of v
                DST2.CS[3][N] = DST2.CS[3][N] + DST2.PV[2][M];                 //sample sum of w
                DST2.CS[4][N] = DST2.CS[4][N] + DST2.PV[0][M] * DST2.PV[0][M]; //sample sum of uu
                DST2.CS[5][N] = DST2.CS[5][N] + DST2.PV[1][M] * DST2.PV[1][M]; //sample sum of vv
                DST2.CS[6][N] = DST2.CS[6][N] + DST2.PV[2][M] * DST2.PV[2][M]; //sample sum of ww
                DST2.CS[7][N] = DST2.CS[7][N] + DST2.PV[0][M] * DST2.PV[1][M]; //sample sum of uv
                DST2.CS[8][N] = DST2.CS[8][N] + DST2.PV[0][M] * DST2.PV[2][M]; //sample sum of uw
                DST2.CS[9][N] = DST2.CS[9][N] + DST2.PV[1][M] * DST2.PV[2][M]; //sample sum of vw
                //third order
                DST2.CS[10][N] = DST2.CS[10][N] + DST2.PV[0][M] * DST2.PV[0][M] * DST2.PV[0][M]; //sample sum of uuu
                DST2.CS[11][N] = DST2.CS[11][N] + DST2.PV[0][M] * DST2.PV[0][M] * DST2.PV[1][M]; //sample sum of uuv
                DST2.CS[12][N] = DST2.CS[12][N] + DST2.PV[0][M] * DST2.PV[0][M] * DST2.PV[2][M]; //sample sum of uuw

                DST2.CS[13][N] = DST2.CS[13][N] + DST2.PV[1][M] * DST2.PV[1][M] * DST2.PV[0][M]; //sample sum of vvu
                DST2.CS[14][N] = DST2.CS[14][N] + DST2.PV[1][M] * DST2.PV[1][M] * DST2.PV[1][M]; //sample sum of vvv
                DST2.CS[15][N] = DST2.CS[15][N] + DST2.PV[1][M] * DST2.PV[1][M] * DST2.PV[2][M]; //sample sum of vvw

                DST2.CS[16][N] = DST2.CS[16][N] + DST2.PV[2][M] * DST2.PV[2][M] * DST2.PV[0][M]; //sample sum of wwu
                DST2.CS[17][N] = DST2.CS[17][N] + DST2.PV[2][M] * DST2.PV[2][M] * DST2.PV[1][M]; //sample sum of wwv
                DST2.CS[18][N] = DST2.CS[18][N] + DST2.PV[2][M] * DST2.PV[2][M] * DST2.PV[2][M]; //sample sum of www

            } //partice loop
        }     //if
    }         //MNC loop
}

///////INDEX/////////////
void INDEX()
{
    int NCELL, M, k;
    //printf("\n INDEX particles according to IR[k]\n");

    for (int j = 0; j < DST1.MNC; j++)
    {
        DST2.IC[1][j] = 0;
    }
    for (int j = 0; j < DST1.NM; j++)
    {
        NCELL = DST2.PP[2][j];
        DST2.IC[1][NCELL] = DST2.IC[1][NCELL] + 1;
    }

    M = 0;
    k = 0;
    for (int j = 0; j < DST1.MNC; j++)
    {
        DST2.IC[0][j] = M;
        M = M + DST2.IC[1][j];
    }

    for (int j = 0; j < DST1.MNC; j++)
    {
        DST2.IC[1][j] = 0;
    }

    for (int j = 0; j < DST1.NM; j++)
    {
        NCELL = DST2.PP[2][j];
        DST2.IC[1][NCELL] = DST2.IC[1][NCELL] + 1;
        k = DST2.IC[0][NCELL] + DST2.IC[1][NCELL];
        k = k - 1;
        DST2.IR[k] = j;
        //particle number j:(0)--(NM-1) has been set in the cross reference array
        //k based on counting cellBYcell referencing, starting in array from 0
    }
}

////////Move Function//////
void Move()
{
    int IFT, N, MC, Ni, Nj, k, kks, NCELL, ks, IFSurf, dir;
    long double AT, X, Y, XI, YI, XC, YC, DX, DY, YS, XS;
    int GT150, GT100;
    dir = 0;
    DST1.REPTN = DST1.REPTN + 1;

    if (DST1.REPTN <= 9000)
    {
        //TRAckParticle();
        TRAckParticle_Groups();
    }

    ////
    IFT = -1;
    ks = 0;
    GT150 = -1;
    GT100 = -1;
    //A negative IFT means that no particles have entered at this step
    N = -1;
    while (GT100 == -1)
    {
        GT150 = -1; //GT100=-1;
        //Move_100: N=N+1;
        N = N + 1;
        if (N < DST1.NM)
        {
            //printf("\n particle %d is selected for move\n", N);

            if (IFT < 0)
                AT = DST1.DTM;
            if (IFT > 0)
                AT = DST1.DTM * RF(Irandom);
            while (GT150 == -1)
            {
                GT150 = -1;
                //Move_150: DST1.MOVT=DST1.MOVT+1;
                DST1.MOVT = DST1.MOVT + 1;
                MC = DST2.PP[2][N]; //initial cell number
                XI = DST2.PP[0][N];
                YI = DST2.PP[1][N];
                if ((XI + 0.00001 * DST2.CG[2][0][0] < 0.0) || (XI - 0.00001 * DST2.CG[2][DST1.NX - 1][0] > DST2.CB[0]))
                {
                    printf("remove x\n");
                    getchar();
                    //printf("\n particle %d is outside X Coordinate: %3.10Lf\n", N,XI);
                    Remove(&N);
                    GT150 = 100;
                    break;
                    //goto Move_100;
                }

                if ((YI + 0.00001 * DST2.CG[5][0][0] < 0.0) || (YI - 0.00001 * DST2.CG[5][0][DST1.NY - 1] > DST2.CB[1]))
                {
                    //printf("\n particle %d is outside Y Coordinate: %3.10Lf\n", N,YI);
                    printf("remove n %d @y: %Lf iter %ld\n", N, DST2.PP[1][N], DST1.REPTN);
                    getchar();
                    Remove(&N);
                    GT150 = 100;
                    break;
                    //goto Move_100;
                }
                DX = DST2.PV[0][N] * AT;
                DY = DST2.PV[1][N] * AT;
                if (fabsl(GravX) > 1e-7)
                    DX = DX + 0.5 * GravX * AT * AT;
                if (fabsl(GravY) > 1e-7)
                    DY = DY + 0.5 * GravY * AT * AT;
                X = XI + DX;
                Y = YI + DY;
                for (ks = 0; ks < 4 + DST1.IFSurf; ks++)
                {
                    // printf("ks %d ibks %d \n", ks, DST2.IB[ks]);
                    // getchar();
                    if ((DST2.IB[ks] == 2) || (DST2.IB[ks] == 3)) //IB=1 Stream; IB=2 Surface; IB=3 Periodic; 10 Vertical surface middle
                    {
                        //printf("\n Side %d is surface: %d\n", ks,DST2.IB[ks]);
                        ///////side 0-1
                        if ((ks == 0) || (ks == 1))
                        {
                            if (DST2.IB[ks] == 2)
                            { //SURFACE X
                                if (ks == 0)
                                    XS = 0.0;
                                if (ks == 1)
                                    XS = DST2.CB[0];
                                if (((ks == 0) && ((XI > XS) && (X < XS))) || ((ks == 1) && ((XI < XS) && (X > XS))))
                                {
                                    //printf("\n Particle %d for surface collision @ XI/XS: %Lf X/XS: %Lf at surface %d\n",N,XI/XS,X/XS,ks);
                                    YC = YI + (XS - XI) * DY / DX;

                                    if ((YC < 0) && (DST2.IB[2] == 2))
                                    {
                                        YC = fmod(YC, DST2.CB[1]);
                                        YC = -YC;
                                        DST2.PV[1][N] = -DST2.PV[1][N];
                                    } //collision place symmetry inside
                                    if ((YC > DST2.CB[1]) && (DST2.IB[3] == 2))
                                    {
                                        YC = fmod(YC, DST2.CB[1]);
                                        YC = DST2.CB[1] - YC;
                                        DST2.PV[1][N] = -DST2.PV[1][N];
                                    } //collision place symmetry inside

                                    if ((YC < 0) && (DST2.IB[2] == 3))
                                    {
                                        YC = fmod(YC, DST2.CB[1]);
                                        YC = YC + DST2.CB[1];
                                    } //consider periodic plane Y=0
                                    if ((YC > DST2.CB[1]) && (DST2.IB[3] == 3))
                                    {
                                        YC = fmod(YC, DST2.CB[1]);
                                    } //consider periodic plane Y=cb1
                                    //calculate cell number of surface collision//
                                    Nj = (int)(YC / DST1.CH);
                                    if (Nj < 0)
                                        Nj = 0;
                                    if (Nj > DST1.NY - 1)
                                        Nj = DST1.NY - 1;
                                    if (ks == 0)
                                    {
                                        MC = Nj * DST1.NX;
                                        MC = DST2.IJtC[0][Nj];
                                    }

                                    if (ks == 1)
                                    {
                                        MC = Nj * DST1.NX + (DST1.NX - 1);
                                        MC = DST2.IJtC[DST1.NX - 1][Nj];
                                    }

                                    AT = AT * (X - XS) / DX; //AT: remained travel-time after particle hits a surface
                                    if (AT < 0)
                                    {
                                        printf("Warning neg time");
                                        getchar();
                                    }
                                    REFLECT2D(N, ks, XS, YC, MC, dir);
                                    GT150 = -1;
                                    break; //reach to loop while GT150
                                    //goto Move_150;
                                }
                            }

                            else if (DST2.IB[ks] == 3)
                            { //periodic boundary X
                                //////-----------

                                if ((DST2.IB[2] == 2) || (DST2.IB[3] == 2))
                                {

                                    if (((YI > 0.) && (Y < 0.)) && (DST2.IB[2] == 2))
                                    {

                                        YS = 0.0;
                                        XC = XI + (YS - YI) * DX / DY;
                                        XC = fmod(XC, DST2.CB[0]);
                                        if (XC < 0.)
                                            XC = XC + DST2.CB[0]; //periodic X=0
                                        if (XC > DST2.CB[0])
                                            XC = XC - DST2.CB[0]; //periodic X=cb0
                                        Ni = (int)(XC / DST1.CW);
                                        Nj = 0;
                                        MC = DST2.IJtC[Ni][0];
                                        AT = AT * (Y - YS) / DY; //AT: remained travel-time after particle hits a surface
                                        if (AT < 0)
                                        {
                                            printf("Warning neg time");
                                            getchar();
                                        }
                                        DST2.PP[2][N] = MC;
                                        DST2.PP[0][N] = XC;
                                        DST2.PP[1][N] = YS;
                                        REFLECT2D(N, 2, XC, YS, MC, dir);
                                        GT150 = -1;
                                        break; //reach to loop while GT150
                                               //goto Move_150;
                                    }
                                    //

                                    if (((YI < DST2.CB[1]) && (Y > DST2.CB[1])) && (DST2.IB[3] == 2))
                                    {

                                        YS = DST2.CB[1];
                                        XC = XI + (YS - YI) * DX / DY;
                                        XC = fmod(XC, DST2.CB[0]);
                                        if (XC < 0.)
                                            XC = XC + DST2.CB[0]; //periodic X=0
                                        if (XC > DST2.CB[0])
                                            XC = XC - DST2.CB[0]; //periodic X=cb0
                                        Ni = (int)(XC / DST1.CW);
                                        Nj = DST1.NY - 1;
                                        MC = DST2.IJtC[Ni][Nj];
                                        AT = AT * (Y - YS) / DY; //AT: remained travel-time after particle hits a surface
                                        if (AT < 0)
                                        {
                                            printf("Warning neg time");
                                            getchar();
                                        }
                                        DST2.PP[2][N] = MC;
                                        DST2.PP[0][N] = XC;
                                        DST2.PP[1][N] = YS;
                                        REFLECT2D(N, 3, XC, YS, MC, dir);
                                        GT150 = -1;
                                        break; //reach to loop while GT150
                                               //goto Move_150;
                                    }
                                }
                                X = fmod(X, DST2.CB[0]);
                                if (X < 0.)
                                    X = X + DST2.CB[0]; //periodic X=0
                                if (X > DST2.CB[0])
                                    X = X - DST2.CB[0]; //periodic X=cb0

                                Ni = (int)(X / DST1.CW);
                                Nj = (int)(Y / DST1.CH);
                                if (Ni < 0)
                                    Ni = 0;
                                if (Ni > DST1.NX - 1)
                                    Ni = DST1.NX - 1;
                                if (Nj < 0)
                                    Nj = 0;
                                if (Nj > DST1.NY - 1)
                                    Nj = DST1.NY - 1;
                                MC = DST2.IJtC[Ni][Nj];
                                DST2.PP[2][N] = MC;
                                DST2.PP[0][N] = X;
                                DST2.PP[1][N] = Y;

                                AT = 0.0; //AT: remained travel-time after particle hits a surface
                                GT150 = 100;
                                continue;
                                //////------
                            }
                        } //    if((ks==0)||(ks==1)) {
                        ///////side 2-3
                        if ((ks == 2) || (ks == 3))
                        {
                            // printf("ks 2 or 3 : %d IBks %d \n", ks, DST2.IB[ks]);
                            // getchar();
                            if (DST2.IB[ks] == 2)
                            { //SURFACE Y
                                if (ks == 2)
                                    YS = 0.0;
                                if (ks == 3)
                                    YS = DST2.CB[1];
                                if (((ks == 2) && ((YI > YS) && (Y < YS))) || ((ks == 3) && ((YI < YS) && (Y > YS))))
                                {

                                    //printf("\n One particle for surface collision %d YI/YS: %Lf Y/YS: %Lf at surface %d\n",N,YI/YS,Y/YS,ks);
                                    XC = XI + (YS - YI) * DX / DY;

                                    if ((XC < 0) && (DST2.IB[0] == 2))
                                    {
                                        XC = fmod(XC, DST2.CB[0]);
                                        XC = -XC;
                                        DST2.PV[0][N] = -DST2.PV[0][N];
                                    } //collision place symmetry inside

                                    if ((XC > DST2.CB[0]) && (DST2.IB[1] == 2))
                                    {
                                        XC = fmod(XC, DST2.CB[0]);
                                        XC = DST2.CB[0] - XC;
                                        DST2.PV[0][N] = -DST2.PV[0][N];
                                    } //collision place symmetry inside

                                    if ((XC < 0) && (DST2.IB[0] == 3))
                                    {
                                        XC = fmod(XC, DST2.CB[0]);
                                        XC = XC + DST2.CB[0];
                                    } //consider periodic plane X=0
                                    if ((XC > DST2.CB[0]) && (DST2.IB[1] == 3))
                                    {
                                        XC = fmod(XC, DST2.CB[0]);
                                    } //consider periodic plane X=cb0
                                    //calculate cell number of surface collision//
                                    Ni = (int)(XC / DST1.CW);
                                    if (Ni < 0)
                                        Ni = 0;
                                    if (Ni > DST1.NX - 1)
                                        Ni = DST1.NX - 1;
                                    if (ks == 2)
                                    {
                                        MC = Ni;
                                        MC = DST2.IJtC[Ni][0];
                                    }
                                    if (ks == 3)
                                    {
                                        MC = Ni + (DST1.NY - 1) * DST1.NX;
                                        MC = DST2.IJtC[Ni][DST1.NY - 1];
                                    }
                                    AT = AT * (Y - YS) / DY; //AT: remained travel-time after particle hits a surface
                                    if (AT < 0)
                                    {
                                        printf("Warning neg time");
                                        getchar();
                                    }
                                    REFLECT2D(N, ks, XC, YS, MC, dir);
                                    GT150 = -1;
                                    break; //reach to loop while GT150
                                    //goto Move_150;
                                }
                            }
                            else if (DST2.IB[ks] == 3)
                            { //periodic boundary Y
                                ///---------------
                                if ((DST2.IB[0] == 2) || (DST2.IB[1] == 2))
                                {

                                    if (((XI > 0.) && (X < 0.)) && (DST2.IB[0] == 2))
                                    {
                                        XS = 0.0;
                                        YC = YI + (XS - XI) * DY / DX;
                                        YC = fmod(YC, DST2.CB[1]);
                                        if (YC < 0.)
                                            YC = YC + DST2.CB[1]; //periodic Y=0
                                        if (YC > DST2.CB[1])
                                            YC = YC - DST2.CB[1]; //periodic Y=cb1
                                        Ni = 0;
                                        Nj = (int)(YC / DST1.CH);
                                        MC = DST2.IJtC[Ni][Nj];
                                        AT = AT * (X - XS) / DX; //AT: remained travel-time after particle hits a surface
                                        if (AT < 0)
                                        {
                                            printf("Warning neg time");
                                            getchar();
                                        }
                                        DST2.PP[2][N] = MC;
                                        DST2.PP[0][N] = XS;
                                        DST2.PP[1][N] = YC;
                                        REFLECT2D(N, 0, XS, YC, MC, dir);
                                        GT150 = -1;
                                        break; //reach to loop while GT150
                                        //goto Move_150;
                                    }
                                    //

                                    if (((XI < DST2.CB[0]) && (X > DST2.CB[0])) && (DST2.IB[1] == 2))
                                    {
                                        XS = DST2.CB[0];
                                        YC = YI + (XS - XI) * DY / DX;
                                        YC = fmod(YC, DST2.CB[1]);
                                        if (YC < 0.)
                                            YC = YC + DST2.CB[1]; //periodic Y=0
                                        if (YC > DST2.CB[1])
                                            YC = YC - DST2.CB[1]; //periodic Y=cb1
                                        Ni = DST1.NX - 1;
                                        Nj = (int)(YC / DST1.CH);
                                        MC = DST2.IJtC[Ni][Nj];
                                        AT = AT * (X - XS) / DX; //AT: remained travel-time after particle hits a surface
                                        if (AT < 0)
                                        {
                                            printf("Warning neg time");
                                            getchar();
                                        }
                                        DST2.PP[2][N] = MC;
                                        DST2.PP[0][N] = XS;
                                        DST2.PP[1][N] = YC;
                                        REFLECT2D(N, 1, XS, YC, MC, dir);
                                        GT150 = -1;
                                        break; //reach to loop while GT150
                                        //goto Move_150;
                                    }
                                    //
                                }
                                Y = fmod(Y, DST2.CB[1]);
                                if (Y < 0.)
                                    Y = Y + DST2.CB[1]; //periodic Y=0
                                if (Y > DST2.CB[1])
                                    Y = Y - DST2.CB[1]; //periodic Y=cb1

                                ///////////////////////////////////////////
                                if (Y < 0.0 || Y > DST2.CB[1])
                                {
                                    printf("Y %Lf %d  %d %Lf", Y, DST2.IB[ks], ks, DST2.CB[1]);
                                    getchar();
                                }
                                ////////////////////////////////////////////

                                Ni = (int)(X / DST1.CW);
                                Nj = (int)(Y / DST1.CH);
                                if (Ni < 0)
                                    Ni = 0;
                                if (Ni > DST1.NX - 1)
                                    Ni = DST1.NX - 1;
                                if (Nj < 0)
                                    Nj = 0;
                                if (Nj > DST1.NY - 1)
                                    Nj = DST1.NY - 1;
                                MC = DST2.IJtC[Ni][Nj];
                                DST2.PP[2][N] = MC;
                                DST2.PP[0][N] = X;
                                DST2.PP[1][N] = Y;
                                AT = 0.0; //AT: remained travel-time after particle hits a surface
                                GT150 = 100;
                                break;
                                //////------

                                ////---------------------
                            }
                        } //if((ks==2)||(ks==3))
                    }     //if ((DST2.IB[ks]==2)||(DST2.IB[ks]==3))

                    if (DST2.IB[ks] == 10) //IB=1 Stream; IB=2 Surface; IB=3 Periodic; IB=10 Vertical surface middle
                    {                      //DST1.SurfHeightL; DST1.SurfDistL
                        if (((XI <= DST1.SurfDistL) && (X >= DST1.SurfDistL)) || ((XI >= DST1.SurfDistL) && (X <= DST1.SurfDistL)))
                        { //the x movement crosses surfaces
                            if (((YI <= DST1.SurfHeightL) && (YI >= 0.0)) || ((Y <= DST1.SurfHeightL) && (Y >= 0.0)))
                            {
                                YS = (Y - YI) / (X - XI) * DST1.SurfDistL - (Y - YI) / (X - XI) * XI + YI;
                                if (((YS <= DST1.SurfHeightL) && (YS >= 0.0)))
                                { //if the particle hits the surface at (DST1.SurfDistL,YS)
                                    XC = DST1.SurfDistL;
                                    AT = AT * fabsl((X - XC) / DX); //AT: remained travel-time after particle hits a surface
                                    if (AT < 0)
                                    {
                                        printf("Warning neg time");
                                        getchar();
                                    }
                                    Ni = (int)(XC / DST1.CW);
                                    Nj = (int)(YS / DST1.CH);
                                    if (Ni < 0)
                                        Ni = 0;
                                    if (Ni > DST1.NX - 1)
                                        Ni = DST1.NX - 1;
                                    if (Nj < 0)
                                        Nj = 0;
                                    if (Nj > DST1.NY - 1)
                                        Nj = DST1.NY - 1;
                                    MC = DST2.IJtC[Ni][Nj];
                                    //determine direction of reflection dir=1 pv[0]=>-, dir=2 pv[0]=>+
                                    if (DX > 0)
                                    {
                                        dir = 1;
                                    } //left to right movement, reflection right to left: pv[0]
                                    else if (DX < 0)
                                    {
                                        dir = 2;
                                    } //right to left movement, reflection left to right
                                    else if ((DX == 0) && (XI < XC))
                                    {
                                        dir = 1;
                                    } //left to right movement, reflection right to left
                                    else if ((DX == 0) && (XI > XC))
                                    {
                                        dir = 2;
                                    } //right to left movement, reflection left to right
                                    else if ((DX == 0) && (XI == XC))
                                    {
                                        dir = (int)(RF(Irandom) * 2 + 1);
                                    } //right to left movement, reflection left to right
                                    if (dir == 3)
                                    {
                                        printf("Warning strange dir");
                                        getchar();
                                    }

                                    REFLECT2D(N, ks, XC, YS, MC, dir);
                                    GT150 = -1;
                                    break; //reach to loop while GT150
                                    //goto Move_150;

                                } //if the particle hits the surface at (DST1.SurfDistL,YS)
                            }     //the y movement in the range
                        }         //the x movement crosses surfaces
                    }             //10 Vertical surface middle

                    //
                } //for ks=0-3
                if (ks >= 4 + DST1.IFSurf)
                    break;
            } //..////////////////////////while(GT150==-1)
            if (GT150 == 100)
                continue; //if the GT150 loop is "continued" before (i.e. GT150>0) then another continue for GT100 loop is desired
            //after checking 4 sides, if the particle is outside domain should be removed in case of openboundary
            //
            DST2.PP[0][N] = X;
            DST2.PP[1][N] = Y;

            if ((X < 0) || (X > DST2.CB[0]))
            {
                if (X < 0)
                    k = 0;
                if (X > DST2.CB[0])
                    k = 1;

                if (DST2.IB[k] == 3)
                { //periodic boundary X
                    X = fmod(X, DST2.CB[0]);
                    if (X < 0)
                        X = X + DST2.CB[0];
                    if (X > DST2.CB[0])
                        X = X - DST2.CB[0];
                    DST2.PP[0][N] = X;
                }
                else if (DST2.IB[k] == 2)
                {
                    X = fmod(X, DST2.CB[0]);
                    if (k == 0)
                        X = -X;
                    if (k == 1)
                        X = DST2.CB[0] - X;
                    DST2.PV[0][N] = -DST2.PV[0][N];
                    DST2.PP[0][N] = X;
                }
                else
                { //particle leaving the flow
                    Remove(&N);
                    GT100 = -1;
                    continue;
                    //goto Move_100;
                }
            }
            //

            if ((Y < 0) || (Y > DST2.CB[1]))
            {
                if (Y < 0)
                    k = 2;
                if (Y > DST2.CB[1])
                    k = 3;

                if (DST2.IB[k] == 3)
                { //periodic boundary Y
                    Y = fmod(Y, DST2.CB[1]);
                    if (Y < 0)
                        Y = Y + DST2.CB[1];
                    if (Y > DST2.CB[1])
                        Y = Y - DST2.CB[1];
                    DST2.PP[1][N] = Y;
                }
                else if (DST2.IB[k] == 2)
                {
                    Y = fmod(Y, DST2.CB[1]);
                    if (k == 2)
                        Y = -Y;
                    if (k == 3)
                        Y = DST2.CB[1] - Y;
                    DST2.PV[1][N] = -DST2.PV[1][N];
                    DST2.PP[1][N] = Y;
                }
                else
                { //particle leaving the flow
                    Remove(&N);
                    GT100 = -1;
                    continue;
                    //goto Move_100;
                }
            }
            //Particle move finished
            Ni = (int)(DST2.PP[0][N] / DST1.CW);
            Nj = (int)(DST2.PP[1][N] / DST1.CH);
            NCELL = DST2.IJtC[Ni][Nj];
            DST2.PP[2][N] = NCELL;
            if (AT < 0)
            {
                printf("\n Time step is negative!! %Lf\n", AT);
                getchar();
            }

            //GT100=-1;
            continue;
            //goto Move_100;
        } //if (N<DST1.NM)
        else
        { //else if N>=DST1.NM

            if (IFT < 0)
            {
                IFT = 1; //Enter new particles from open-boundaries
                ENTER();
                N = N - 1;
                //GT100=-1;//it needs to be corrected for STREAMING!!!!
                //goto Move_100;
            } //end else if (IFT<0)
            else if (IFT > 0)
            {
                break;
            }

        } //else { //else if N>=DST1.NM

    } //while(GT100==-1)
    if ((fabsl(GravX) > 1e-7) || (fabsl(GravY) > 1e-7))
    {
        AT = DST1.DTM;
        N = 0;
        for (N = 0; N < DST1.NM; N++)
        {
            //printf("\n N is %d\n",N);

            if (fabsl(GravX) > 1e-7)
            {
                DST2.PV[0][N] = DST2.PV[0][N] + GravX * AT;
            }
            if (fabsl(GravY) > 1e-7)
            {
                DST2.PV[1][N] = DST2.PV[1][N] + GravY * AT;
            }
        }
    }

} //end Move
////////

//collision
void Collision()
{
    int n, NSEL, i, j, ISEL, L, M;
    long double ASEL, AVN, SN, CVM, A, VR;
    //printf("\n Select particles for collision and change of velocities\n");

    for (n = 0; n < DST1.MNC; n++)
    {
        i = DST2.CtIJ[0][n]; //row
        j = DST2.CtIJ[1][n]; //column
        SN = DST2.CS[0][n];  //cell number of sampled done (will use it in AVN calculation)

        if (SN > 1)
        {
            AVN = SN / (float)DST1.NSMP;
        }
        else
        {
            AVN = DST2.IC[1][n];
        }
        //AVN is the average number of particles in the cell
        ASEL = 0.5 * DST2.IC[1][n] * AVN * DST1.FNUM * DST2.CCG[0][i][j] * DST1.DTM / DST2.CC[i][j] +
               DST2.CCG[1][i][j];
        //ASEL is the number of pairs to be selected, see eqn (11.5)
        NSEL = ASEL;
        DST2.CCG[1][i][j] = ASEL - NSEL;

        if (NSEL > 0)
        { //if there is any pair to select
            if (DST2.IC[1][n] < 2)
            { //if there is not enough particle for collision
                //then the number NSEL is added to the remainder of CCG1
                DST2.CCG[1][i][j] = DST2.CCG[1][i][j] + NSEL;
            }
            else
            { //if there are enough particles for collision, now start pair selection
                CVM = DST2.CCG[0][i][j];
                DST1.SELT = DST1.SELT + NSEL;
                for (ISEL = 1; ISEL <= NSEL; ISEL++)
                {
                    SELECT(n, i, j, &L, &M, &VR);
                    //printf("\n  %d %d %d %d %d\n", n,i,j,L,M);
                    if (DST1.CVR > CVM)
                        CVM = DST1.CVR; //update the maximum CVM

                    if (DST1.CVR / DST2.CCG[0][i][j] > RF(Irandom))
                    { //check the acceptance of a collision

                        DST1.NCOL = DST1.NCOL + 1;
                        DST2.ICol[0][n] = DST2.ICol[0][n] + 1; //collision counter of cell
                        A = pow(DST2.PP[0][L] - DST2.PP[0][M], 2.) + pow(DST2.PP[1][L] - DST2.PP[1][M], 2.);
                        A = sqrtl(A);
                        DST1.SEPT = DST1.SEPT + A;
                        DST2.DCol[0][n] = DST2.DCol[0][n] + A; //separation counter of cell
                        ELASTIC(L, M, VR);
                    }
                }
                DST2.CCG[0][i][j] = CVM;
            } //if there are enough particles for collision, now start pair selection
        }     //if there is any pair to select

    } // for all the cells
} //end collision
  ///////
  //////
////GAMA FUNCTION////
long double GAM(long double X)
{
    long double gama, A, Y, powY2, powY3, powY4, powY5;
    A = 1.0;
    Y = X;

    if (Y < 1.0)
    {
        A = A / Y;
    }
    else
    {
    GAM_50:
        Y = Y - 1;
        if (Y > 1.0)
        {
            A = A * Y;
            goto GAM_50;
        }
    }

    powY2 = pow(Y, 2.0);
    powY3 = pow(Y, 3.0);
    powY4 = pow(Y, 4.0);
    powY5 = pow(Y, 5.0);
    gama = A * (1. - 0.5748646 * Y + 0.9512363 * powY2 - 0.6998588 * powY3 + 0.4245549 * powY4 - 0.1010678 * powY5);
    return gama;
}
//

////Random FUNCTION////long double X
long double RF(int X)
{
    long double r = 0.0;
    //long double s = 1.0;
    ///do {
    ///s /= RAND_MAX + 1.0;
    ///r += rand() * s;
    ///} while (s > DBL_EPSILON);
    if (X == 1)
    {
        r = drand48();
    } // generate uniformly distributed pseudo-random numbers using a linear congruential algorithm and 48-bit integer arithmetic.
    if (X == 2)
    {
        r = genrand64_real2();
    } //Mersenne Twister [0,1) //http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/C-LANG/mt19937-64.c
    return r;
}
//
///////RVELC Function///////
//Generate two random velocity components U and V in an equilibrium
//gas with most probable speed VMP  (based on eqns (C10) and (C12))
void RVELC(long double *U, long double *V, long double VMP)

{
    long double A, B;
    A = sqrtl(-logl(RF(Irandom)));
    B = 6.283185308 * RF(Irandom);
    *U = A * sinl(B) * VMP;
    *V = A * cosl(B) * VMP;
}
////ERF FUNCTION////long double X
long double ERF(long double S)
{
    long double erf, B, D, C, T;
    B = fabs(S);

    if (B > 4.)
    {
        D = 1.0;
    }

    else
    {
        C = exp(-B * B);
        T = 1. / (1. + 0.3275911 * B);
        D = 1. - (0.254829592 * T - 0.284496736 * T * T + 1.421413741 * T * T * T -
                  1.453152027 * T * T * T * T + 1.061405429 * T * T * T * T * T) *
                     C;
    }

    if (S < 0.)
    {
        D = -D;
    }
    erf = D;
    return erf;
}
/////////

void Remove(int *N)
{

    //printf("\n particle %d is removed here\n", N);
    DST2.PP[0][*N] = DST2.PP[0][DST1.NM - 1];
    DST2.PP[1][*N] = DST2.PP[1][DST1.NM - 1];
    DST2.PP[2][*N] = DST2.PP[2][DST1.NM - 1];
    //
    DST2.PV[0][*N] = DST2.PV[0][DST1.NM - 1];
    DST2.PV[1][*N] = DST2.PV[1][DST1.NM - 1];
    DST2.PV[2][*N] = DST2.PV[2][DST1.NM - 1];
    //set end out
    DST2.PP[0][DST1.NM - 1] = -9000000;
    DST2.PP[1][DST1.NM - 1] = -9000000;
    DST2.PP[2][DST1.NM - 1] = -9000000;
    //
    DST2.PV[0][DST1.NM - 1] = -9000000;
    DST2.PV[1][DST1.NM - 1] = -9000000;
    DST2.PV[2][DST1.NM - 1] = -9000000;
    //
    DST1.NM = DST1.NM - 1;
    *N = *N - 1;
}
//

void REFLECT2D(int n, int k, long double X1, long double Y1, int m, int dir)
{
    long double VMP, a1, b1, c1;
    int Specular;
    Specular = -10;
    //printf("\n input %d %d %Lf %Lf %d\n", n,k,X1,Y,m);
    //n: the number of particle N=0,NM-1
    //k: number of surface k=0,1,2,3  x=0,x=cb0,y=0,y=cb1
    //X1: X location of surface hitting point
    //Y: Y location of surface hitting point
    //m: Cell number of hitiing point
    ///////////////////////////////////////
    if (k == 3)
    {
        DST2.PV[0][n] = DST2.PV[0][n] - WallvelX;
    }
    //if the upper wall is selected, check wall movement X direction
    //-WallvelX is added to PV0 if surface sampling is desired
    //Surface sampling of upper wall should be collected here
    ///////////////////////////////////////
    //check whether the bottom specular relection required
    if ((k == 2) && (SpecDist > 0))
    {
        if (Y1 == 0)
        {
            if ((X1 >= 0) && (X1 < SpecDist * DST2.CB[0])) //specular reflection is needed [0,SpecDist)
                Specular = 100;
        }
    }
    if ((DST2.TSURF[k] < 0) || (Specular > 0))
    { //TSURF>0 diffuse, TSURF<0 specular

        //printf("\n Specular reflection\n");
        if ((k == 2) || (k == 3))
        {
            DST2.PV[1][n] = -DST2.PV[1][n];
        }
        if ((k == 0) || (k == 1))
        {
            DST2.PV[0][n] = -DST2.PV[0][n];
        }
        if (k == 4)
        {
            if (dir == 1)
            {
                DST2.PV[0][n] = -DST2.PV[0][n];
            }
            if (dir == 2)
            {
                DST2.PV[0][n] = DST2.PV[0][n];
            }
        }
    }

    else
    { //Diffuse reflection
        // printf("\n Diffuse reflection\n");
        //printf("Boltz%Le,SP4 mass: %Le \n",BOLTZ,DST2.sp[4]);
        VMP = sqrtl(2. * BOLTZ * DST2.TSURF[k] / DST2.sp[4]);
        if (k == 0)
        {
            DST2.PV[0][n] = sqrtl(-logl(RF(Irandom))) * VMP;
            RVELC(&DST2.PV[1][n], &DST2.PV[2][n], VMP);
        }
        if (k == 1)
        {
            DST2.PV[0][n] = -sqrtl(-logl(RF(Irandom))) * VMP;
            RVELC(&DST2.PV[1][n], &DST2.PV[2][n], VMP);
        }
        if (k == 2)
        {
            DST2.PV[1][n] = sqrtl(-logl(RF(Irandom))) * VMP;
            RVELC(&DST2.PV[0][n], &DST2.PV[2][n], VMP);
        }
        if (k == 3)
        {
            DST2.PV[1][n] = -sqrtl(-logl(RF(Irandom))) * VMP;
            RVELC(&DST2.PV[0][n], &DST2.PV[2][n], VMP);
        }

        if (k == 4)
        {
            if (dir == 1)
            {
                DST2.PV[0][n] = -sqrtl(-logl(RF(Irandom))) * VMP;
                RVELC(&DST2.PV[1][n], &DST2.PV[2][n], VMP);
            }
            if (dir == 2)
            {
                DST2.PV[0][n] = sqrtl(-logl(RF(Irandom))) * VMP;
                RVELC(&DST2.PV[1][n], &DST2.PV[2][n], VMP);
            }
        }
    }
    //new particle positions at reflection point
    if (k == 0)
    {
        DST2.PP[0][n] = X1 + 0.0000001 * DST1.CW;
        DST2.PP[1][n] = Y1;
    }
    if (k == 1)
    {
        DST2.PP[0][n] = X1 - 0.0000001 * DST1.CW;
        DST2.PP[1][n] = Y1;
    }
    if (k == 2)
    {
        DST2.PP[0][n] = X1;
        DST2.PP[1][n] = Y1 + 0.0000001 * DST1.CH;
    }
    if (k == 3)
    {
        DST2.PP[0][n] = X1;
        DST2.PP[1][n] = Y1 - 0.0000001 * DST1.CH;
    }

    if ((k == 4) && ((dir == 1) || (dir == 2)))
    {
        if (dir == 1)
        {
            DST2.PP[0][n] = X1 - 0.0000001 * DST1.CW;
            DST2.PP[1][n] = Y1;
        }
        if (dir == 2)
        {
            DST2.PP[0][n] = X1 + 0.0000001 * DST1.CW;
            DST2.PP[1][n] = Y1;
        }
    }

    DST2.PP[2][n] = m;
    // printf("\n cell number: %Lf, m= %d,\n", DST2.PP[2][n],m);
    //+WallvelX is added to PV0 to create wall movement, here at upper wall only
    if (k == 3)
    {
        DST2.PV[0][n] = DST2.PV[0][n] + WallvelX;
    }
}
//

void ENTER()
{
    int N, NC, NCS, M, k, MC, Ni, Nj, NCELL;
    long double A, VMP, SC, FS1, FS2, QA, U, UN;
    for (N = 0; N < 4; N++)
    {
        if (DST2.IB[N] == 1)
        { //it is an open-boundary
            if ((N == 0) || (N == 1))
                NCS = DST1.NY;
            if ((N == 2) || (N == 3))
                NCS = DST1.NX;
            for (NC = 0; NC < NCS; NC++)
            { //for all the column/row spaces
                SC = 0;
                VMP = sqrtl(2. * BOLTZ * FTMP / DST2.sp[4]);
                if ((N == 0) || (N == 1))
                    A = DST2.BME[N] * DST1.CH / DST1.FH + DST2.BMR[N];
                if ((N == 2) || (N == 3))
                    A = DST2.BME[N] * DST1.CW / DST1.FW + DST2.BMR[N];
                M = A;
                DST2.BMR[N] = A - M;
                if (M > 0)
                {
                    if ((N == 0) || (N == 1))
                    {
                        if (fabsl(VFX) > 1.E-6)
                        {
                            if (N == 0)
                                SC = VFX / VMP;
                            if (N == 1)
                                SC = -VFX / VMP;
                        }
                    }

                    if ((N == 2) || (N == 3))
                    {
                        if (fabsl(VFY) > 1.E-6)
                        {
                            if (N == 2)
                                SC = VFY / VMP;
                            if (N == 3)
                                SC = -VFY / VMP;
                        }
                    }
                    //
                    FS1 = SC + sqrtl(SC * SC + 2.);
                    FS2 = 0.5 * (1. + SC * (2. * SC - FS1));
                    //the above constants are required for the entering distn. of eqn (12.5)
                    for (k = 0; k < M; k++)
                    {
                        if (DST1.NM < DST1.MNM)
                        {
                            DST1.NM = DST1.NM + 1;
                            if (((N < 2) && (fabsl(VFX) > 1.E-6)) || ((N > 1) && (fabsl(VFY) > 1.E-6)))
                            {
                                QA = 3.;
                                if (SC < -3.0)
                                    QA = fabsl(SC) + 1.;
                                UN = -0.01;
                                A = -2;

                                do
                                {
                                    do
                                    {
                                        U = -QA + 2.0 * QA * RF(Irandom); // U is a potential normalised thermal velocity component
                                        UN = U + SC;                      //UN is a potential inward velocity component
                                    } while (UN < 0);

                                    A = (2. * UN / FS1) * exp(FS2 - U * U);
                                } while (A < RF(Irandom));
                                //the inward normalised vel. component has been selected (eq 12-5)
                                if (N == 0)
                                    DST2.PV[0][DST1.NM - 1] = UN * VMP;
                                if (N == 1)
                                    DST2.PV[0][DST1.NM - 1] = -UN * VMP;
                                if (N == 2)
                                    DST2.PV[1][DST1.NM - 1] = UN * VMP;
                                if (N == 3)
                                    DST2.PV[1][DST1.NM - 1] = -UN * VMP;
                            } //if the particles are velocity driven at open-boundaries
                            else
                            { //if velocity does not have any role in driving particles
                                if (N == 0)
                                    DST2.PV[0][DST1.NM - 1] = sqrtl(-logl(RF(Irandom))) * VMP;
                                if (N == 1)
                                    DST2.PV[0][DST1.NM - 1] = -sqrtl(-logl(RF(Irandom))) * VMP;
                                if (N == 2)
                                    DST2.PV[1][DST1.NM - 1] = sqrtl(-logl(RF(Irandom))) * VMP;
                                if (N == 3)
                                    DST2.PV[1][DST1.NM - 1] = -sqrtl(-logl(RF(Irandom))) * VMP;
                            }
                            //get other velocity directions (rather than inflow one)
                            if ((N == 0) || (N == 1))
                            {
                                RVELC(&DST2.PV[1][DST1.NM - 1], &DST2.PV[2][DST1.NM - 1], VMP);
                                DST2.PV[1][DST1.NM - 1] = DST2.PV[1][DST1.NM - 1] + VFY; //tangential velocity
                            }
                            if ((N == 2) || (N == 3))
                            {
                                RVELC(&DST2.PV[0][DST1.NM - 1], &DST2.PV[2][DST1.NM - 1], VMP);
                                DST2.PV[0][DST1.NM - 1] = DST2.PV[0][DST1.NM - 1] + VFX; //tangential velocity
                            }
                            //set positions of entering particles
                            //perpendicular to boundary positions
                            if (N == 0)
                                DST2.PP[0][DST1.NM - 1] = 0.0 + 0.0001 * DST1.CW;
                            if (N == 1)
                                DST2.PP[0][DST1.NM - 1] = DST2.CB[0] - 0.0001 * DST1.CW;
                            if (N == 2)
                                DST2.PP[1][DST1.NM - 1] = 0.0 + 0.0001 * DST1.CH;
                            if (N == 3)
                                DST2.PP[1][DST1.NM - 1] = DST2.CB[1] - 0.0001 * DST1.CH;
                            //Tangential to boundary positions
                            if ((N == 0) || (N == 1))
                            {
                                //if (N==0) MC=NC*DST1.NX;
                                //if (N==1) MC=NC*DST1.NX+(DST1.NX-1);
                                //DST2.PP[2][DST1.NM]=MC;
                                DST2.PP[1][DST1.NM - 1] = DST2.CG[3][0][NC] + RF(Irandom) * DST1.CH;
                            }

                            if ((N == 2) || (N == 3))
                            {
                                //if (N==2) MC=NC;
                                //if (N==3) MC=NC+(DST1.NY-1)*DST1.NX;
                                //DST2.PP[2][DST1.NM]=MC;
                                DST2.PP[0][DST1.NM - 1] = DST2.CG[0][NC][0] + RF(Irandom) * DST1.CW;
                            }
                            Ni = (int)(DST2.PP[0][DST1.NM - 1] / DST1.CW);
                            Nj = (int)(DST2.PP[1][DST1.NM - 1] / DST1.CH);
                            NCELL = DST2.IJtC[Ni][Nj];
                            DST2.PP[2][DST1.NM - 1] = NCELL;

                        } //
                        else
                        {
                            printf("ExceedParticle %d > %d:: increase Max Number of Molecules(MNM)\n", DST1.NM, DST1.MNM);
                            getchar();
                        }

                    } //for each of the entered particles
                }     //if there are any entered particle
            }         //for all the column/row spaces
        }             //if open-boundary
    }                 //for all sides

} //end ENTER
////////SELECT pair Function//////
void SELECT(int N, int i, int j, int *L, int *M, long double *vr)
{
    int k;
    long double VR0, VR1, VR2, VR, VRR;
    *L = -1;
    *M = -2;
    //printf("\n  %d %d %d\n", N,i,j);
    //Select the first and second particle randomly
    k = (int)(RF(Irandom) * (DST2.IC[1][N] - 0.001)) + DST2.IC[0][N] + 1;
    //k a random number between 1 till IC[1]
    k = k - 1;
    *L = DST2.IR[k]; //L is the first randomly selected particle

    do
    {
        k = (int)(RF(Irandom) * (DST2.IC[1][N] - 0.001)) + DST2.IC[0][N] + 1;
        //k a random number between 1 till IC[1]
        k = k - 1;
        *M = DST2.IR[k]; //M is the second randomly selected particle
    } while (*L == *M);
    VR0 = DST2.PV[0][*L] - DST2.PV[0][*M];
    VR1 = DST2.PV[1][*L] - DST2.PV[1][*M];
    VR2 = DST2.PV[2][*L] - DST2.PV[2][*M];
    VR0 = VR0 * VR0;
    VR1 = VR1 * VR1;
    VR2 = VR2 * VR2;
    VRR = VR0 + VR1 + VR2;
    VR = sqrtl(VRR);
    *vr = VR;
    //
    DST1.CVR = ((VR * DST1.SPM0) / (DST1.SPM5)) * pow(((2. * BOLTZ * DST1.SPM1) / (VRR * DST1.SPM4)), (DST1.SPM2 - 0.5));
    //CVR the collision cross-section is based on VR* eqn (4.63)

} //end function

/////////Elastic collision (VHS) of SELECTed particles////
void ELASTIC(int M1, int M2, long double VR)
{
    int i1;
    long double RM1, RM2, VCCM[3], VRCP[3];
    long double A, B, C;

    //VCCM=malloc(3 *sizeof(long double));
    //VCCM defines the components of the centre-of-mass velocity, eqn (2.1)
    //VRCP=malloc(3 *sizeof(long double));
    //VRCP are the post-collision relative velocities

    /*weight of each particle mass: as particles have same type, the weights are equal
but in general: [reduced_mass of pairs]/[the particle mass]
*/
    RM1 = DST1.SPM4 / DST2.sp[4]; //for M1
    RM2 = DST1.SPM4 / DST2.sp[4]; //for M2
    for (i1 = 0; i1 < 3; i1++)
    {
        VCCM[i1] = RM1 * DST2.PV[i1][M1] + RM2 * DST2.PV[i1][M2];
    }
    //VHS logic
    B = 2. * RF(Irandom) - 1.;
    //B is the cosine of a random elevation angle
    A = sqrtl(1. - B * B);
    VRCP[0] = B * VR;
    C = 2. * M_PI * RF(Irandom);
    VRCP[1] = A * cosl(C) * VR;
    VRCP[2] = A * sinl(C) * VR;
    for (i1 = 0; i1 < 3; i1++)
    {
        DST2.PV[i1][M1] = VCCM[i1] + VRCP[i1] * RM1;
        DST2.PV[i1][M2] = VCCM[i1] - VRCP[i1] * RM2;
    }
}

//////////////////
///////
void OUTPUT()
{
    int N, NN, SN, ii, seconds, minute, hour;
    long double A, SM, SM_NN, SMCC, SMCC_B, DENN, DEN, SUM7, UU, UU_NN, TT, Pressure, XC, YC, CF_Theo, CF_Num, CF_Ratio, NP_Avg, TIME_registered, Mach;
    long double MFP_VHS_Th_V, MFP_VHS_RMS_V, MCS, SOF, Kn_thV_chrL, G_L_den, G_L_Vel, Kn_thV_GLvel, Kn_RMS_V_chrL, Kn_RMS_V_GLvel;
    long double MCT, MColT, MTrTX, MTrTY, MColTmin, MTrTXmin, MTrTYmin, SN_NN, A_NN, DEN_NN, Abs_Vel, Abs_Vel_NN;
    long double SMU[3], VEL[3], SMU_NN[3], VEL_NN[3];
    long double Heat_x, Heat_y, P_XX, P_YY, P_XY, P_XZ, P_YZ, SMCC_U, SMCC_V;
    long double X_Vor_C, Y_Vor_C, V_Vor_C, X_Vor_C1, Y_Vor_C1, V_Vor_C1;
    int MIMCT = 5; // number of moves in a mean collision time (standard is 5)
    int MIMTT = 2; // number of moves in  mean transit time (standard is 2)
    MColTmin = MTrTXmin = MTrTYmin = 1e6 * DST1.DTM;
    //
    endTime = clock() - startTime;                          //an end in counting time
    DST1.time_taken = ((double)(endTime)) / CLOCKS_PER_SEC; // in seconds
    hour = DST1.time_taken / 3600;
    minute = (DST1.time_taken - hour * 3600) / 60;
    seconds = DST1.time_taken - hour * 3600 - minute * 60;

    //long double S_MUU;
    /*
MFP_VHS_Th_V: Mean Free Path based on thermal velocities at Variable Hard Sphere Gas
MFP_VHS_RMS_V:Mean Free Path based on Root-Mean-Square Velocities of Gaseous Particles
G_L_den:Length based on density gradient
*/
    //File Definition//
    V_Vor_C = 1e9;
    V_Vor_C1 = 1e9;
    FILE *fp, *fpp; //file with header suitable in tecplot
    char buf[65];
    //sprintf(buf, "../Output/Flowfield.dat");
    sprintf(buf, "../Output/Flowfield%d.dat", DST1.NPR);
    fp = fopen(buf, "w");
    fprintf(fp, "TITLE = \"Number of Particles: %d, Time_Steps: %lu, %dhours,%dminutes,%dseconds\"\n", DST1.NM, DST1.REPTN, hour, minute, seconds);
    fprintf(fp, "VARIABLES = \"X\",\"Y\",\"Density\",\"Temp(Trans)\",\"U\",\"V\",\"W\",\"Abs_Vel\",\"Mach\",\"Pressure\",\"CF_Ratio\",\"SOF\",\"Kn_thV_chrL\",\"Kn_thV_GLvel\",\"Kn_RMS_V_chrL\",\"Kn_RMS_V_GLvel\",\"HEAT_X\",\"HEAT_Y\",\"Ndensity\"\n");
    fprintf(fp, "ZONE T=\"KN=%Lf MAchi:%Lf Mols:%d FND=%Le DTM=%Le\" I =%d J=%d \n", Kn, DST1.Machi, DST1.NM, DST1.FND, DST1.DTM, DST1.NX, DST1.NY);
    fprintf(fp, "F= POINT\n");

    sprintf(buf, "../Output/Vort_Center.dat");
    fpp = fopen(buf, "w");
    //File Definition//
    //printf("\n Output the results of simulation\n");
    TIME_registered = DST1.TIME - DST1.TIMI;
    //Flowfield properties
    for (int j = 0; j < DST1.NY; j++)
    {
        for (int i = 0; i < DST1.NX; i++)
        {

            N = DST2.IJtC[i][j]; //cell number
                                 //ii is one horizental neighbor
                                 //if (i==0) ii=1; //right

            ii = i + 1;
            if (i == DST1.NX - 1)
                ii = DST1.NX - 2; //left
            if (i == 0)
                ii = 1;
            if (DST1.NX == 1)
                ii = 0;
            NN = DST2.IJtC[ii][j]; //cell number

            A = DST1.FNUM / (DST2.CC[i][j] * DST1.NSMP);
            A_NN = DST1.FNUM / (DST2.CC[ii][j] * DST1.NSMP);
            //for (int species=1,L)
            SN = DST2.CS[0][N];     //the number sum (single species)
            SN_NN = DST2.CS[0][NN]; //for HZ neighbor cell
            NP_Avg = SN / ((float)DST1.NSMP);
            SM = DST2.sp[4] * DST2.CS[0][N]; //sum of molecular mass (single species)
            SM_NN = DST2.sp[4] * DST2.CS[0][NN];
            //SMU:
            for (int k = 0; k <= 2; k++)
            {
                SMU[k] = DST2.sp[4] * DST2.CS[k + 1][N];
            } //SMU(1-3) are the sum of mu, mv, mw
            for (int k = 0; k <= 2; k++)
            {
                SMU_NN[k] = DST2.sp[4] * DST2.CS[k + 1][NN];
            } //SMU_NN(1-3) are the sum of mu, mv, mw

            SMCC = DST2.sp[4] * (DST2.CS[4][N] + DST2.CS[5][N] + DST2.CS[6][N]); //sum of m(u**2+v**2+w**2)
            SUM7 = (DST2.CS[4][N] + DST2.CS[5][N] + DST2.CS[6][N]);              //sum of (u**2+v**2+w**2)
            //End for (int species=1,L)
            DENN = SN * A;                                   //:Number Density, see eqn (1.34)
            DEN = DENN * DST2.sp[4];                         //: Density, see eqn (1.42)
            DEN_NN = SN_NN * A_NN * DST2.sp[4];              //density of neighbor
            G_L_den = fabsl(DEN_NN - DEN) / (DST1.CW * DEN); // Length based on density gradient
            //in case of different species, the mass is obtained through a division (averaging): {SM/SN}
            for (int k = 0; k <= 2; k++)
            {
                VEL[k] = SMU[k] / SM;
            }
            for (int k = 0; k <= 2; k++)
            {
                VEL_NN[k] = SMU_NN[k] / SM_NN;
            }

            UU = pow(VEL[0], 2.) + pow(VEL[1], 2.) + pow(VEL[2], 2.);
            UU_NN = pow(VEL_NN[0], 2.) + pow(VEL_NN[1], 2.) + pow(VEL_NN[2], 2.);

            TT = (SMCC - SM * UU) / (3.0 * BOLTZ * SN); //TT is the translational temperature, see eqn (1.51)
            //because only single monatomic species is selected then
            //only translational temp is considered,
            DST2.CT[i][j] = TT;
            XC = DST2.CG[6][i][j];
            YC = DST2.CG[7][i][j];
            ////////////////////////////////////Collision frequency//////////////////////////////////////
            //CF_Theor=4 n_d d_ref**2 sqrt(pi boltz T_ref/mass) (Temp/T_ref)**(1-omega)
            CF_Theo = 4 * DENN * (pow(DST2.sp[0], 2.)) * sqrtl(M_PI * BOLTZ * DST2.sp[1] / DST2.sp[4]) * pow((DST2.CT[i][j] / DST2.sp[1]), (1. - DST2.sp[2]));
            //CF_Numerical=NCOLL/(0.5*NPPCELL_AVG*TIME)
            CF_Num = DST2.ICol[0][N] / (0.5 * NP_Avg * TIME_registered);
            if (0.5 * NP_Avg * TIME_registered <= 0.0)
            {
                if (DST1.NPR > 0)
                {
                    printf("CF_Num has denominator zero\n");
                }
                CF_Num = 1;
            }
            CF_Ratio = CF_Num / CF_Theo;
            MCS = DST2.DCol[0][N] / DST2.ICol[0][N];
            MCT = (0.5 * NP_Avg * TIME_registered) / (DST2.ICol[0][N]); //Mean Collision Time
            //Local Mean Free Path based on equation 4-65: //VHS scheme based on thermal speeds
            MFP_VHS_Th_V = 1.0 / (sqrtl(2) * M_PI * (pow(DST2.sp[0], 2.)) * DENN * pow((DST2.sp[1] / DST2.CT[i][j]), (DST2.sp[2] - 0.50)));
            //MFP based on Root-Mean-Square Velocities of Gaseous Particles:
            MFP_VHS_RMS_V = 0.92132 * sqrtl(fabsl(SUM7 / SN - UU)) * MCT; //mean free path (based on r.m.s speed with correction factor for equilib.) DS2V definition
            MColT = MCT / (2.0 * (float)(MIMCT));
            if (MColT < DST1.DTM)
            {
                printf("Warning! MCT (%Le) smaller than DTM (%Le)\n", MColT, DST1.DTM);
                DST1.DTM = MColT;
            }
            //mean transit time along X//
            MTrTX = fabsl(DST1.CW / (VEL[0] * 2 * (float)(MIMTT)));
            MTrTY = fabsl(DST1.CH / (VEL[1] * 2 * (float)(MIMTT)));
            if (MTrTX < DST1.DTM)
            {
                printf("Warning! MTrTX (%Le) smaller than DTM (%Le)\n", MTrTX, DST1.DTM);
                DST1.DTM = MTrTX;
            }
            if (MTrTY < DST1.DTM)
            {
                printf("Warning! MTrTY (%Le) smaller than DTM (%Le)\n", MTrTY, DST1.DTM);
                DST1.DTM = MTrTY;
            }

            if (MColT < MColTmin)
            {
                MColTmin = MColT;
            }
            if (MTrTX < MTrTXmin)
            {
                MTrTXmin = MTrTX;
            }
            if (MTrTY < MTrTYmin)
            {
                MTrTYmin = MTrTY;
            }

            Abs_Vel = sqrtl(UU);
            //find minimum vel position
            if (V_Vor_C > Abs_Vel)
            {
                V_Vor_C = Abs_Vel;
                X_Vor_C = XC;
                Y_Vor_C = YC;
            }
            if ((XC > 0.12 * DST2.CB[0]) && (XC < 0.88 * DST2.CB[0]) && (YC > 0.12 * DST2.CB[1]) && (YC < 0.88 * DST2.CB[1]))
            {
                if (V_Vor_C1 > Abs_Vel)
                {
                    V_Vor_C1 = Abs_Vel;
                    X_Vor_C1 = XC;
                    Y_Vor_C1 = YC;
                }
            }
            //
            Abs_Vel_NN = sqrtl(UU_NN);
            G_L_den = fabsl(DEN_NN - DEN) / (DST1.CW * DEN);             // Length based on density gradient
            G_L_Vel = fabsl(Abs_Vel_NN - Abs_Vel) / (DST1.CW * Abs_Vel); // Length based on velocity gradient

            SOF = MCS / MFP_VHS_Th_V; //separation on free path based on thermal velocity
            Kn_thV_chrL = MFP_VHS_Th_V / DST1.CharLength;
            Kn_RMS_V_chrL = MFP_VHS_RMS_V / DST1.CharLength;

            Kn_thV_GLvel = MFP_VHS_Th_V / G_L_Vel; //Knudsen number based on local gradient of density
            Kn_RMS_V_GLvel = MFP_VHS_RMS_V / G_L_Vel;
            //MAch calculation
            Mach = Abs_Vel / sqrtl(5. / 3. * BOLTZ * DST2.CT[i][j] / DST2.sp[4]);
            Pressure = DENN * BOLTZ * TT; //p=nkT
            ////////heat_FLUX X&Y////////
            //directional pressures
            P_XX = DENN * DST2.sp[4] * (DST2.CS[4][N] / (float)(SN)-VEL[0] * VEL[0]);
            P_YY = DENN * DST2.sp[4] * (DST2.CS[5][N] / (float)(SN)-VEL[1] * VEL[1]);
            P_XY = DENN * DST2.sp[4] * (DST2.CS[7][N] / (float)(SN)-VEL[0] * VEL[1]);
            P_XZ = DENN * DST2.sp[4] * (DST2.CS[8][N] / (float)(SN)-VEL[0] * VEL[2]);
            P_YZ = DENN * DST2.sp[4] * (DST2.CS[9][N] / (float)(SN)-VEL[1] * VEL[2]);
            SMCC_U = DST2.sp[4] * (DST2.CS[10][N] + DST2.CS[13][N] + DST2.CS[16][N]) / ((float)SN); //sum of m(uuu+vvu+wwu)
            SMCC_V = DST2.sp[4] * (DST2.CS[11][N] + DST2.CS[14][N] + DST2.CS[17][N]) / ((float)SN); //sum of m(uuv+vvv+wwv)
            SMCC_B = SMCC / ((float)SN);
            Heat_x = (DENN / 2. * (SMCC_U - SMCC_B * VEL[0])) - P_XX * VEL[0] - P_XY * VEL[1] - P_XZ * VEL[2]; //Bird 2013 equ(1-22)
            Heat_y = (DENN / 2. * (SMCC_V - SMCC_B * VEL[1])) - P_XY * VEL[0] - P_YY * VEL[1] - P_YZ * VEL[2];

            //third order moments

            ///
            fprintf(fp, "%Le %Le %Le %Le %Le %Le %Le %Le %Lf %Le %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Le\n", XC / DST2.CB[0], YC / DST2.CB[1],
                    DEN, DST2.CT[i][j], VEL[0], VEL[1], VEL[2], Abs_Vel, Mach, Pressure,
                    CF_Ratio, SOF, Kn_thV_chrL, Kn_thV_GLvel, Kn_RMS_V_chrL, Kn_RMS_V_GLvel, Heat_x, Heat_y, DENN);

        } //NX
    }     //NY

    fclose(fp);
    fprintf(fpp, "\n%Lf %Lf %Lf\n", X_Vor_C / DST2.CB[0], Y_Vor_C / DST2.CB[1], V_Vor_C);
    fprintf(fpp, "\n%Lf %Lf %Lf\n", X_Vor_C1 / DST2.CB[0], Y_Vor_C1 / DST2.CB[1], V_Vor_C1);
    fclose(fpp);
    printf("DTM, MColT/DTM, MTrTX/DTM, MTrTY/DTM:  %Le  %Lf  %Lf  %Lf\n", DST1.DTM, MColTmin / DST1.DTM, MTrTXmin / DST1.DTM, MTrTYmin / DST1.DTM);
    printf("%d hours %d minutes %d seconds\n\n\n", hour, minute, seconds);
}
///////

void TRAckParticle()
{
    int Ncell, Ni, Nj, i;
    FILE *fp;
    char buf[65];
    sprintf(buf, "../Output/Track_particle/PartTrack_iter(%ld).dat", DST1.REPTN);
    fp = fopen(buf, "w");
    fprintf(fp, "Particle_Number_idendity     vx(m/s)           vy(m/s)          position_x(m)       position_y(m) \n");
    for (i = 0; i < DST1.NM; i++)
    {
        fprintf(fp, "%d                          %Lf          %Lf          %.8Le                %.8Le\n", i, DST2.PV[0][i], DST2.PV[1][i], DST2.PP[0][i], DST2.PP[1][i]);
    }

    fclose(fp);
    ////////
}
void TRAckParticle_Groups()
{
    int Ncell, Ni, Nj, i;
    int detG;
    long double pos_x[groups], pos_y[groups], v_x[groups], v_y[groups], v_z[groups];
    long double v_xx[groups], v_yy[groups], v_zz[groups];
    long double temp[groups];
    double group_length = (float)(DST1.NM) / (float)(groups);
    double CC, UU;
    for (i = 0; i < groups; i++)
    {
        pos_x[i] = 0.0;
        pos_y[i] = 0.0;
        v_x[i] = 0.0;
        v_y[i] = 0.0;
        v_z[i] = 0.0;
        v_xx[i] = 0.0;
        v_yy[i] = 0.0;
        v_zz[i] = 0.0;
        temp[i] = 0.0;
    }
    for (i = 0; i < DST1.NM; i++)
    {
        detG = (int)(i / div_group);
        pos_x[detG] += DST2.PP[0][i] / group_length;
        pos_y[detG] += DST2.PP[1][i] / group_length;
        v_x[detG] += DST2.PV[0][i] / group_length;
        v_y[detG] += DST2.PV[1][i] / group_length;
        v_z[detG] += DST2.PV[2][i] / group_length;
        v_xx[detG] += DST2.PV[0][i] * DST2.PV[0][i] / group_length;
        v_yy[detG] += DST2.PV[1][i] * DST2.PV[1][i] / group_length;
        v_zz[detG] += DST2.PV[2][i] * DST2.PV[2][i] / group_length;
    }
    CC = 0.0;
    UU = 0.0;
    for (i = 0; i < groups; i++)
    {
        CC = v_xx[i] + v_yy[i] + v_zz[i];
        UU = v_x[i] * v_x[i] + v_y[i] * v_y[i] + v_z[i] * v_z[i];
        temp[i] = DST2.sp[4] * (CC - UU) / (3.0 * BOLTZ);
    }
    printf("Groups %f \n", group_length);

    FILE *fp;
    char buf[65];
    sprintf(buf, "../Output/Track_particle/GroupPartTrack_iter(%ld).dat", DST1.REPTN);
    fp = fopen(buf, "w");
    fprintf(fp, "Group_Number_idendity         Temp(K)          vx(m/s)           vy(m/s)          position_x(m)       position_y(m) \n");
    for (i = 0; i < groups; i++)
    {
        fprintf(fp, "%d                       %Lf                %Lf               %Lf                %Lf                %Lf\n", i, temp[i], v_x[i], v_y[i], pos_x[i], pos_y[i]);
    }

    fclose(fp);
    ////////
}

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
void initAftRead()
{
    int MM, NCELL, Nj, Ni;
    long double A, REM, SC;
    DST1.IFSurf = 0;
    if ((SurfDist >= 0) && (SurfDist <= 1.0))
    {
        DST1.IFSurf = 1;
    } //if there is additional surface inside

    if (DST1.IFSurf == 1)
    {
        DST2.IB[4] = ib5;
        DST2.TSURF[4] = Tsurf5;
    }
    DST1.SurfHeightL = SurfHeight * DST2.CB[1];
    DST1.SurfDistL = SurfDist * DST2.CB[0];
    DST2.TSURF[0] = Tsurf1;
    DST2.TSURF[1] = Tsurf2;
    DST2.TSURF[2] = Tsurf3;
    DST2.TSURF[3] = Tsurf4;

    DST2.sp[0] = sp1D;
    DST2.sp[1] = sp2T;
    DST2.sp[2] = sp3VT;
    DST2.sp[3] = sp4ReVSS;
    DST2.sp[4] = sp5M;

    DST1.SPM0 = M_PI * DST2.sp[0] * DST2.sp[0];                      //collision cross section (1-35)
    DST1.SPM1 = DST2.sp[1];                                          //the reference temperature
    DST1.SPM2 = DST2.sp[2];                                          //the viscosity-temperature power law
    DST1.SPM3 = DST2.sp[3];                                          //the reciprocal of the VSS scattering parameter
    DST1.SPM4 = DST2.sp[4] * DST2.sp[4] / (DST2.sp[4] + DST2.sp[4]); //the reduced mass is defined in eqn (2.7)
    DST1.SPM5 = GAM(2.5 - DST1.SPM2);                                //the Gamma function of (5/2 - viscosity-temperature power law)

    DST1.C_Soundi = sqrtl(5. / 3. * BOLTZ * FTMP / DST2.sp[4]); //initial sound speed:sqrtl(gama k T / m):322.6752
    DST1.VMP = sqrtl(2. * BOLTZ * FTMP / DST2.sp[4]);

    //calculate the number of particles that enter at each time step
    //4 borders: Xbtm/upp Ybtm/upp IB0123
    for (int k = 0; k < 4; k++)
    {
        if (DST2.IB[k] == 1)
        {
            printf("\n Side %d stream set up\n", k);
            DST1.VMP = sqrtl(2. * BOLTZ * FTMP / DST2.sp[4]);
            if (k == 0)
                SC = VFX / DST1.VMP;
            if (k == 1)
                SC = -VFX / DST1.VMP;
            if (k == 2)
                SC = VFY / DST1.VMP;
            if (k == 3)
                SC = -VFY / DST1.VMP;
            if (fabsl(SC) < 10.1)
            {
                A = (exp(-SC * SC) + SPI * SC * (1. + ERF(SC))) / (2. * SPI);
            }
            if (SC > 10.)
            {
                A = SC;
            }
            if (SC < -10.)
            {
                A = 0.;
            } //A is the non-dimensional flux of eqn (4.22)

            if ((k == 0) || (k == 1))
            {
                DST2.BME[k] = DST1.FND * A * DST1.VMP * DST1.DTM * DST1.FH / DST1.FNUM;
            }
            else
            {
                DST2.BME[k] = DST1.FND * A * DST1.VMP * DST1.DTM * DST1.FW / DST1.FNUM;
            }
            printf("entering mols %Lf\n", DST2.BME[k]);
        }
    }
}
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

///////////////////MT Random Generator////////////////////////
#define RNN 312
#define RMM 156
#define MATRIX_A 0xB5026F5AA96619E9ULL
#define UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define LM 0x7FFFFFFFULL         /* Least significant 31 bits */

/* The array for the state vector */
static unsigned long long mt[RNN];
/* mti==NN+1 means mt[NN] is not initialized */
static int mti = RNN + 1;

/* initializes mt[NN] with a seed */
void init_genrand64(unsigned long long seed)
{
    if (Iseeding == 0)
        mt[0] = seed; //A constant SEED
    if (Iseeding == 1)
        mt[0] = seed * drand48() + 1ULL; //Varying seed for each run
    DST1.RandSeed = mt[0];
    for (mti = 1; mti < RNN; mti++)
        mt[mti] = (6364136223846793005ULL * (mt[mti - 1] ^ (mt[mti - 1] >> 62)) + mti);
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
void init_by_array64(unsigned long long init_key[],
                     unsigned long long key_length)
{
    unsigned long long i, j, k;
    init_genrand64(19650218ULL);
    i = 1;
    j = 0;
    k = (RNN > key_length ? RNN : key_length);
    for (; k; k--)
    {
        mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 62)) * 3935559000370003845ULL)) + init_key[j] + j; /* non linear */
        i++;
        j++;
        if (i >= RNN)
        {
            mt[0] = mt[RNN - 1];
            i = 1;
        }
        if (j >= key_length)
            j = 0;
    }
    for (k = RNN - 1; k; k--)
    {
        mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 62)) * 2862933555777941757ULL)) - i; /* non linear */
        i++;
        if (i >= RNN)
        {
            mt[0] = mt[RNN - 1];
            i = 1;
        }
    }

    mt[0] = 1ULL << 63; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0, 2^64-1]-interval */
unsigned long long genrand64_int64(void)
{
    int i;
    unsigned long long x;
    static unsigned long long mag01[2] = {0ULL, MATRIX_A};

    if (mti >= RNN)
    { /* generate NN words at one time */

        /* if init_genrand64() has not been called, */
        /* a default initial seed is used     */
        if (mti == RNN + 1)
            init_genrand64(15489ULL); //init_genrand64(5489ULL);

        for (i = 0; i < RNN - RMM; i++)
        {
            x = (mt[i] & UM) | (mt[i + 1] & LM);
            mt[i] = mt[i + RMM] ^ (x >> 1) ^ mag01[(int)(x & 1ULL)];
        }
        for (; i < RNN - 1; i++)
        {
            x = (mt[i] & UM) | (mt[i + 1] & LM);
            mt[i] = mt[i + (RMM - RNN)] ^ (x >> 1) ^ mag01[(int)(x & 1ULL)];
        }
        x = (mt[RNN - 1] & UM) | (mt[0] & LM);
        mt[RNN - 1] = mt[RMM - 1] ^ (x >> 1) ^ mag01[(int)(x & 1ULL)];

        mti = 0;
    }

    x = mt[mti++];

    x ^= (x >> 29) & 0x5555555555555555ULL;
    x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
    x ^= (x << 37) & 0xFFF7EEE000000000ULL;
    x ^= (x >> 43);

    return x;
}

/* generates a random number on [0, 2^63-1]-interval */
long long genrand64_int63(void)
{
    return (long long)(genrand64_int64() >> 1);
}

/* generates a random number on [0,1]-real-interval */
double genrand64_real1(void)
{
    return (genrand64_int64() >> 11) * (1.0 / 9007199254740991.0);
}

/* generates a random number on [0,1)-real-interval */
double genrand64_real2(void)
{
    return (genrand64_int64() >> 11) * (1.0 / 9007199254740992.0);
}

/* generates a random number on (0,1)-real-interval */
double genrand64_real3(void)
{
    return ((genrand64_int64() >> 12) + 0.5) * (1.0 / 4503599627370496.0);
}
///////////////////MT Random Generator////////////////////////
