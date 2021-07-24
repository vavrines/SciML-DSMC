//******************************************************************************
//Title:
//        Direct Simulation Monte Carlo (DSMC) program for two-dimensional flows
//Version:
//        Version 2.3.05
//History:
//        Written by Bijan Gosaheyshi, 7/2021.
//******************************************************************************

#include "dsmcINCLUDE.h"
/////
int main()
{
    dsmcINITIALIZATION();

    //New_start or continue
    if (NEWsimulation == 1)
    {
        printf("\n=============================================================================\n");
        printf("\n New Simulation, varialbes are initialized\n");
        printf("\n=============================================================================\n");

        INIT();
    }
    else if (NEWsimulation == 0)
    {
        printf("\n=============================================================================\n");
        printf("\n Continue Simulation, variables are read from \"dsmcBinary.bin\" \n");
        Read();
        printf("\n Start at NPR: %d\n", DST1.NPR);
        printf("\n Please Confirm, otherwise stop\n");
        printf("\n=============================================================================\n");
    }
    //New_Sample or continue
    if (NEWsample == 1)
    {
        printf("\n=============================================================================\n");
        printf("\n New Samples requested: available SAMPLES are INITIALIZED >> SAMPI(); \n");
        SAMPI();
        printf("\n=============================================================================\n");
    }
    else if (NEWsample == 0)
    {
        printf("\n=============================================================================\n");
        printf("\n Available SAMPLES are CONTINUED\n");
        ContinueSample();
        printf("\n=============================================================================\n");
    }
    printf("\n START@ Iter:%d, Samp:%d \n", DST1.NPR, DST1.NSMP);
    printf("\n TOTAL iteration will BE: %d\n", NPT);
    printf("\n\n Loops treated by NIS:%d, NSP:%d, NPS(SAMPI):%d\n", NIS, NSP, NPS);

    /////////////////
    do
    {

        DST1.NPR = DST1.NPR + 1;
        //total number of <moves-Index-Collision>:(NPR)*(NSP)*(NIS)
        if (DST1.NPR <= NPS)
            SAMPI();

        for (int JJJ = 1; JJJ <= NSP; JJJ++)
        {
            for (int III = 1; III <= NIS; III++)
            {
                DST1.TIME = DST1.TIME + DST1.DTM;
                Move();
                INDEX();
                Collision();

            } //time loop
            Sample();
        } //sample loop

        //if (    ((DST1.NPR%200)==0)||(DST1.NPR==1) ) {
        printf("\nOUTPUT@ Iter:%d of %d, Samp:%d, Colls: %lu, Movs: %lu for %d Particls\n", DST1.NPR, NPT, DST1.NSMP, DST1.NCOL, DST1.MOVT, DST1.NM);
        printf("Total number of <MOVE_INDEX_COLLISION>: %d, %ld\n", DST1.NPR * NSP * NIS, DST1.REPTN);

        WriteOut();
        OUTPUT();
        //}
    } while (DST1.NPR < NPT);

    printf("End of the Program at iteration: %d\n", DST1.NPR);
    return 0;
}
