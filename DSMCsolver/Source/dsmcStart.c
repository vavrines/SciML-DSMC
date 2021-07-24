//******************************************************************************
//Title:
//        Direct Simulation Monte Carlo (DSMC) program for two-dimensional flows
//Version:
//        Version2 2.3.05
//History:
//        Written by Bijan Gosaheyshi, 7/2021.
//******************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

long double Kn, ppcel, cb1, cb2, sp1D, sp2T, sp3VT, sp4ReVSS, sp5M;
long double WallvelX, WallvelY, VFX, VFY, GravX, GravY, FTMP, TUSRF1, TUSRF2, TUSRF3, TUSRF4, TUSRF5, SpecDist, SurfDist, SurfHeight;
int NEWsample, NEWsimulation, Init_AftRead, charL, MinDevisionX, MinDevisionY;
int IB1, IB2, IB3, IB4, NIS, NSP, NPS, Irandom, Iseeding;
unsigned int NPT;
int main()
{
    char string[200];

    FILE *fpm;
    fpm = fopen("../Data/pre-program_DSMC.txt", "r");
    fgets(string, 200, fpm);
    sscanf(string, "%d", &NEWsimulation); //  1/0 New/Continue
    //fscanf(fpm, "%d \n", &NEWsimulation);       //  1/0 New/Continue
    fgets(string, 200, fpm);
    sscanf(string, "%d \n", &NEWsample); //  1/0 New/Continue
    //fscanf(fpm, "%d \n", &NEWsample);           //  1/0 New/Continue
    fgets(string, 200, fpm);
    sscanf(string, "%d \n", &Init_AftRead);
    fgets(string, 200, fpm);
    sscanf(string, "%Lf \n", &Kn);
    //fscanf(fpm, "%Le \n", &Kn);
    fgets(string, 200, fpm);
    sscanf(string, "%Lf \n", &ppcel);
    //fscanf(fpm, "%Le \n", &ppcel);
    fgets(string, 200, fpm);
    sscanf(string, "%Le \n", &cb1);
    fgets(string, 200, fpm);
    sscanf(string, "%Le \n", &cb2);
    fgets(string, 200, fpm);
    sscanf(string, "%d \n", &IB1);
    fgets(string, 200, fpm);
    sscanf(string, "%d \n", &IB2);

    //fscanf(fpm, "%Le \n", &cb1);
    //fscanf(fpm, "%Le \n", &cb2);
    //fscanf(fpm, "%d \n", &IB1);
    //fscanf(fpm, "%d \n", &IB2);
    if ((IB1 == 3) && (IB2 != 3))
    {
        IB2 = IB1;
        printf("\nAt periodic Boundary:0-1_X_direction, IB values are not equal to 3!\n");
        printf("\nCorrection is then made so that IBs are equal to 3!\n");
    }
    fgets(string, 200, fpm);
    sscanf(string, "%d \n", &IB3);
    fgets(string, 200, fpm);
    sscanf(string, "%d \n", &IB4);
    //fscanf(fpm, "%d \n", &IB3);
    //fscanf(fpm, "%d \n", &IB4);
    if ((IB3 == 3) && (IB4 != 3))
    {
        IB4 = IB3;
        printf("\nAt periodic Boundary:2-3_Y_direction, IB values are not equal to 3!\n");
        printf("\nCorrection is then made so that IBs are equal to 3!\n");
    }
    fgets(string, 200, fpm);
    sscanf(string, "%Lf \n", &TUSRF1);
    fgets(string, 200, fpm);
    sscanf(string, "%Lf \n", &TUSRF2);
    fgets(string, 200, fpm);
    sscanf(string, "%Lf \n", &TUSRF3);
    fgets(string, 200, fpm);
    sscanf(string, "%Lf \n", &TUSRF4);
    fgets(string, 200, fpm);
    sscanf(string, "%Le \n", &sp1D);
    fgets(string, 200, fpm);
    sscanf(string, "%Le \n", &sp2T);
    fgets(string, 200, fpm);
    sscanf(string, "%Le \n", &sp3VT);
    fgets(string, 200, fpm);
    sscanf(string, "%Le \n", &sp4ReVSS);
    fgets(string, 200, fpm);
    sscanf(string, "%Le \n", &sp5M);
    fgets(string, 200, fpm);
    sscanf(string, "%Le \n", &WallvelX);
    fgets(string, 200, fpm);
    sscanf(string, "%Le \n", &WallvelY);
    fgets(string, 200, fpm);
    sscanf(string, "%Le \n", &VFX);
    fgets(string, 200, fpm);
    sscanf(string, "%Le \n", &VFY);
    fgets(string, 200, fpm);
    sscanf(string, "%Le \n", &GravX);
    fgets(string, 200, fpm);
    sscanf(string, "%Le \n", &GravY);
    fgets(string, 200, fpm);
    sscanf(string, "%Le \n", &FTMP);
    fgets(string, 200, fpm);
    sscanf(string, "%d \n", &NIS);
    fgets(string, 200, fpm);
    sscanf(string, "%d \n", &NSP);
    fgets(string, 200, fpm);
    sscanf(string, "%d \n", &NPS);
    fgets(string, 200, fpm);
    sscanf(string, "%d \n", &NPT);
    fgets(string, 200, fpm);
    sscanf(string, "%d \n", &Irandom);
    fgets(string, 200, fpm);
    sscanf(string, "%d \n", &Iseeding);
    fgets(string, 200, fpm);
    sscanf(string, "%Le \n", &SpecDist);
    fgets(string, 200, fpm);
    sscanf(string, "%Le \n", &SurfDist);
    fgets(string, 200, fpm);
    sscanf(string, "%Le \n", &SurfHeight);
    if (SurfDist < 0)
        SurfHeight = 0;
    fgets(string, 200, fpm);
    sscanf(string, "%d \n", &charL);
    fgets(string, 200, fpm);
    sscanf(string, "%d \n", &MinDevisionX);
    fgets(string, 200, fpm);
    sscanf(string, "%d \n", &MinDevisionY);
    fgets(string, 200, fpm);
    sscanf(string, "%Lf \n", &TUSRF5);
    //fscanf(fpm, "%Lf \n", &TUSRF5);
    printf("Knudsen= %Le,  ppcel= %Le %Le \n%Le %Le %Le %Le %Le %Le \n", Kn, ppcel, cb1, cb2, sp1D, sp2T, sp3VT, sp4ReVSS, sp5M);
    printf("NIS,NSP,NPS,NPT= %d %d %d %d\n", NIS, NSP, NPS, NPT);

    //getchar();
    fclose(fpm);
    FILE *f1pm;
    f1pm = fopen("dsmcBasicInput.h", "w");

    fprintf(f1pm, "#ifndef _STDIO_H \n");
    fprintf(f1pm, "#define _STDIO_H \n");
    fprintf(f1pm, "#endif \n");

    fprintf(f1pm, "#ifndef _STDLIB_H \n");
    fprintf(f1pm, "#define _STDLIB_H \n");
    fprintf(f1pm, "#endif \n");

    fprintf(f1pm, "#ifndef _STRING_H \n");
    fprintf(f1pm, "#define _STRING_H \n");
    fprintf(f1pm, "#endif \n");

    fprintf(f1pm, "#ifndef _MATH_H \n");
    fprintf(f1pm, "#define _MATH_H \n");
    fprintf(f1pm, "#endif \n");

    fprintf(f1pm, "#ifndef _TGMATH_H \n");
    fprintf(f1pm, "#define _TGMATH_H \n");
    fprintf(f1pm, "#endif \n");

    fprintf(f1pm, "#ifndef _FLOAT_H \n");
    fprintf(f1pm, "#define _FLOAT_H \n");
    fprintf(f1pm, "#endif \n");

    fprintf(f1pm, "#ifndef _TIME_H \n");
    fprintf(f1pm, "#define _TIME_H \n");
    fprintf(f1pm, "#endif \n");

    fprintf(f1pm, "#define NEWsimulation  %d\n", NEWsimulation);
    fprintf(f1pm, "#define NEWsample %d\n", NEWsample);
    fprintf(f1pm, "#define Init_AftRead  %d\n", Init_AftRead);
    fprintf(f1pm, "#define Kn (long double)  %Le\n", Kn);
    fprintf(f1pm, "#define ppcel (long double)  %Le\n", ppcel);
    fprintf(f1pm, "#define cb1 (long double)  %Le\n", cb1);
    fprintf(f1pm, "#define cb2 (long double)  %Le\n", cb2);
    fprintf(f1pm, "#define ib1 %d   //x leftt border    /// IB=1 Stream; IB=2 Surface; IB=3 Periodic:\n", IB1);
    fprintf(f1pm, "#define ib2 %d   //x right border\n", IB2);
    fprintf(f1pm, "#define ib3 %d   //y downn border\n", IB3);
    fprintf(f1pm, "#define ib4 %d   //y upppp border\n", IB4);

    fprintf(f1pm, "#define Tsurf1 %Lf   //x leftt border  //TUSRF (+: diffuse reflection), (-: specular reflection)\n", TUSRF1);
    fprintf(f1pm, "#define Tsurf2 %Lf   //x right border\n", TUSRF2);
    fprintf(f1pm, "#define Tsurf3 %Lf  //y downn border\n", TUSRF3);
    fprintf(f1pm, "#define Tsurf4 %Lf   //y downn border\n", TUSRF4);
    fprintf(f1pm, "#define sp1D (long double)  %Le //the reference cross-section (diameter in the data)\n", sp1D);
    fprintf(f1pm, "#define sp2T (long double)  %Le //the reference temperature \n", sp2T);
    fprintf(f1pm, "#define sp3VT (long double)  %Le //the viscosity-temperature power law (Omega) \n", sp3VT);
    fprintf(f1pm, "#define sp4ReVSS (long double)  %Le //the reciprocal of the VSS scattering parameter, here only VHS:1 \n", sp4ReVSS);
    fprintf(f1pm, "#define sp5M (long double)  %Le //the molecular mass\n", sp5M);
    fprintf(f1pm, "#define WallvelX (long double)  %Le\n", WallvelX);
    fprintf(f1pm, "#define WallvelY (long double)  %Le\n", WallvelY);
    fprintf(f1pm, "#define VFX (long double)  %Le\n", VFX); //side 0-3
    fprintf(f1pm, "#define VFY (long double)  %Le\n", VFY);
    fprintf(f1pm, "#define GravX (long double)  %Le\n", GravX); //side 0-3
    fprintf(f1pm, "#define GravY (long double)  %Le\n", GravY);

    fprintf(f1pm, "#define FTMP (long double)  %Le\n", FTMP);
    fprintf(f1pm, "#define NIS %d\n", NIS);
    fprintf(f1pm, "#define NSP %d\n", NSP);
    fprintf(f1pm, "#define NPS %d\n", NPS);
    fprintf(f1pm, "#define NPT %d\n", NPT);
    fprintf(f1pm, "#define Irandom %d\n", Irandom);
    fprintf(f1pm, "#define Iseeding %d\n", Iseeding);
    fprintf(f1pm, "#define SpecDist %Lf   //fraction where ib0 is specular\n", SpecDist);
    fprintf(f1pm, "#define SurfDist %Lf   //fraction X where a vert surface is defined\n", SurfDist);
    fprintf(f1pm, "#define SurfHeight %Lf //fraction Y where height of vert surface is defined\n", SurfHeight);
    fprintf(f1pm, "#define ib5 10   //vert surface possible presence if ib5=10\n");
    fprintf(f1pm, "#define charL %d\n", charL);
    fprintf(f1pm, "#define MinDevisionX %d\n", MinDevisionX);
    fprintf(f1pm, "#define MinDevisionY %d\n", MinDevisionY);
    fprintf(f1pm, "#define Tsurf5 %Lf  \n", TUSRF5);

    fclose(f1pm);

    return 0;
}
