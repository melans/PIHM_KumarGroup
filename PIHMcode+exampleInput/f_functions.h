# ifndef FUNCTIONS_H_
# define FUNCTIONS_H_
#include "pihm.h"
realtype returnVal(realtype rArea, realtype rPerem, realtype eqWid,realtype ap_Bool);
realtype CS_AreaOrPerem(int rivOrder, realtype rivDepth, realtype rivCoeff, realtype a_pBool);
void OverlandFlow(realtype **flux, int loci, int locj, realtype avg_y, realtype grad_y, realtype avg_sf, realtype crossA, realtype avg_rough);
realtype avgY(realtype diff, realtype yi, realtype yinabr);
realtype effKV(realtype ksatFunc,realtype gradY,realtype macKV,realtype KV,realtype areaF);
realtype effKH(int mp,realtype tmpY, realtype aqDepth, realtype MacD, realtype MacKsatH, realtype areaF, realtype ksatH);
realtype Interpolation(TSD *Data, realtype t);
realtype GradCalc(realtype yNabr, realtype zNabr, realtype y, realtype z, realtype dist, int  interp, int interpNabr, realtype coeff, realtype coeffNabr, int booltype);
realtype avgKH(int macp,realtype y,realtype aqD,realtype macD,realtype macKsatH,realtype vAreaF,realtype KsatH,int macpnabr,realtype ynabr,realtype aqDnabr,realtype macDnabr,realtype macKsatHnabr,realtype vAreaFnabr,realtype KsatHnabr);
void fluxCalc_Ele(Model_Data DS, int i,realtype t);
void fluxCalc_Riv(Model_Data DS, int i,realtype t);
void flux_cal(realtype, N_Vector, void *);
#endif
