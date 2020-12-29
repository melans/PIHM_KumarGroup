#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "nvector_serial.h"
#include "sundialstypes.h"   
#include "pihm.h"

realtype effKV_var(realtype ksatFunc, realtype gradY, realtype macKV, realtype KV, realtype areaF)
{
	if (ksatFunc >= 0.98)
	{
		return (macKV*areaF + KV*(1 - areaF)*ksatFunc);
	}
	else
	{
		if (fabs(gradY)*ksatFunc*KV <= 1 * KV*ksatFunc)
		{
			return KV*ksatFunc;
		}
		else
		{
			if (fabs(gradY)*ksatFunc*KV<(macKV*areaF + KV*(1 - areaF)*ksatFunc))
			{
				return (macKV*areaF*ksatFunc + KV*(1 - areaF)*ksatFunc);
			}
			else
			{
				return (macKV*areaF + KV*(1 - areaF)*ksatFunc);
			}
		}
	}
}

void var(Model_Data DS, realtype t, N_Vector CV_Y)
	{
		realtype elemSatn,TotalY_unsat,effK_unsat,satKfunc,tmpVar,effKI,AquiferDepth;
		int i, i_plus_ele, i_plus_2ele, i_plus_3ele;
		realtype *Y;
	//	Y = NV_DATA_S(CV_Y);

		for (i = 0; i<DS->totele; i++)
		{
			DS->DummyY[i] = (NV_Ith_S(CV_Y,i) >= 0) ? NV_Ith_S(CV_Y,i) : 0;
		}

#ifndef NO_UNSAT
#ifdef SUB_SURF_RIV
		for (i = 0; i < DS->NumEle; i++)
		{
			i_plus_ele = i + DS->NumEle;
			i_plus_2ele = i + 2 * DS->NumEle;
			AquiferDepth = (DS->Ele[i].zmax - DS->Ele[i].zmin);

			/* Note: Assumption is OVL flow depth less than EPS_3 is immobile water */
			/* Calculation of saturation of top layer */
			elemSatn = DS->DummyY[i_plus_ele] / DS->Ele[i].infD;
			elemSatn = ((elemSatn<0.1) ? 0.1 : ((elemSatn>1.0) ? 1 : elemSatn));

			Avg_Y = -(pow(pow(1 / elemSatn, DS->Ele[i].Beta / (DS->Ele[i].Beta - 1)) - 1, 1 / DS->Ele[i].Beta) / DS->Ele[i].Alpha);
			TotalY_unsat = ((elemSatn >= 0.5) ? (DS->DummyY[i_plus_ele] - 0.5*DS->Ele[i].infD) : 0) + ((Avg_Y < MINpsi) ? MINpsi : Avg_Y) + DS->Ele[i].zmin + AquiferDepth - 0.5*DS->Ele[i].infD;                /* Calculation of gradient between overland flow and top unsaturated zone */
			Grad_Y = (DS->DummyY[i] + DS->Ele[i].zmax - TotalY_unsat) / (0.5*DS->Ele[i].infD);
			Grad_Y = (Grad_Y>0) ? ((DS->DummyY[i]<EPS_3) ? 0 : Grad_Y) : ((DS->DummyY[i_plus_ele]>EPS_3) ? Grad_Y : 0);
			/* The effective conductivity for infiltration is calculated based on saturated condition */
			/*
							Following are the expressions being evaluated in the next line
							elemSatn=1.0;
							satKfunc=pow(elemSatn,0.5)*pow(-1+pow(1-pow(elemSatn,DS->Ele[i].Beta/(DS->Ele[i].Beta-1)),(DS->Ele[i].Beta-1)/DS->Ele[i].Beta),2);
							satKfunc=1.0;
							effK_unsat=(DS->Ele[i].Macropore==1)?effKV(satKfunc,Grad_Y,DS->Ele[i].macKsatV,DS->Ele[i].infKsatV,DS->Ele[i].hAreaF):DS->Ele[i].infKsatV;
							*/
			effK_unsat = (DS->Ele[i].Macropore == 1) ? DS->Ele[i].macKsatV*DS->Ele[i].hAreaF + (1 - DS->Ele[i].hAreaF)*DS->Ele[i].infKsatV : DS->Ele[i].infKsatV;
			/* Infiltration/exfiltration rate calculated */
		
			DS->EleViR[i] = ((DS->EleNetPrep[i] <= (effK_unsat*Grad_Y)) && (DS->DummyY[i] <= EPS_3)) ? DS->EleNetPrep[i] : (effK_unsat*Grad_Y);
		
			//		DS->EleViR[i] = (effK_unsat)*Grad_Y;
			/* Calculation of effective conductivity of top layer */
			elemSatn = DS->DummyY[i_plus_ele] / DS->Ele[i].infD;
			elemSatn = (elemSatn < 0.1) ? 0.1 : (elemSatn>1 ? 1 : elemSatn);
			//             satKfunc=pow(elemSatn,0.5)*pow(-1+pow(1-pow(elemSatn,DS->Ele[i].Beta/(DS->Ele[i].Beta-1)),(DS->Ele[i].Beta-1)/DS->Ele[i].Beta),2);
			effK_unsat = (DS->Ele[i].Macropore == 1) ? effKV_var(elemSatn, Grad_Y, DS->Ele[i].macKsatV, DS->Ele[i].infKsatV, DS->Ele[i].hAreaF) : DS->Ele[i].infKsatV*elemSatn;
			//              effK_unsat=(DS->Ele[i].Macropore==1)?effKV(satKfunc,Grad_Y,DS->Ele[i].macKsatV,DS->Ele[i].infKsatV,DS->Ele[i].hAreaF):DS->Ele[i].infKsatV*satKfunc;        
#ifdef LAYER3
			/* Calculation of saturation of intermediate layer */
			i_plus_3ele = i + 3*DS->NumEle;
			tmpVar = ((AquiferDepth - DS->Ele[i].infD - DS->DummyY[i_plus_2ele]) < DS->Ele[i].infD) ? DS->Ele[i].infD : (AquiferDepth - DS->Ele[i].infD - DS->DummyY[i_plus_2ele]);
			elemSatn = DS->DummyY[i_plus_3ele] / tmpVar;
			elemSatn = (elemSatn < 0.1) ? 0.1 : (elemSatn>1 ? 1 : elemSatn);
			Avg_Y = -(pow(pow(1 / elemSatn, DS->Ele[i].Beta / (DS->Ele[i].Beta - 1)) - 1, 1 / DS->Ele[i].Beta) / DS->Ele[i].Alpha);
			TotalY = ((Avg_Y<MINpsi) ? MINpsi : Avg_Y) + DS->Ele[i].zmin - 0.5*tmpVar + AquiferDepth - DS->Ele[i].infD + (((elemSatn>0.5) && (tmpVar > DS->Ele[i].infD)) ? (DS->DummyY[i_plus_3ele] - 0.5*(AquiferDepth - DS->Ele[i].infD - DS->DummyY[i_plus_2ele])) : 0) + ((tmpVar == DS->Ele[i].infD) ? ((DS->DummyY[i_plus_3ele] - 0.5*tmpVar) > 0 ? (DS->DummyY[i_plus_3ele] - 0.5*tmpVar) : 0.0*tmpVar) : 0);
			/* Calculation of gradient between top and intermediate layer */
			Grad_Y = (TotalY_unsat - TotalY) / (0.5*(DS->Ele[i].infD + tmpVar));

			/* Calculation of effective conductivity of the intermediate layer */
			//                satKfunc=pow(elemSatn,0.5)*pow(-1+pow(1-pow(elemSatn,DS->Ele[i].Beta/(DS->Ele[i].Beta-1)),(DS->Ele[i].Beta-1)/DS->Ele[i].Beta),2);
			effKI = (DS->Ele[i].Macropore == 1) ? effKV_var(elemSatn, Grad_Y, DS->Ele[i].macKsatV, DS->Ele[i].KsatV, DS->Ele[i].hAreaF) : DS->Ele[i].KsatV*elemSatn;
			//              effKI=(DS->Ele[i].Macropore==1)?effKV(satKfunc,Grad_Y,DS->Ele[i].macKsatV,DS->Ele[i].KsatV,DS->Ele[i].hAreaF):DS->Ele[i].KsatV*satKfunc;
			/* Calculation of exchange flux between top and intermediate layer */
			/* Using upwind method */
			DS->RechargeI[i] = (Grad_Y > 0) ? ((DS->DummyY[i_plus_ele] > EPS_3) ? effK_unsat*Grad_Y : 0) : ((DS->DummyY[i_plus_3ele] > EPS_3) ? effKI*Grad_Y : 0);
			/* Calculation of gradient between intermediate and bottom layer */
			Grad_Y = (TotalY - (DS->DummyY[i_plus_2ele] + DS->Ele[i].zmin)) / (0.5*(AquiferDepth - DS->Ele[i].infD));
			/* Calculation of effective conductivity of the bottom layer */
			effKI = (DS->Ele[i].Macropore == 1) ? effKV_var(elemSatn, Grad_Y, DS->Ele[i].macKsatV, DS->Ele[i].KsatV, DS->Ele[i].hAreaF) : DS->Ele[i].KsatV*elemSatn;
			tmpVar = (DS->DummyY[i_plus_2ele] > AquiferDepth - DS->Ele[i].macD) ? effKI : DS->Ele[i].KsatV*elemSatn;
			//                tmpVar=(DS->DummyY[i_plus_2ele]>AquiferDepth-DS->Ele[i].macD)?effKI:DS->Ele[i].KsatV*satKfunc;
			/*
							Following statements will be rewritted to avoid repetetive computations
							effKI calculates top layer effective conductivity
							satKfunc=1.0;
							effKI=(DS->Ele[i].Macropore==1)?effKV(satKfunc,Grad_Y,DS->Ele[i].macKsatV,DS->Ele[i].KsatV,DS->Ele[i].hAreaF):DS->Ele[i].KsatV*satKfunc;
							effK=(DS->Ele[i].Macropore==1)?((DS->DummyY[i_plus_2ele]>AquiferDepth-DS->Ele[i].macD)?effKI:DS->Ele[i].KsatV*satKfunc):DS->Ele[i].KsatV*satKfunc;
							*/
			effKI = (DS->Ele[i].Macropore == 1) ? DS->Ele[i].macKsatV*DS->Ele[i].hAreaF + (1 - DS->Ele[i].hAreaF)*DS->Ele[i].KsatV : DS->Ele[i].KsatV;
			effK = (DS->Ele[i].Macropore == 1) ? ((DS->DummyY[i_plus_2ele] > AquiferDepth - DS->Ele[i].macD) ? effKI : DS->Ele[i].KsatV) : DS->Ele[i].KsatV;
			/* Calculation of exchange flux between top and bottom layer */
			/* Using upwind method */
			DS->Recharge[i] = (Grad_Y > 0) ? ((DS->DummyY[i_plus_3ele] > EPS_3) ? tmpVar*Grad_Y : 0) : ((DS->DummyY[i_plus_2ele] > EPS_3) ? effK*Grad_Y : 0);
#else
			/*                Effective conductivity of the top layer. The conductivity should be infKsatV assuming that this is the conductivity of the top layer. That is why the line below is commented

			*/
			//              effK_unsat=(DS->Ele[i].Macropore==1)?effKV(satKfunc,Grad_Y,DS->Ele[i].macKsatV,DS->Ele[i].KsatV,DS->Ele[i].hAreaF):DS->Ele[i].KsatV*satKfunc;
			elemSatn = DS->DummyY[i_plus_2ele] / (AquiferDepth - DS->Ele[i].infD);
			elemSatn = (elemSatn<0.1) ? 0.1 : (elemSatn>1.0 ? 1.0 : elemSatn);
			//                satKfunc=pow(elemSatn,0.5)*pow(-1+pow(1-pow(elemSatn,DS->Ele[i].Beta/(DS->Ele[i].Beta-1)),(DS->Ele[i].Beta-1)/DS->Ele[i].Beta),2);
			Avg_Y = -(pow(pow(1 / elemSatn, DS->Ele[i].Beta / (DS->Ele[i].Beta - 1)) - 1, 1 / DS->Ele[i].Beta) / DS->Ele[i].Alpha);
			TotalY = ((Avg_Y<MINpsi) ? MINpsi : Avg_Y) + DS->Ele[i].zmin + DS->DummyY[i_plus_2ele] + (DS->DummyY[i_plus_2ele]>(AquiferDepth - DS->Ele[i].infD) ? 0 : 0.5*(AquiferDepth - DS->Ele[i].infD - DS->DummyY[i_plus_2ele]));
			/* Calculation of gradient between intermediate and bottom layer */
			Grad_Y = (TotalY_unsat - TotalY) / (0.5*DS->Ele[i].infD);
			tmpVar = (DS->Ele[i].Macropore == 1) ? effKV_var(elemSatn, Grad_Y, DS->Ele[i].macKsatV, DS->Ele[i].KsatV, DS->Ele[i].hAreaF) : DS->Ele[i].KsatV*elemSatn;
			//		tmpVar=(DS->Ele[i].Macropore==1)?effKV(satKfunc,Grad_Y,DS->Ele[i].macKsatV,DS->Ele[i].KsatV,DS->Ele[i].hAreaF):DS->Ele[i].KsatV*satKfunc;
			effK = (DS->Ele[i].Macropore == 1) ? ((DS->DummyY[i_plus_2ele] > AquiferDepth - DS->Ele[i].macD) ? tmpVar : DS->Ele[i].KsatV*elemSatn) : DS->Ele[i].KsatV*elemSatn;
			//               effK=(DS->Ele[i].Macropore==1)?((DS->DummyY[i_plus_2ele]>AquiferDepth-DS->Ele[i].macD)?tmpVar:DS->Ele[i].KsatV*satKfunc):DS->Ele[i].KsatV*satKfunc;
			/* Calculation of exchange flux between top and bottom layer */
			/* Using upwind method */
			DS->Recharge[i] = (Grad_Y > 0) ? ((DS->DummyY[i_plus_ele] > EPS_3) ? effK_unsat*Grad_Y : 0) : ((DS->DummyY[i_plus_2ele] > EPS_3) ? effK*Grad_Y : 0);
			DS->RechargeI[i] = DS->Recharge[i];
#endif
		}
#endif
#endif
	}
