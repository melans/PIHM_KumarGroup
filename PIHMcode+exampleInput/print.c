/*******************************************************************************
 * File        : print.c	                                               *
 * Function    : print out model results output files                          *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 * Actively maintained and developed by : Mukesh Kumar (mukesh.kumar@ua.edu)   *
 * This code is free for research purposes only.                               *
 * Please provide relevant references if you use this code in your research    *
 *******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cvode/cvode.h"
#include "cvode/cvode_dense.h"
#include "pihm.h"
/*Temporal average of State vectors */
void avgResults_NV(int fileCounter,FILE **fpin,realtype *varToPrint,N_Vector tmpNV,Model_Data tmpDS,int tmpIntv,realtype tmpt,int tmpCounter,int modIntv)
        {
        int j,tmpNumObj;
	tmpNumObj=(fileCounter+tmpCounter)<12?tmpDS->NumEle:tmpDS->NumRiv;
	switch(fileCounter+tmpCounter)
		{
		case 0:
			for(j=0;j<tmpNumObj;j++) 
				{
				varToPrint[j]=varToPrint[j]+tmpDS->EleNetPrep[j];
				}
                        break;
		case 1:
			for(j=0;j<tmpNumObj;j++) 
				{
				varToPrint[j]=varToPrint[j]+tmpDS->EleIS[j];
				}
                        break;
		case 2:
	#ifdef SUB_SURF_RIV
		case 3:
		case 4:
	#endif
			for(j=0;j<tmpNumObj;j++) 
				{
				varToPrint[j]=varToPrint[j]+NV_Ith_S(tmpNV,j+(fileCounter-2)*tmpDS->NumEle);
				}
                        break;
#ifdef SUB_SURF_RIV
		case 5:
			for(j=0;j<tmpNumObj;j++) 
				{
				varToPrint[j]=varToPrint[j]+tmpDS->EleET[j][0];
				}
                        break;
		case 6:
			for(j=0;j<tmpNumObj;j++) 
				{
				varToPrint[j]=varToPrint[j]+tmpDS->EleET[j][1];
				}
                        break;
		case 7:
			for(j=0;j<tmpNumObj;j++) 
				{
				varToPrint[j]=varToPrint[j]+tmpDS->EleET[j][2];
				}
                        break;
#endif
		case 8:
			for(j=0;j<tmpNumObj;j++) 
				{
                                varToPrint[j]=varToPrint[j]+tmpDS->ElePrep[j];
				}
                        break;
#ifdef SUB_SURF_RIV
		case 9:
			for(j=0;j<tmpNumObj;j++) 
				{
				varToPrint[j]=varToPrint[j]+tmpDS->Recharge[j];
				}
                        break;
	#ifdef LAYER3
		case 10:	
			for(j=0;j<tmpNumObj;j++) 
				{
				varToPrint[j]=varToPrint[j]+NV_Ith_S(tmpNV,j+3*tmpDS->NumEle);
				}
                        break;
		case 11:
			for(j=0;j<tmpNumObj;j++) 
				{
				varToPrint[j]=varToPrint[j]+tmpDS->RechargeI[j];
				}
                        break;
	#endif
#endif
		case 12:
			tmpNumObj=tmpDS->NumRiv;
			for(j=0;j<tmpNumObj;j++) 
				{
				varToPrint[j]=varToPrint[j]+NV_Ith_S(tmpNV,j+tmpDS->totele);
				}
                        break;
	#ifdef SUB_SURF_RIV
		case 13:
			tmpNumObj=tmpDS->NumRiv;
			for(j=0;j<tmpNumObj;j++) 
				{
				varToPrint[j]=varToPrint[j]+NV_Ith_S(tmpNV,j+tmpDS->totele+tmpDS->NumRiv);
				}
                        break;
	#endif
                case 14:
		case 15:
			tmpNumObj=tmpDS->NumRiv;
/*	               	for(j=0;j<tmpNumObj;j++)
        		       {
		               tmpDS->FluxRiv[j][0]=0;
              		       }	
*/			for(j=0;j<tmpNumObj;j++) 
				{
				varToPrint[j]=varToPrint[j]+tmpDS->FluxRiv[j][fileCounter+tmpCounter-14];
				}
                        break;
		case 16:
			tmpNumObj=tmpDS->NumRiv;
			for(j=0;j<tmpNumObj;j++) 
				{
				varToPrint[j]=varToPrint[j]+tmpDS->FluxRiv[j][fileCounter+tmpCounter-14];
				}
                        break;
		case 17:
			tmpNumObj=tmpDS->NumRiv;
			for(j=0;j<tmpNumObj;j++) 
				{
				varToPrint[j]=varToPrint[j]+tmpDS->FluxRiv[j][fileCounter+tmpCounter-14];
				}
                        break;
	#ifdef SUB_SURF_RIV
		case 18:
			tmpNumObj=tmpDS->NumRiv;
			for(j=0;j<tmpNumObj;j++) 
				{
				varToPrint[j]=varToPrint[j]+tmpDS->FluxRiv[j][fileCounter+tmpCounter-14];
				}
                        break;
		case 19:
			tmpNumObj=tmpDS->NumRiv;
			for(j=0;j<tmpNumObj;j++) 
				{
				varToPrint[j]=varToPrint[j]+tmpDS->FluxRiv[j][fileCounter+tmpCounter-14];
				}
                        break;
		case 20:
			tmpNumObj=tmpDS->NumRiv;
			for(j=0;j<tmpNumObj;j++) 
				{
				varToPrint[j]=varToPrint[j]+tmpDS->FluxRiv[j][fileCounter+tmpCounter-14];
				}
                        break;
		case 21:
			tmpNumObj=tmpDS->NumRiv;
			for(j=0;j<tmpNumObj;j++) 
				{
				varToPrint[j]=varToPrint[j]+tmpDS->FluxRiv[j][fileCounter+tmpCounter-14];
				}
                        break;
		case 22:
			tmpNumObj=tmpDS->NumRiv;
			for(j=0;j<tmpNumObj;j++) 
				{
				varToPrint[j]=varToPrint[j]+tmpDS->FluxRiv[j][fileCounter+tmpCounter-14];
				}
                        break;
		case 23:
		case 24:
		       tmpNumObj=tmpDS->NumRiv;
                       for(j=0;j<tmpNumObj;j++)
                               	{
				tmpDS->FluxRiv[j][10]=0;
                               	}

                 	tmpNumObj=tmpDS->NumRiv;
			for(j=0;j<tmpNumObj;j++) 
				{
				varToPrint[j]=varToPrint[j]+tmpDS->FluxRiv[j][fileCounter+tmpCounter-14];
				}
                        break;
	#endif
                default:
                        break;
        	}	
       	if(((int)tmpt%tmpIntv)==0)      
               	{
                fprintf(fpin[fileCounter],"%lf\t",tmpt);
       	        for(j=0;j<tmpNumObj;j++)
               	        {               
                       	fprintf(fpin[fileCounter],"%lf\t",varToPrint[j]/(tmpIntv/modIntv));
                        varToPrint[j]=0; 
       	                }
               	fprintf(fpin[fileCounter],"\n");     
                fflush(fpin[fileCounter]);           
       	        } 
	} 
/* print individual states */
void PrintData(FILE **outp,Control_Data *cD, Model_Data DS, N_Vector CV_Y, realtype t)
	{
	int k,m=0;
	for(k=0;k<cD->NumFilesToPrint;k++)
		{
		avgResults_NV(cD->FileNoToPrint[k],outp,DS->PrintVar[cD->FileNoToPrint[k]],CV_Y,DS,cD->fileInt[cD->FileNoToPrint[k]],t,m,cD->b);
		}
	}
  
