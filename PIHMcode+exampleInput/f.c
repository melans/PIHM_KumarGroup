/*******************************************************************************
 *-----------------------------------------------------------------------------*
 * File        : f.c   (PIHM v.A.0)                                            *
 * Function    : Model Kernel: Building ODE system for each physical process   *
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
#include "nvector/nvector_serial.h"
#include "sundials/sundials_types.h"

#include "f_functions.h"

int f(realtype t, N_Vector CV_Y, N_Vector CV_Ydot, void *DS)
	{
  	int i, j,ieleBC,inabr,inabr_left,inabr_right,inabr_plus_2ele,ileft,iright,i_plus_ele,i_plus_2ele,i_plus_3ele,iright_plus_2ele, ileft_plus_2ele, i_plus_totele, i_plus_totele1riv;
  	realtype Delta, Gamma,tmpVar,tmpVar1;
  	realtype Rn, T, Vel, RH, VP,P,LAI,zero_dh,cnpy_h,rl,r_a,r_s,alpha_r,f_r,eta_s,beta_s,Rmax;
	realtype Avg_Sf,Distance;
   	realtype Cwr,Perem, Perem_down,Avg_Rough,Avg_Perem,Avg_Y_Riv,Wid,Wid_down,Avg_Wid;
   	realtype nabrAqDepth,AquiferDepth, Deficit,elemSatn,satKfunc,effKI,effK_unsat,TotalY_Ele,TotalY_Ele_down,TotalY_unsat;
	realtype elemSatn1,elemSatn2,beta_s1,beta_s2,r_s1,r_s2;
  	realtype *Y, *DY;
  	Model_Data MD;
  	Y = NV_DATA_S(CV_Y);
  	DY = NV_DATA_S(CV_Ydot);
  	MD = (Model_Data) DS;

	flux_cal(t,CV_Y,MD);
	for(i=0; i<MD->totele+MD->totriv; i++)
		{
		DY[i]=0;
		MD->DummyY[i]=(Y[i]>=0)?Y[i]:0;
		}
	for(i=0;i<MD->NumEle;i++)
		{
                i_plus_ele=i+MD->NumEle;
                i_plus_2ele=i+2*MD->NumEle;
                AquiferDepth=(MD->Ele[i].zmax-MD->Ele[i].zmin);
#ifdef SUB_SURF_RIV 
        #ifdef LAYER3
                i_plus_3ele=i+3*MD->NumEle;
		if(MD->DummyY[i_plus_ele]+MD->DummyY[i_plus_2ele]+MD->DummyY[i_plus_3ele]<AquiferDepth)
			{
			DY[i_plus_3ele] = MD->RechargeI[i]-MD->Recharge[i];
			}
		tmpVar=MD->DummyY[i_plus_3ele];
        #else 
		tmpVar=0;
	#endif
                if(MD->DummyY[i_plus_ele]+MD->DummyY[i_plus_2ele]+tmpVar<AquiferDepth)
                	{
                	DY[i] = MD->EleNetPrep[i] - MD->EleViR[i]-((MD->DummyY[i]<EPS_3)?0:MD->EleET[i][2]);
                	DY[i_plus_ele] = MD->EleViR[i]-MD->RechargeI[i]-((MD->DummyY[i]<EPS_3)?MD->EleET[i][2]:0);
                	DY[i_plus_2ele]= MD->Recharge[i];
                	}   
                else
                	{
                	DY[i] =MD->EleNetPrep[i]-((MD->DummyY[i]<EPS_3)?0:MD->EleET[i][2]);
			DY[i_plus_ele] =-((MD->DummyY[i]<EPS_3)?MD->EleET[i][2]:0);
                	}   

/*                DY[i] = MD->EleNetPrep[i] - MD->EleViR[i]-((MD->DummyY[i]<EPS_3)?0:MD->EleET[i][2]);
                DY[i_plus_ele] = MD->EleViR[i]-MD->RechargeI[i]-((MD->DummyY[i]<EPS_3)?MD->EleET[i][2]:0);
                DY[i_plus_2ele]= MD->Recharge[i];
*/                if(MD->DummyY[i_plus_2ele]>AquiferDepth-MD->Ele[i].RzD)
                        {   
                        DY[i_plus_2ele]=DY[i_plus_2ele]-MD->EleET[i][1];
                        }   
                else 
                        {   
                #ifdef LAYER2
                        DY[i_plus_ele] = DY[i_plus_ele]-MD->EleET[i][1];
                #endif 
/*      xc 20130529     The following block is to incorporate ET1 calculation from second unsaturated zone. 
 *                              Purpose is to extract ET1 from the 2 unsaturated zones based on water availability.                             */        
		#ifdef LAYER3
                        elemSatn1 =MD->DummyY[i_plus_ele]/MD->Ele[i].infD;
                        elemSatn2 =MD->DummyY[i_plus_3ele]/(((AquiferDepth-MD->Ele[i].infD-MD->DummyY[i_plus_2ele])<MD->Ele[i].infD)?MD->Ele[i].infD:(AquiferDepth-MD->Ele[i].infD-MD->DummyY[i_plus_2ele]));
                        beta_s1= ((elemSatn1<EPS_4)?EPS_4:((elemSatn1>1.0)?1:elemSatn1));
                        beta_s2= ((elemSatn2<EPS_4)?EPS_4:((elemSatn2>1.0)?1:elemSatn2));
                        s1 = beta_s1;
                        s2 = beta_s2;

                        if(s1<=EPS_4)
                                {   
                                if (s2<=EPS_4)
                                        {   
                                        }   
                                else
                                        {   
                                        DY[i_plus_3ele]=DY[i_plus_3ele]-MD->EleET[i][1];
                                        }   
                                }   
                        else
                                {   
                                if (s2<=EPS_4)
                                        {   
                                        DY[i_plus_ele] = DY[i_plus_ele]-MD->EleET[i][1];
                                        }   
                                else
                                        {   
                                        tmpVar1=s1/(s1+s2);
                                        DY[i_plus_ele] = DY[i_plus_ele]-tmpVar1*MD->EleET[i][1];
                                        DY[i_plus_3ele] = DY[i_plus_3ele]-(1-tmpVar1)*MD->EleET[i][1];
                                        }   
                                }   
                #endif
                        }  
#elif SURF_RIV
                DY[i] = MD->EleNetPrep[i];
#endif
		tmpVar1=0;
     		for(j=0; j<3; j++)
      			{
        		DY[i] =  DY[i] - MD->FluxSurf[i][j]/MD->Ele[i].area;
		#ifdef SUB_SURF_RIV
			tmpVar1=tmpVar1+MD->FluxSub[i][j]/MD->Ele[i].area;
		#endif
      			}
	#ifdef SUB_SURF_RIV
                if(MD->DummyY[i_plus_ele]+MD->DummyY[i_plus_2ele]+tmpVar>=AquiferDepth)
                	{
			if(tmpVar1<0)
				{
                		DY[i] =DY[i]-tmpVar1;
				#ifdef LAYER3
				if(MD->DummyY[i_plus_3ele]>MD->Ele[i].infD)
					{
					DY[i_plus_3ele]= DY[i_plus_3ele]+tmpVar1;
        	        		DY[i_plus_2ele]= DY[i_plus_2ele]-tmpVar1;
					}
				else if(MD->DummyY[i_plus_ele]>MD->Ele[i].infD)
					{
                                        DY[i_plus_ele]= DY[i_plus_ele]+tmpVar1;
                                        DY[i_plus_2ele]= DY[i_plus_2ele]-tmpVar1;
					}
				#elif LAYER2	
				if(MD->DummyY[i_plus_ele]>MD->Ele[i].infD)
					{
                                        DY[i_plus_ele]= DY[i_plus_ele]+tmpVar1;
                                        DY[i_plus_2ele]= DY[i_plus_2ele]-tmpVar1;
					}
				#endif
                		}
                	else
                		{
                		DY[i_plus_2ele]= DY[i_plus_2ele]-tmpVar1;
                		}
			}
		else
			{
			DY[i_plus_2ele]= DY[i_plus_2ele]-tmpVar1;
			}
        	MD->FluxSource[i] =(MD->Ele[i].source>0)?Interpolation(&MD->TSD_Source[MD->Ele[i].source - 1], t):0;           
		DY[i_plus_2ele] = DY[i_plus_2ele] - ((MD->DummyY[i_plus_2ele]>0)?MD->FluxSource[i]/MD->Ele[i].area:0);	
	#endif
		DY[i]=DY[i]/(UNIT_C);
	#ifdef NO_UNSAT
		DY[i_plus_2ele] = DY[i_plus_2ele] + MD->ElePrep[i];
	#endif
#ifdef SUB_SURF_RIV
	      	DY[i_plus_ele] = DY[i_plus_ele]/(MD->Ele[i].Porosity*UNIT_C);
      		DY[i_plus_2ele] = DY[i_plus_2ele]/(MD->Ele[i].Porosity*UNIT_C);
	#ifdef LAYER3
        	i_plus_3ele=i+3*MD->NumEle;
      		DY[i_plus_3ele] = DY[i_plus_3ele]/(MD->Ele[i].Porosity*UNIT_C);
	#endif
#endif
		}
   	for(i=0; i<MD->NumRiv; i++)
    		{
        	i_plus_ele=i+MD->NumEle;
        	i_plus_totele=i+MD->totele;
        	i_plus_totele1riv=i+MD->totele+MD->NumRiv;
		for(j=0;j<=6;j++)
			{
			/* Note the limitation due to d(v)/dt=a*dy/dt+y*da/dt for CS other than rectangle */
			DY[i_plus_totele] = DY[i_plus_totele]-MD->FluxRiv[i][j]/(MD->Riv[i].Length*CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,MD->Riv[i].depth,MD->Riv[i].coeff,3));
			}
		DY[i_plus_totele] = DY[i_plus_totele]/(UNIT_C);
	#ifdef SUB_SURF_RIV
		DY[i_plus_totele1riv] = DY[i_plus_totele1riv] -MD->FluxRiv[i][7] -MD->FluxRiv[i][8]-MD->FluxRiv[i][9] -MD->FluxRiv[i][10]+MD->FluxRiv[i][6];
      		DY[i_plus_totele1riv] = DY[i_plus_totele1riv]/(MD->Ele[i_plus_ele].Porosity*MD->Riv[i].Length*CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,MD->Riv[i].depth,MD->Riv[i].coeff,3)*UNIT_C);
	#endif
    		}
//	printf("\nf  %lf, %lf, %lf, %lf, %lf, %lf",MD->FluxRiv[0][0],MD->FluxRiv[1][0],MD->FluxRiv[0][1],MD->FluxRiv[1][1],MD->FluxRiv[1][2],MD->FluxRiv[1][3]);
	return(0);
	}

