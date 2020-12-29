/*******************************************************************************
 *-----------------------------------------------------------------------------*
 * File        : flux_cal.c   (PIHM v.A.0)                                     *
 * Function    : Update flux after CVODE convergence                           *
 *-----------------------------------------------------------------------------*
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

void flux_cal(realtype t,N_Vector CV_Y,void *DS)
	{
  	int i, j,ieleBC,inabr,inabr_left,inabr_right,inabr_plus_2ele,ileft,iright,i_plus_ele,i_plus_2ele,i_plus_3ele,iright_plus_2ele, ileft_plus_2ele, i_plus_totele, i_plus_totele1riv;
  	realtype Delta, Gamma,tmpVar,tmpVar1;
  	realtype Rn, T, Vel, RH, VP,P,LAI,zero_dh,cnpy_h,rl,r_a,r_s,alpha_r,f_r,eta_s,beta_s,Rmax;
	realtype Avg_Sf,Distance;
   	realtype Cwr,Perem, Perem_down,Avg_Rough,Avg_Perem,Avg_Y_Riv,Wid,Wid_down,Avg_Wid;
   	realtype nabrAqDepth,AquiferDepth, Deficit,elemSatn,satKfunc,effKI,effK_unsat,TotalY_Ele,TotalY_Ele_down,TotalY_unsat;
  	realtype *Y;
  	Model_Data MD;
  	Y = NV_DATA_S(CV_Y);
  	MD = (Model_Data) DS;

	/* Initialization of temporary state variables */
       	MD->totriv=2*MD->NumRiv;
        #ifdef SURF_RIV
        MD->totriv=MD->NumRiv;
        #endif
	for(i=0; i<MD->totele+MD->totriv; i++)
  		{
		MD->DummyY[i]=(Y[i]>=0)?Y[i]:0;
                if(i<MD->NumRiv)
                        {   
                        MD->FluxRiv[i][0]=0;
                        MD->FluxRiv[i][10]=0;
                        }   
   	#ifdef DIFFUSION
		if(i<MD->NumEle)
			{
			for(j=0;j<3;j++)
				{
                		inabr=MD->Ele[i].nabr[j]-1;
                		ieleBC=-(MD->Ele[i].BC[j]/4)-1;
			      	MD->Ele[i].surfH[j]=(inabr>-1)?((-ieleBC>0)?(MD->Ele[inabr].zmax+MD->DummyY[inabr]):((MD->DummyY[ieleBC+MD->totele]>MD->Riv[ieleBC].depth)?MD->Riv[ieleBC].zmin+MD->DummyY[ieleBC+MD->totele]:MD->Riv[ieleBC].zmax)):((MD->Ele[i].BC[j]!=1)?(MD->Ele[i].zmax+MD->DummyY[i]):Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j])-1], t));
			   //  MD->Ele[i].surfH[j]=0;
                        	}
                	MD->Ele[i].dhBYdx=-1*(MD->Ele[i].surfY[2]*(MD->Ele[i].surfH[1]-MD->Ele[i].surfH[0])+MD->Ele[i].surfY[1]*(MD->Ele[i].surfH[0]-MD->Ele[i].surfH[2])+MD->Ele[i].surfY[0]*(MD->Ele[i].surfH[2]-MD->Ele[i].surfH[1]))/(MD->Ele[i].surfX[2]*(MD->Ele[i].surfY[1]-MD->Ele[i].surfY[0])+MD->Ele[i].surfX[1]*(MD->Ele[i].surfY[0]-MD->Ele[i].surfY[2])+MD->Ele[i].surfX[0]*(MD->Ele[i].surfY[2]-MD->Ele[i].surfY[1]));
                	MD->Ele[i].dhBYdy=-1*(MD->Ele[i].surfX[2]*(MD->Ele[i].surfH[1]-MD->Ele[i].surfH[0])+MD->Ele[i].surfX[1]*(MD->Ele[i].surfH[0]-MD->Ele[i].surfH[2])+MD->Ele[i].surfX[0]*(MD->Ele[i].surfH[2]-MD->Ele[i].surfH[1]))/(MD->Ele[i].surfY[2]*(MD->Ele[i].surfX[1]-MD->Ele[i].surfX[0])+MD->Ele[i].surfY[1]*(MD->Ele[i].surfX[0]-MD->Ele[i].surfX[2])+MD->Ele[i].surfY[0]*(MD->Ele[i].surfX[2]-MD->Ele[i].surfX[1]));
			}
   	#endif
  		}	
	/* Lateral Flux Calculation between Triangular elements Follows  */
	for(i=0; i<MD->NumEle; i++)
  		{
        	i_plus_ele=i+MD->NumEle;
        	i_plus_2ele=i+2*MD->NumEle;
		AquiferDepth=(MD->Ele[i].zmax-MD->Ele[i].zmin);
    		for(j=0; j<3; j++)
    			{
      			if(MD->Ele[i].nabr[j] > 0)
      				{
               			inabr=MD->Ele[i].nabr[j]-1;   
			#ifdef SUB_SURF_RIV
                                if(MD->Ele[i].Calc[j]>0)
                                MD->FluxSub[i][j] = - MD->FluxSub[inabr][MD->Ele[i].Calc[j]-1];
                                else
                                {
                        	inabr_plus_2ele=inabr+2*MD->NumEle;
			        /***************************************************************************/
			        /* Subsurface Lateral Flux Calculation between Triangular elements Follows */
			        /***************************************************************************/
                    		GradCalc(MD->DummyY[inabr_plus_2ele], MD->Ele[inabr].zmin, MD->DummyY[i_plus_2ele], MD->Ele[i].zmin, MD->Ele[i].Dist[j],0,0,0,0,ELE_ELE);
        			/* take care of macropore effect */
                    		nabrAqDepth=(MD->Ele[inabr].zmax-MD->Ele[inabr].zmin);
                    		avgKH(MD->Ele[i].Macropore,MD->DummyY[i_plus_2ele],AquiferDepth,MD->Ele[i].macD,MD->Ele[i].macKsatH,MD->Ele[i].vAreaF,MD->Ele[i].KsatH,MD->Ele[inabr].Macropore,MD->DummyY[inabr_plus_2ele],nabrAqDepth,MD->Ele[inabr].macD,MD->Ele[inabr].macKsatH,MD->Ele[inabr].vAreaF,MD->Ele[inabr].KsatH);
        			/* groundwater flow modeled by Darcy's law */
        			MD->FluxSub[i][j] = Avg_Ksat*Grad_Y*Avg_Y*MD->Ele[i].edge[j];
                               }
			#endif
		#ifndef NO_UNSAT
			    	/***************************************************************************/
			    	/* Surface Lateral Flux Calculation between Triangular elements Follows    */
			    	/***************************************************************************/
                        if(MD->Ele[i].Calc[j]>0)
                                MD->FluxSurf[i][j] = - MD->FluxSurf[inabr][MD->Ele[i].Calc[j]-1];
                         else
                        {
                 	#ifdef DIFFUSION
                    		GradCalc(MD->DummyY[inabr], MD->Ele[inabr].zmax, MD->DummyY[i], MD->Ele[i].zmax, MD->Ele[i].Dist[j],0,0,0,0,ELE_ELE);
                    		Avg_Sf=sqrt(pow(MD->Ele[i].dhBYdx,2)+pow(MD->Ele[i].dhBYdy,2));
                    		Avg_Sf=(Avg_Sf>EPS_5)?Avg_Sf:EPS_5;
                    	#else
                    		GradCalc(MD->DummyY[inabr], MD->Ele[inabr].zmax-MD->DummyY[inabr], MD->DummyY[i], MD->Ele[i].zmax-MD->DummyY[i], MD->Ele[i].Dist[j],0,0,0,0,ELE_ELE);
                    		Avg_Sf=(fabs(Grad_Y)>EPS_5)?Grad_Y:EPS_5;
                    	#endif
                    		/* Weighting needed */
        			Avg_Rough = 0.5*(MD->Ele[i].Rough + MD->Ele[inabr].Rough);
        			CrossA = Avg_Y*MD->Ele[i].edge[j];
                    		/* INCLUDE CROSSADOWN */
				OverlandFlow(MD->FluxSurf,i,j, Avg_Y,Grad_Y,Avg_Sf,CrossA,Avg_Rough);
	                }
            	#endif
      				}
			/************************************************/
			/* Boundary condition Flux Calculations Follows */
			/************************************************/
      			else
      				{
        			/*  No flow (natural) boundary condition is default */
        			if(MD->Ele[i].BC[j] == 0)
        				{
          				MD->FluxSurf[i][j] = 0;
          				MD->FluxSub[i][j] = 0;
        				}
        			else if(MD->Ele[i].BC[j] == 1)	/* Note: ideally different boundary conditions need to be incorporated	for surf and subsurf respectively */
					/* Note: the formulation assumes only dirichlet TS right now */
        				{
          				MD->FluxSurf[i][j] = 0;	/* Note the assumption here is no flow for surface*/ 
                        		effK=effKH(MD->Ele[i].Macropore,MD->DummyY[i_plus_2ele],AquiferDepth,MD->Ele[i].macD,MD->Ele[i].macKsatH,MD->Ele[i].vAreaF,MD->Ele[i].KsatH);
           				Avg_Ksat = effK;
                        		GradCalc(Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j])-1], t)- MD->Ele[i].zmin, MD->Ele[i].zmin, MD->DummyY[i_plus_2ele], MD->Ele[i].zmin, MD->Ele[i].Dist[j],0,0,0,0,ELE_ELE);
           				MD->FluxSub[i][j] = Avg_Ksat*Grad_Y*Avg_Y*MD->Ele[i].edge[j];
          				}
          			else       /* Neumann BC (Note: MD->Ele[i].BC[j] value have to be = 2+(index of neumann boundary TS)*/
          				{
            				MD->FluxSurf[i][j] = Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j])-1], t);
            				MD->FluxSub[i][j] = Interpolation(&MD->TSD_EleBC[(-MD->Ele[i].BC[j])-1], t);
          				}
      				}
    			}
#ifndef NO_UNSAT
	#ifdef SUB_SURF_RIV		
		/* Calculation of ET0, ET1, ViR, Recharge, RechargeI is performed in fluxCalc_Ele() */
		fluxCalc_Ele(MD,i,t); 
	#endif
#endif
  		} 
	/* Lateral Flux Calculation between River-River and River-Triangular elements Follows */   
	for(i=0; i<MD->NumRiv; i++)
  		{
		fluxCalc_Riv(MD,i,t);
  		}
	}

