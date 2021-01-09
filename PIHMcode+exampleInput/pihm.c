/*******************************************************************************
 * File        : pihm.c                                                        *
 * Version     : Mar, 2019 (A.0)                                               *
 * Actively maintained and developed by : Mukesh Kumar (mukesh.kumar@ua.edu)   *
 * Contributors to the PIHM code development: Chris Duffy, Yizhong Qu, Mukesh  *
 * 		Kumar, Gopal Bhatt, Xing Chen, Yanlan Liu, Dongdong Wang       * 
 * This code is free for research purposes only.                               *
 * Please provide relevant references if you use this code in your research    *
 *-----------------------------------------------------------------------------*
 * DEVELOPMENT RELATED REFERENCES:					       *
 * parallel-PIHM:							       *
 * 	a) M. Kumar and C. Duffy, Exploring the Role of Domain Partitioning on *
 * 	Efficiency of Parallel Distributed Hydrologic Model Simulations,       *
 * 	Journal of Hydrogeology and Hydrologic Engineering, 2015.	       *
 * PIHMgis:								       *
 * 	a) G. Bhatt, M. Kumar, and C. Duffy, A tightly coupled GIS and distrib-*
 * 	uted hydrologic modeling framework, Environmental Modeling and Software*
 * 	, 2014.								       *
 *	b) M. Kumar, G. Bhatt and C. Duffy, An Object Oriented Shared Data     *
 *	Model for GIS and Distributed Hydrologic Models, International Journal *
 *	of GIS, pp. 1061-1079, v.24-7, 2010.				       *
 *	c) M. Kumar, G. Bhatt and C. Duffy, An efficient domain decomposition  *
 *	framework for accurate representation of geodata in distributed hydrol-*
 *	ogic models, International Journal of GIS, pp. 1569-1596, v.23-12, 2009* 
 * PIHM2.0:								       *
 *	a) Kumar, M., 2009, "Toward a Hydrologic Modeling System", PhD Thesis, *
 *	https://etda.libraries.psu.edu/files/final_submissions/4388	       *
 * PIHM1.0:								       *
 *	a) Qu, Y., 2005, "An Integrated hydrologic model for multiproces       *
 *	simulation using semi-discrete finite volume approach".PhD Thesis, PSU *
 *	b) Qu, Y. & C. Duffy, 2007, "A semidiscrete finite volume formulation  *
 *	for multiprocess watershed simulation". Water Resources Research       *
 *-----------------------------------------------------------------------------*
 * SOME APPLICATION RELATED REFERENCES FROM KUMAR GROUP AND COLLABORATORS      * 
 *	y) S. Seo, R. Bhowmik, S. Arumugam, G. Mahinthakumar, and M. Kumar, The*
 *	 role of cross-correlation between precipitation and temperature in    *
 *	 basin-scale simulations of hydrologic variables, Journal of Hydrology,*
 *	  2019.								       *
 *	z) S. Seo, K. Mahinthakumar, S. Arumugam, and M. Kumar, Conjunctive    *
 *	Management of Surface Water and Groundwater Resources under Drought    *
 *	Conditions Using a Fully Coupled Hydrological Model, Journal of Water  *
 *	Resources Planning and Management, 2018.			       *
 *	a) D. Wang*, Y. Liu*, and M. Kumar, Using nested discretization for a  * 
 *	detailed yet computationally efficient simulation of local hydrology in*
 *	 a distributed hydrologic model, Nature Scientific Reports, 2018       *
 *	b) M. Zhang*, X. Chen*, M. Kumar, M. Marani, and M. Goralczyk, Hurri-  *
 *	canes and Tropical Storms, a Necessary Evil to Ensure Water Supply?,   *
 *	Hydrological Processes, 2017.					       *
 *	c) Y. Liu* and M. Kumar, Role of meteorological controls on the inter- *
 *	annual variations in wet-period characteristics of wetlands, Water     *
 *	Resour. Res., 2016.						       *
 *	d) S. Seo, T. Sinha, K. Mahinthakumar, S. Arumugam, and M. Kumar,      *
 *	Identification of dominant source of errors in developing streamflow   *
 *	and groundwater projection under near-term climate change, Journal of  *
 *	Geophysical Research-Atmospheres, 2016.				       *
 *	e) X. Chen*, M. Kumar, R. Wang*, A. Winstral, and D.Marks, Assessment  * 
 *	of the timing of daily peak streamflow during the melt season in a     *
 *	snow-dominated watershed, Journal of Hydrometeorology, 2016.           *
 *	f) R. Wang*, M. Kumar, and T. Link, Potential Trends in Snowmelt       * 
 *	Generated Peak Streamflows in a Warming Climate, Geophysical Research  *
 *	Letters, 2016.							       *
 *	i) X. Chen*, M. Kumar, and B. McGlynn, Variations in streamflow respo- *
 *	nse to large hurricane-season storms in a southeastern US watershed,   *
 *	Journal of Hydrometeorology, 2015.				       *
 *	j) M. Kumar, D. Marks, J. Dozier, M. Reba, and A. Winstral, Evaluation *
 *	of distributed hydrologic impacts of temperature-index and energy-based* 
 *	snow model simulations, Advances in Water Resources, 2013.	       *
 *	k) R. Wang*, M. Kumar and D. Marks, Anomalous trend in soil evaporation*
 *	in semi-arid, snow-dominated watersheds, Advances in Water Resources,  *
 *	2013.								       *
 *	l) M. Kumar, R. Wang*, and T. Link, Effects of more extreme precipita- *
 *	tion regimes on maximum seasonal snow water equivalent, Geophysical    * 
 *	Research Letters, 2012. 					       *
 *******************************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

/* SUNDIAL Header Files */
#include "sundials/sundials_types.h"    /* realtype, integertype, booleantype defination */
#include "cvode/cvode.h"                /* CVODE header file                             */
#include "cvode/cvode_spgmr.h"          /* CVSPGMR linear header file                    */
#include "sundials/sundials_dense.h"    /* use generic DENSE linear solver for "small"   */
#include "nvector/nvector_serial.h"     /* contains the definition of type N_Vector      */
#include "sundials/sundials_math.h"     /* contains UnitRoundoff, RSqrt, SQR functions   */
#include "cvode/cvode_dense.h"          /* CVDENSE header file                           */
#include "pihm.h"                       /* Data Model and Variable Declarations          */

/* Function Declarations */
void read_alloc(char *, Model_Data, Control_Data *);	/* Variable definition */
void initialize(char *, Model_Data, Control_Data *, N_Vector);
void is_sm_et(realtype, realtype, Model_Data, N_Vector);	
/* Calculate Recharge & RechargeI and maintain the value in the current time step	*/	//xchen_20170305
void f_decouple(Model_Data, realtype, N_Vector); 
/* Function to calculate right hand side of ODE systems */
int f(realtype, N_Vector, N_Vector, void *);
/* Flux update/recalculation after CVODE convergence*/        //xchen_20170305
void flux_cal(realtype, N_Vector, void *);          
void update(realtype, Model_Data);	 
void PrintData(FILE **,Control_Data *, Model_Data, N_Vector, realtype);

/* Main Function */
int main(int argc, char *argv[])
	{  
  	Model_Data mData;               /* Model Data                */
  	Control_Data cData;             /* Solver Control Data       */
  	N_Vector CV_Y,atol;             /* State Variables Vector    */
  	void *cvode_mem;                /* Model Data Pointer        */
  	int flag;                       /* flag to test return value */
  	FILE *Ofile[25];           	/* Output file     */
	char *ofn[25];
	FILE *iproj,*oproj;		/* Project File */
  	int N;                          /* Problem size              */
  	int i,j,k;                      /* loop index                */
  	realtype t,expl_t;    			/* simulation time           */
  	realtype NextPtr, StepSize;     /* stress period & step size */
  	clock_t start, end_r, end_s;    /* system clock at points    */
  	realtype cputime_r, cputime_s;  /* for duration in realtype  */
	char *filename,*outfilename;
        char tmpLname[25][12]={".NetP",".IS",".surf",".unsat",".GW",".et0",".et1",".et2",".P",".Rech",".unsatT",".RechT",".stage",".subriv"};
	/* Project Input Name */




// modified by https://github.com/melans on 01/09/2021 from the original code 

printf("argc => %s 2 \n",(argc!=2)?"NE":"EQ");
        switch (argc)
        {
        case 1: // if 0 args
            // read input file
            iproj=fopen("projectName.txt","r");
            if(iproj==NULL){
                printf("\t\nUsage ./pihm project_name");
                printf("\t\n         OR              ");
                printf("\t\nUsage ./pihm, and have a file in the current directory named projectName.txt with the project name in it\n");
                exit(0);
            }else{
                filename = (char *)malloc(15*sizeof(char));
                fscanf(iproj,"%s",filename);
            }

            // read output file
            oproj=fopen("outputPath.txt","r");
            if(oproj==NULL){
                printf("\t\nUsage ./pihm outputPath");
                printf("\t\n         OR              ");
                printf("\t\nUsage ./pihm, and have a file in the current directory named outputPath.txt with the project name in it\n");
                exit(0);
            }else{
                outfilename = (char *)malloc(100*sizeof(char));
                fscanf(oproj,"%s",outfilename);
			}

            break;
        case 2: // if 1 arg
            iproj=fopen("projectName.txt","r");
            if(iproj==NULL){    // if no input file then read arg as input
                filename = (char *)malloc(strlen(argv[1])*sizeof(char));
                strcpy(filename,argv[1]);
                oproj=fopen("outputPath.txt","r");
                if(oproj==NULL){
                    printf("\t\nUsage ./pihm outputPath");
                    printf("\t\n         OR              ");
                    printf("\t\nUsage ./pihm, and have a file in the current directory named outputPath.txt with the project name in it\n");
                    exit(0);
                }else{
                    outfilename = (char *)malloc(100*sizeof(char));
                    fscanf(oproj,"%s",outfilename);
                }
            }else{  // if there's input file then read arg as output
                filename = (char *)malloc(15*sizeof(char));
                fscanf(iproj,"%s",filename);
                outfilename = (char *)malloc(strlen(argv[1])*sizeof(char));
                strcpy(outfilename,argv[1]);
                // fscanf(oproj,"%s",outfilename);
            }

            break;
        case 3: // if 2 args
            // arg 1 = input
            filename = (char *)malloc(strlen(argv[1])*sizeof(char));
            strcpy(filename,argv[1]);
            // arg 2 = output
            outfilename = (char *)malloc(strlen(argv[2])*sizeof(char));
            strcpy(outfilename,argv[2]);
            break;
        }
printf("infilename = %s \n",filename);
printf("outfilename = %s \n",outfilename);
        /* Create the output directory if it does not exist*/
        struct stat st = {0};
        if (stat(outfilename, &st) == -1) {
            mkdir(outfilename, 0700);
        }

// end of https://github.com/melans modifications

	
	

  	/* allocate memory for model data structure */
  	mData = (Model_Data)malloc(sizeof *mData);
  
  	printf("\n ...  PIHM is starting ... \n");
	printf("\n Version     : Mar 2019 (A.0)\n");                                               
	printf("Contact Mukesh Kumar mukesh.kumar@ua.edu if you have any questions\n");
	printf("This code is free for research purposes only.\n");                             
	printf("Please provide relevant references if you use this code in your research"); 
 	/* read in 9 input files with "filename" as prefix */
  	read_alloc(filename, mData, &cData); 
#ifdef SUB_SURF_RIV
	#ifdef LAYER2
		mData->totele=3*mData->NumEle;
  		atol = N_VNew_Serial(mData->totele + 2*mData->NumRiv);
	#else
		mData->totele=4*mData->NumEle;
  		atol = N_VNew_Serial(mData->totele + 2*mData->NumRiv);
		for(j=0;j<mData->NumEle;j++)
			{
			NV_Ith_S(atol,j+3*mData->NumEle)=cData.abstol[1];
			}
	#endif 	
	for(j=0;j<mData->NumEle;j++)
		{
		NV_Ith_S(atol,j)=cData.abstol[0];
		NV_Ith_S(atol,j+mData->NumEle)=cData.abstol[1];
		NV_Ith_S(atol,j+2*mData->NumEle)=cData.abstol[2];
		}
	for(j=0;j<mData->NumRiv;j++)
		{
		NV_Ith_S(atol,j+mData->totele)=cData.abstol[0];
		NV_Ith_S(atol,j+mData->totele+mData->NumRiv)=cData.abstol[2];
		}
  	/* problem size */
  	N = mData->totele + 2*mData->NumRiv;
#elif SURF_RIV
	mData->totele=mData->NumEle;
	atol = N_VNew_Serial(mData->totele+mData->NumRiv);
        for(j=0;j<mData->NumEle;j++)
                {
                NV_Ith_S(atol,j)=cData.abstol[0];
		}	
        for(j=0;j<mData->NumRiv;j++)
                {
                NV_Ith_S(atol,j+mData->totele)=cData.abstol[0];
		}
  	/* problem size */
  	N = mData->totele+mData->NumRiv;
#endif

	mData->DummyY=(realtype *)malloc(N*sizeof(realtype));
  	/* initial state variable depending on machine*/
  	CV_Y = N_VNew_Serial(N);
  
  	/* initialize mode data structure */
  	initialize(filename, mData, &cData, CV_Y);
  	printf("\nSolving ODE system ... \n");
  
  	/* allocate memory for solver */
  	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  	if(cvode_mem == NULL) {printf("CVodeMalloc failed. \n"); return(1);}
  
  	flag = CVodeSetUserData(cvode_mem, mData);
  	flag = CVodeSetInitStep(cvode_mem,cData.InitStep);
  	flag = CVodeSetStabLimDet(cvode_mem,TRUE);  
  	flag = CVodeSetMaxStep(cvode_mem,cData.MaxStep); 
	flag = CVodeSetMaxNumSteps(cvode_mem, 4000);
	flag = CVodeInit(cvode_mem, f, cData.StartTime, CV_Y);
        flag = CVodeSVtolerances(cvode_mem, cData.reltol, atol);
	flag = CVSpgmr(cvode_mem, PREC_NONE, 0);
  	flag = CVSpilsSetGSType(cvode_mem, MODIFIED_GS);
  
  	/* set start time */
  	t=expl_t = cData.StartTime;
  	start = clock();

        for(i=0;i<11;i++)
                {
                sprintf(tmpLname[i+cData.totFiles-11],".rivFlx%d",i);
                }
        /* Open Output Files */
        for(i=0;i<cData.totFiles;i++)
                {
                ofn[i] = (char *)malloc((strlen(outfilename)+strlen(filename)+12)*sizeof(char));
                strcat(ofn[i], outfilename);
		strcat(ofn[i], filename);
                Ofile[i]=fopen(strcat(ofn[i], tmpLname[i]),"w");
		//printf("%s",strcat(ofn[i], tmpLname[i]));
                //Ofile[i]=fopen(strcat(outfilename,strcat(ofn[i], tmpLname[i])),"w");
                }
  	/* start solver in loops */
  	for(i=0; i<cData.NumSteps; i++)
  		{
    		/* inner loops to next output points with ET step size control */
    		while(expl_t < cData.Tout[i+1])
    			{
      			if (expl_t + cData.ETStep >= cData.Tout[i+1])
      				{
        			NextPtr = cData.Tout[i+1];
      				}
      			else
      				{
        			NextPtr = expl_t + cData.ETStep;
      				}
      			StepSize = NextPtr - expl_t; 
     			 
      			/* calculate Interception Storage */
			is_sm_et(expl_t, StepSize, mData,CV_Y);
			expl_t=expl_t+StepSize;
    			}
		NextPtr=cData.Tout[i+1];
		printf("\n Tsteps = %f ",t); 
     		
	#ifdef DCPL_VFLUX
                f_decouple(mData,t,CV_Y);
        #endif

		flag = CVode(cvode_mem, NextPtr, CV_Y, &t, CV_NORMAL);  
                while(flag==CV_TOO_MUCH_WORK)
                	{
                        flag = CVode(cvode_mem, NextPtr, CV_Y, &t, CV_NORMAL);
                        printf("\nCV TOO MUCH WORK");
			}
		flux_cal(t,CV_Y,mData);                      //xchen 20170305update flux after calculating Y
		PrintData(Ofile,&cData,mData, CV_Y,t);
		update(t,mData);
		}
        /* Free memory */
  	N_VDestroy_Serial(CV_Y);
  	/* Free integrator memory */
  	CVodeFree(cvode_mem);
  	free(mData);
  	return 0;
	}

