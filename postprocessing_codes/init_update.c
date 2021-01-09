// Made by Junehyeong Park, Dec 19 2018
// 80% of code was copied from 'analysisbalance.c' that Dr. Kumar Mukesh gave it to me
// Org. code was for calculating representative values of whole domain based on results
// I added precipitation and simple water budget calculation

#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
//#include <string.h>
#define TargetTS 354
#define NumEle 969 // Num. of elements
#define NumRiv 128 // Num. of river segments
main()
	{
	double val,val_unsat,val_unsatT,val_GW,val_stage,val_subriv,init[NumEle][6],init_new[NumEle][6],initriv[NumRiv][2],initriv_new[NumRiv][2];
	int i,j;
	FILE *f_unsat,*f_unsatT,*f_GW,*f_stage,*f_subriv,*f_init,*f_result;

	// reading model result files
	f_unsat=fopen("test02.unsat","r");
	f_unsatT=fopen("test02.unsatT","r");
	f_GW=fopen("test02.GW","r");
	f_stage=fopen("test02.stage","r");
	f_subriv=fopen("test02.subriv","r");

	// reading model input file
	f_init=fopen("test02.init","r");

	// opening result file
	f_result=fopen("test02.init_new","w");

	for(i=0;i<NumEle;i++)
		{
		fscanf(f_init,"%lf%lf%lf%lf%lf%lf",&init[i][0],&init[i][1],&init[i][2],&init[i][3],&init[i][4],&init[i][5]);
		}
	for(i=0;i<NumRiv;i++)
		{     
		fscanf(f_init,"%lf%lf",&initriv[i][0],&initriv[i][1]);
		}

	for(i=0;i<TargetTS-1;i++)
		{
		fscanf(f_unsat,"%lf",&val);
		fscanf(f_unsatT,"%lf",&val);
		fscanf(f_GW,"%lf",&val);

		for(j=0;j<NumEle;j++)
			{
			fscanf(f_unsat,"%lf",&val);
			fscanf(f_unsatT,"%lf",&val);
			fscanf(f_GW,"%lf",&val);
                        }
                 
		fscanf(f_subriv,"%lf",&val);
		fscanf(f_stage,"%lf",&val);

		for(j=0;j<NumRiv;j++)
			{
			fscanf(f_subriv,"%lf",&val);
			fscanf(f_stage,"%lf",&val);
			}
		}

	for(i=TargetTS-1;i<TargetTS;i++)
		{
		fscanf(f_unsat,"%lf",&val);
                fscanf(f_unsatT,"%lf",&val);
                fscanf(f_GW,"%lf",&val);

                for(j=0;j<NumEle;j++)
                        {
                        fscanf(f_unsat,"%lf",&val_unsat);
                        fscanf(f_unsatT,"%lf",&val_unsatT);
                        fscanf(f_GW,"%lf",&val_GW);

			init_new[j][3]=val_unsat;
			init_new[j][4]=val_GW;
			init_new[j][5]=val_unsatT;
                        }

                fscanf(f_subriv,"%lf",&val);
                fscanf(f_stage,"%lf",&val);

                for(j=0;j<NumRiv;j++)
                        {
                        fscanf(f_subriv,"%lf",&val_subriv);
                        fscanf(f_stage,"%lf",&val_stage);

			initriv_new[j][0]=val_stage;
			initriv_new[j][1]=val_subriv;
                        }
		}

	for(i=0;i<NumEle;i++)
		{
		fprintf(f_result,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",init[i][0],init[i][1],init[i][2],init_new[i][3],init_new[i][4],init_new[i][5]);
		}
	for(i=0;i<NumRiv;i++)
		{
		fprintf(f_result,"%lf\t%lf\n",initriv_new[i][0],initriv_new[i][1]);
		}

	fclose(f_unsat);
	fclose(f_unsatT);
	fclose(f_GW);
	fclose(f_subriv);
	fclose(f_stage);
	fclose(f_init);
	fclose(f_result);
	}
