// Made by Junehyeong Park, Dec 19 2018
// 80% of code was copied from 'analysisbalance.c' that Dr. Kumar Mukesh gave it to me
// Org. code was for calculating representative values of whole domain based on results
// I added precipitation and simple water budget calculation

#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
//#include <string.h>
#define NumEle 969 // Num. of elements
#define NumTS 8760 // Num. of records in time-series (model duration 30 (days) * 24 (hours) * 2 (records per hour) 
#define Interval 60 // Time interval length in minute recorded in result files, which is indicated in .para file
#define NumRiv 128 // Num. of river segments
#define RivOutlet 1 // River id of the outlet river segment
main()
	{
	double sumArea=0.0,val,totVar0=0,totVar1=0,totVar2=0,totVar3=0,totVar4=0,totVar5=0,totVar6=0,avget0=0,avget1=0,avget2=0,totet0=0,totet1=0,totet2=0,totetsum=0,avgprec=0,totprec=0,wb=0,area[NumEle],rivarea[NumRiv],val0[NumTS][NumRiv],val1,val2,val3,val4,val5,val6,et0,et1,et2,prec,poros[NumEle][3],init[NumEle][6],porosriv[NumRiv],initriv[NumRiv][2],tvar;
	int i,j;
	FILE *fpin0,*fpin1,*fpin2,*fpin3,*fpin4,*fpin5,*fpin6,*fparea,*fprivarea,*fpresult,*fpporos,*fpinit,*fpet0,*fpet1,*fpet2,*fpprec;

	// reading model result files
	fpin0=fopen("test02.rivFlx1","r");
	fpin1=fopen("test02.unsat","r");
	fpin2=fopen("test02.unsatT","r");
	fpin3=fopen("test02.GW","r");
        fpin4=fopen("test02.stage","r");
        fpin5=fopen("test02.subriv","r");
        fpin6=fopen("test02.surf","r");
	fpet0=fopen("test02.et0","r");
        fpet1=fopen("test02.et1","r");
        fpet2=fopen("test02.et2","r");
	fpprec=fopen("test02.P","r");

	// reading model input file
        fpinit=fopen("test02.init","r");

	// reading files including initial information of the model (can be obtained by modifying 'initialize.c')
	fparea=fopen("/home/jpark102/pihm/test02/result/ancillarypihmcodespostprocessing/back/area.txt","r");
        fprivarea=fopen("/home/jpark102/pihm/test02/result/ancillarypihmcodespostprocessing/back/rivarea.txt","r");
	fpporos=fopen("/home/jpark102/pihm/test02/result/ancillarypihmcodespostprocessing/back/poros.txt","r");

	// opening the balance analysis result file
	fpresult=fopen("test02.wb","w");

       for(i=0;i<NumEle;i++)
          {
                fscanf(fparea,"%lf",&area[i]);
	        fscanf(fpporos,"%lf%lf%lf",&poros[i][0],&poros[i][1],&poros[i][2]);
	        fscanf(fpinit,"%lf%lf%lf%lf%lf%lf",&init[i][0],&init[i][1],&init[i][2],&init[i][3],&init[i][4],&init[i][5]);
         	sumArea+=area[i];
	        //poros[i][0]=0.09647;
		//poros[i][1]=0.09647;
		//poros[i][2]=0.09647;
	   }
       for(i=0;i<NumRiv;i++)
          {     fscanf(fprivarea,"%lf",&rivarea[i]); 
               // fscanf(fpporos,"%lf",&porosriv[i]);
               fscanf(fpinit,"%lf%lf",&initriv[i][0],&initriv[i][1]);
                //porosriv[i]=0.09647;
                porosriv[i]=0.44;
          }
      fprintf(fpresult,"rivFlx1 unsat unsatT GW stage subriv surf et_sum et0 et1 et2 prec wb\n");
       for(i=0;i<NumTS;i++)
	   {
		fscanf(fpin0,"%lf",&val);
                fscanf(fpin4,"%lf",&val);
                fscanf(fpin5,"%lf",&val);
		
                totVar4=totVar5=0.0;  
                for(j=0;j<NumRiv;j++)
                       {
                       fscanf(fpin0,"%lf",&val0[i][j]);
                       fscanf(fpin4,"%lf",&val4);
                       fscanf(fpin5,"%lf",&val5);
                       
                       //val0[i][j]=(val4<=0)?0:val0[i][j]; 
                       //val4=(val4<0)?0:val4;
                       //val5=(val5<0)?0:val5;
//                       totVar0+=val0*Interval/1440/sumArea;
                       totVar4+=(val4-initriv[j][0])*rivarea[j]/sumArea;
                       totVar5+=(val5-initriv[j][1])*porosriv[j]*rivarea[j]/sumArea;
                       }
		
                fscanf(fpin1,"%lf",&val);
		fscanf(fpin2,"%lf",&val);
		fscanf(fpin3,"%lf",&val);
                fscanf(fpin6,"%lf",&val);
		fscanf(fpet0,"%lf",&val);
                fscanf(fpet1,"%lf",&val);
                fscanf(fpet2,"%lf",&val);
		fscanf(fpprec,"%lf",&val);

                totVar1=totVar2=totVar3=totVar6=avget0=avget1=avget2=avgprec=0.0;
		for(j=0;j<NumEle;j++)
			{
			fscanf(fpin1,"%lf",&val1);
			fscanf(fpin2,"%lf",&val2);
			fscanf(fpin3,"%lf",&val3);
                        fscanf(fpin6,"%lf",&val6);
			fscanf(fpet0,"%lf",&et0);
                        fscanf(fpet1,"%lf",&et1);
                        fscanf(fpet2,"%lf",&et2);
			fscanf(fpprec,"%lf",&prec);
			if(j==0)
				{
				tvar=prec;
				}
//			printf("%1f\t",prec);

			//val1=(val1<0)?0:val1;
			//val2=(val2<0)?0:val2;
			//val3=(val3<0)?0:val3;
                        //val6=(val6<0)?0:val6;

		        totVar1+=(val1-init[j][3])*poros[j][0]*area[j]/sumArea; 
	         	totVar2+=(val2-init[j][5])*poros[j][1]*area[j]/sumArea;
		        totVar3+=(val3-init[j][4])*poros[j][2]*area[j]/sumArea;
                        totVar6+=(val6-init[j][2])*area[j]/sumArea;
			avgprec+=tvar*Interval/1440*area[j]/sumArea;
//			avgprec+=prec*area[j];
			avget0+=et0*Interval/1440*area[j]/sumArea;
                        avget1+=et1*Interval/1440*area[j]/sumArea;
                        avget2+=et2*Interval/1440*area[j]/sumArea;

                        }
                 
		totVar0+=val0[i][RivOutlet-1]*Interval/1440/sumArea;
		totprec+=avgprec;
		totet0+=avget0;
		totet1+=avget1;
		totet2+=avget2;
		totetsum=totet0+totet1+totet2;
		wb=totprec-totVar0-totVar1-totVar2-totVar3-totVar4-totVar5-totVar6-totet0-totet1-totet2;
//		wb=totprec-totVar0-totVar4-totVar6;
		fprintf(fpresult,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",totVar0,totVar1,totVar2,totVar3,totVar4,totVar5,totVar6,totetsum,totet0,totet1,totet2,totprec,wb);
	}
	fclose(fpin0);
	fclose(fpin1);
	fclose(fpin2);
	fclose(fpin3);
	fclose(fpin4);
	fclose(fpin5);
	fclose(fpin6);
	fclose(fpporos);
	fclose(fpresult);
	fclose(fparea);
	fclose(fprivarea);
	fclose(fpinit);	
	fclose(fpet1);
	fclose(fpet2);
	fclose(fpprec);
	//printf("Water budget analysis is done.\n");
	}
