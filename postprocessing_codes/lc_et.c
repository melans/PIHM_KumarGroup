// Written by Junehyeong Park, Dec 20 2018
// Most of codes were from 'classSpat_modi.c' for calculating ViR for each land cover type
// Calcualte ET for each land cover type
#include <stdio.h>
#define NumEle 969 // Num. of elements 
#define NumTS 8760 // Num. of records in time-series (model duration 30 (days) * 24 (hours) * 2 (records per hour) 
#define NumLC 14 // Num. of classes of land cover in .forc or .att
#define TSforDay 24 // record TS is 30 minutes, so 48 TS is needed to record daily information
#define Interval 60 // record TS is 30 minutes
main()
	{
	float sumArea[NumLC]={0},val0,val1,val2,et0,et1,et2,avgVar0[NumLC]={0},avgVar1[NumLC]={0},avgVar2[NumLC]={0},area[NumEle];
	int i,j,k,counter=0,lc[NumEle],tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11,tmp12,tmp13,tmp14,tmp15,tmp16,tmp17,tmp18,tmp19,tmp20;
	FILE *fpet0,*fpet1,*fpet2,*fparea,*fpClass,*fpres0,*fpres1,*fpres2,*fpres3;
	fpet0=fopen("test02.et0","r");
	fpet1=fopen("test02.et1","r");
	fpet2=fopen("test02.et2","r");
	fparea=fopen("/home/jpark102/pihm/test02/result/ancillarypihmcodespostprocessing/back/area.txt","r");
	fpClass=fopen("test02.att","r");
	fpres0=fopen("lc_et0.txt","w");
	fpres1=fopen("lc_et1.txt","w");
	fpres2=fopen("lc_et2.txt","w");
	fpres3=fopen("lc_etsum.txt","w");
	for(i=0;i<NumEle;i++)
		{
		fscanf(fparea,"%f",&area[i]);
		fscanf(fpClass,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",&tmp0,&tmp1,&tmp2,&lc[i],&tmp3,&tmp4,&tmp5,&tmp6,&tmp7,&tmp8,&tmp9,&tmp10,&tmp11,&tmp12,&tmp13,&tmp14,&tmp15,&tmp16,&tmp17,&tmp18,&tmp19,&tmp20);
		sumArea[lc[i]-1]=sumArea[lc[i]-1]+area[i]/1000000.0;
		}
	for(i=1;i<=NumTS;i++)
		{
		fscanf(fpet0,"%f",&val0);
		fscanf(fpet1,"%f",&val1);
		fscanf(fpet2,"%f",&val2);
		//printf("%d\n",i);
		if(i%TSforDay==0)
			{
			counter=0;
			fprintf(fpres0,"\n");
			fprintf(fpres1,"\n");
			fprintf(fpres2,"\n");
			fprintf(fpres3,"\n");
			}
		else
			{
			counter++;
			}
		for(j=0;j<NumEle;j++)
			{
			fscanf(fpet0,"%f",&et0);
			fscanf(fpet1,"%f",&et1);
			fscanf(fpet2,"%f",&et2);
			avgVar0[lc[j]-1]+=et0*Interval/60.0/24.0*area[j]/1000000.0;
			avgVar1[lc[j]-1]+=et1*Interval/60.0/24.0*area[j]/1000000.0;
			avgVar2[lc[j]-1]+=et2*Interval/60.0/24.0*area[j]/1000000.0;
			}
		if(counter==0)
			{
			for(k=0;k<NumLC;k++)
				{
				fprintf(fpres0,"%f\t",avgVar0[k]/sumArea[k]);
				fprintf(fpres1,"%f\t",avgVar1[k]/sumArea[k]);
				fprintf(fpres2,"%f\t",avgVar2[k]/sumArea[k]);
				fprintf(fpres3,"%f\t",avgVar0[k]/sumArea[k]+avgVar1[k]/sumArea[k]+avgVar2[k]/sumArea[k]);
				avgVar0[k]=0;
				avgVar1[k]=0;
				avgVar2[k]=0;
				}
			}
		}
	}

