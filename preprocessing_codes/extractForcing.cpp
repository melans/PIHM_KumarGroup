// ******************************************************************
// * EXTRACT NLDAS 2 CLIMATE DATA                                   *
// * Version 1.0 .1                                                 *
// * Developer:GOPAL BHATT(gopal.bhatt @ psu.edu)                   *
// * Contact:Prof.CHRIS DUFFY(cxd11 @ psu.edu)
// * Modifying:Junehyeong Park(sai0259@gmail.com)                   *
// ******************************************************************


#include <fstream>
#include <iomanip>
#include <ctime>
#include <string>

#include "pickGridVal.h"

#include <gdal.h>
#include <gdal_priv.h>

#define DEBUG if(0)

using namespace std;

	typedef struct {
		int             row;
		int             col;
	}               GRIDS;

void            stepOneHour(struct tm * t);

int mainFunction(int syyyy,	int	smm,	int sdd,	int shh,	int eyyyy,	int emm,	int edd,	int ehh,	int NumGrid,	int *col, int *row);
int main(){

	//int col[12]={367, 367, 367, 368, 368, 368, 369, 369, 369, 370, 370, 370};
	//int row[12]={90, 91, 92, 90, 91, 92, 90, 91, 92, 90, 91, 92};

	//mainFunction( 1980,		1,	 1,	 5,	 2012,	 1,	 1,	 5,	 12,	 row,	 col);

/*	int col[6]; 
	int row[6];
	
	row[0]=91;
	col[0]=369;
	row[1]=91;
	col[1]=370;
	row[2]=90;
	col[2]=370;
	row[3]=90;
	col[3]=369;
	row[4]=90;
	col[4]=368;
	row[5]=91;
	col[5]=368;
*/
/*	int col[8]; 
	int row[8];
	
	row[0]=367;
	col[0]=90;
	row[1]=368;
	col[1]=90;
	row[2]=369;
	col[2]=90;
	row[3]=370;
	col[3]=90;
	row[4]=368;
	col[4]=89;
	row[5]=369;
	col[5]=89;
	row[6]=370;
	col[6]=89;
	row[7]=356;
	col[7]=85;

	//mainFunction( 2005,	1,	 1,	0,	 2014,	 9,	 23,	 12,	 8,	 row,	 col);
	mainFunction( 2006,	1,	1,	0,	 2007,	 1,	 1,  23,	 8,	 row,	 col);
*/
	int col[6]; 
	int row[6];
	//int col[1];
        //int row[1];
	
        //row[0]=87;
        //col[0]=356;
	row[0]=354;
	col[0]=85;
	row[1]=354;
	col[1]=86;
        row[2]=355;
        col[2]=85;
        row[3]=355;
        col[3]=86;
        row[4]=356;
        col[4]=85;
        row[5]=356;
        col[5]=86;
        //row[6]=355;
        //col[6]=85;
        //row[7]=356;
        //col[8]=85;

	//mainFunction( 1979, 1, 2, 0, 2017, 12, 31, 23, 8, row,col);
	//mainFunction( 1980, 1, 1, 0, 2017, 12, 31, 23, 6, row, col);
	mainFunction( 2015, 1, 1, 0, 2017, 12, 31, 23, 6, row, col);
	return 0;
}

	void            stepOneHour(struct tm * t)
{
	int             yyyy, mm, dd, dd2, hh;
	yyyy = t->tm_year;
	mm = t->tm_mon + 1;
	dd = t->tm_mday;
	dd2 = t->tm_yday + 1;
	hh = t->tm_hour;

	int             dmax, lpyear = 0;
	if (yyyy % 4 == 0 || yyyy % 100 == 0)
		lpyear = 1;
	switch (mm) {
	case 1:
		dmax = 31;
		break;
	case 2:
		if (lpyear == 1)
			dmax = 29;
		else
			dmax = 28;
		break;
	case 3:
		dmax = 31;
		break;
	case 4:
		dmax = 30;
		break;
	case 5:
		dmax = 31;
		break;
	case 6:
		dmax = 30;
		break;
	case 7:
		dmax = 31;
		break;
	case 8:
		dmax = 31;
		break;
	case 9:
		dmax = 30;
		break;
	case 10:
		dmax = 31;
		break;
	case 11:
		dmax = 30;
		break;
	case 12:
		dmax = 31;
		break;
	}

	if (hh < 23) {
		t->tm_hour++;
		return;
	} else if (dd < dmax) {
		t->tm_mday++;
		t->tm_yday++;
		t->tm_hour = 0;
		return;
	} else if (mm < 12) {
		t->tm_mon++;
		t->tm_mday = 1;
		t->tm_yday++;
		t->tm_hour = 0;
		cout << " " << t->tm_mon + 1;
		flush(cout);
		return;
	} else {
		t->tm_year++;
		t->tm_mon = 0;
		t->tm_mday = 1;
		t->tm_yday = 0;
		t->tm_hour = 0;
//		cout << "\nYear = " << (t->tm_year < 100 ? 19 : 20) << t->tm_year << " Month = 1 ";
		cout << "\nYear = " << 1900+ t->tm_year  << " Month = 1 ";
		flush(cout);
		return;
	}

	return;
}


int mainFunction(int syyyy,	int	smm,	int sdd,	int shh,	int eyyyy,	int emm,	int edd,	int ehh,	int NumGrid,	int *Col,	int *Row)
{
	/*
	int syyyy=1980;
	int	smm=1;
	int sdd=1;
	int shh=5;
	int eyyyy=2012;
	int emm=1;
	int edd=1;
	int ehh=5;
	int NumGrid=1;
	int row=90;
	int col=367;
	*/
	char NLDAS_FOLDER[1000]="/home/ualjxp/NLDAS/back";
	char OUT_FOLDER[1000]="/home/ualjxp/NLDAS/ext_yan";

	GRIDS *grid;
	grid = (GRIDS *) malloc(NumGrid * sizeof(GRIDS));
	cout << "\nmalloc = " << NumGrid * sizeof(GRIDS);
	for(int GridNo=0;GridNo<NumGrid;GridNo++){
		grid[GridNo].row=Row[GridNo];
		grid[GridNo].col=Col[GridNo];
	}
	//June	NLDAS_ID	xx356y87 x357y87 x354y86 x355y86 x356y86 x354y85 x355y85 x356y85	(from 01/01/1979 to 12/31/2017)		1979 01  01  00  2017 12  31  23  8  356  87  357  87  354  86  355  86  356  86  354  85  355  85  356  85



	/*
	char            NLDAS_FOLDER[1000];
	char            OUT_FOLDER[1000];
	cin.getline(NLDAS_FOLDER, 1000, '\n');
	cin.getline(OUT_FOLDER, 1000, '\n');
	//cin >> NLDAS_FOLDER;
	//cin >> OUT_FOLDER;

	cout << "NLDAS DIR  = " << NLDAS_FOLDER << "NLDAS_FORA0125_H.002\n";
	cout << "OUTPUT DIR = " << OUT_FOLDER << "\n";

	int             syyyy, smm, sdd, shh;
	int             eyyyy, emm, edd, ehh;

	GRIDS          *grid;
	//ifstream rowsxcols;
	//rowsxcols.open("IN/rowsxcols.txt");

	//rowsxcols >> syyyy;
	//rowsxcols >> smm;
	//rowsxcols >> sdd;
	//rowsxcols >> shh;
	//rowsxcols >> eyyyy;
	//rowsxcols >> emm;
	//rowsxcols >> edd;
	//rowsxcols >> ehh;

	cin >> syyyy;
	cin >> smm;
	cin >> sdd;
	cin >> shh;

	cin >> eyyyy;
	cin >> emm;
	cin >> edd;
	cin >> ehh;

	int             NumGrid;
	//rowsxcols >> NumGrid;
	cin >> NumGrid;
	DEBUG           cout << NumGrid;

	grid = (GRIDS *) malloc(NumGrid * sizeof(GRIDS));

	for (int i = 0; i < NumGrid; i++) {
		//rowsxcols >> grid[i].col;
		//rowsxcols >> grid[i].row;
		cin >> grid[i].col;
		cin >> grid[i].row;
		//cout << "\nRunning... " << grid[i].col << " x " << grid[i].row << "\n";
	}
	*/
	char            row[10], col[10], filename[2500];
	ofstream      **files = new ofstream *[NumGrid];
	for (int i = 0; i < NumGrid; i++)
		files[i] = new ofstream[7];

	for (int i = 0; i < NumGrid; i++) {
		sprintf(row, "%d", grid[i].row);
		sprintf(col, "%d", grid[i].col);
		sprintf(filename, "%s/x%dy%dzPP.txt", OUT_FOLDER, grid[i].col, grid[i].row);
		DEBUG           cout << "\n" << filename;
		files[i][0].open(filename, ios::out);
		files[i][0] << fixed;

		sprintf(filename, "%s/x%dy%dzTT.txt", OUT_FOLDER, grid[i].col, grid[i].row);
		files[i][1].open(filename, ios::out);
		files[i][1] << fixed;
		sprintf(filename, "%s/x%dy%dzRH.txt", OUT_FOLDER, grid[i].col, grid[i].row);
		files[i][2].open(filename, ios::out);
		files[i][2] << fixed;
		sprintf(filename, "%s/x%dy%dzWD.txt", OUT_FOLDER, grid[i].col, grid[i].row);
		files[i][3].open(filename, ios::out);
		files[i][3] << fixed;
		sprintf(filename, "%s/x%dy%dzRN.txt", OUT_FOLDER, grid[i].col, grid[i].row);
		files[i][4].open(filename, ios::out);
		files[i][4] << fixed;
		sprintf(filename, "%s/x%dy%dzVP.txt", OUT_FOLDER, grid[i].col, grid[i].row);
		files[i][5].open(filename, ios::out);
		files[i][5] << fixed;
	}



	double          gridValue[7], V, VP, VPsat;
	GDALDataset    *temp;
	GDALAllRegister();

	struct tm       stime, etime;
	stime.tm_year = syyyy - 1900;
	stime.tm_mon = smm - 1;
	stime.tm_mday = sdd;
	stime.tm_sec = 0;
	stime.tm_min = 0;
	stime.tm_hour = shh;
	stime.tm_isdst = -1;
	etime.tm_year = eyyyy - 1900;
	etime.tm_mon = emm - 1;
	etime.tm_mday = edd;
	etime.tm_sec = 0;
	etime.tm_min = 0;
	etime.tm_hour = ehh + 1;
	etime.tm_isdst = -1;

	//DEBUG           cout << "\nStart time is: " << asctime(&stime);
	//DEBUG           cout << "\nEnd   time is: " << asctime(&etime);
	//cout << "\nEnd   time is: " << asctime(&etime);
	//cout << "\nStart time is: " << asctime(&stime);

	time_t          s_time, e_time;
	s_time = mktime(&stime);
	e_time = mktime(&etime);
	int             timediffsec = difftime(e_time, s_time);
	int             timediffhrs = timediffsec / 3600;
	DEBUG           cout << "\nTime diff= " << timediffsec;
	//temp = (GDALDataset *) GDALOpen("NLDAS_FORA0125_H.A19790101.1300.002.grb", GA_ReadOnly);

	char            grbFileName[200], fileA[100], fileB[20], stry[5],
	                strm[3], str3d[4], str2d[3], strh[5];
	int             yyyy, mm, dd, dd2, hh;

	cout << "Running...\n";
	double          est = ((23.0 / (2 * 60)) / 8760) * NumGrid * timediffhrs;
	cout << "Estimated Run Time: " << est/2 << " Hrs. \n";
	if (est > 1.0)
		cout << "\nTime for some coffee?\n";
	else if (est > 24.0)
		cout << "\nTry breaking the jobs into pieces!!\n";
	else if (est > 24.0 * 7)
		cout << "\nMore than a week? seriously?\nConsider breaking the input into smaller groups\n";
//	cout << "\nYear = " << (stime.tm_year < 100 ? 19 : 20) << stime.tm_year << " Month = " << stime.tm_mon + 1 << " ";

	cout << "\nYear = " << 1900+ stime.tm_year  << " Month = " << stime.tm_mon + 1 << " ";
	flush(cout);
	ifstream        grbFile;
	for (int i = 0; i < timediffhrs; i++) {

		yyyy = stime.tm_year;
		mm = stime.tm_mon + 1;
		dd = stime.tm_mday;
		dd2 = stime.tm_yday + 1;
		hh = stime.tm_hour;

		if (yyyy < 100)
			sprintf(stry, "19%2d", yyyy);
		else
			sprintf(stry, "20%02d", yyyy - 100);
		sprintf(strm, "%02d", mm);
		sprintf(str2d, "%02d", dd);
		sprintf(str3d, "%03d", dd2);
		sprintf(strh, "%04d", hh * 100);
		sprintf(fileA, "%s/%s/%s", stry, str3d, "NLDAS_FORA0125_H.A");
		sprintf(fileB, "%s", ".002.grb");
		//sprintf(grbFileName, "%s/NLDAS_FORA0125_H.002/%s%s%s%s.%s%s", NLDAS_FOLDER, fileA, stry, strm, str2d, strh, fileB);
		sprintf(grbFileName, "%s/%s%s%s%s.%s%s", NLDAS_FOLDER, "NLDAS_FORA0125_H.A", stry, strm, str2d, strh, fileB);

		DEBUG           cout << "\n" << grbFileName;


		temp = (GDALDataset *) GDALOpen(grbFileName, GA_ReadOnly);
		for (int j = 0; j < NumGrid; j++) {
			DEBUG           cout << "$#" << getRasterValue(temp, 10, grid[j].row, grid[j].col) << "\n";
			gridValue[0] = getRasterValue(temp, 10, 224 - grid[j].row, grid[j].col - 1);
			//Precip
				gridValue[1] = getRasterValue(temp, 1, 224 - grid[j].row, grid[j].col - 1);
			//Temp
				gridValue[2] = getRasterValue(temp, 2, 224 - grid[j].row, grid[j].col - 1);
			//Sp.RH
				gridValue[3] = getRasterValue(temp, 4, 224 - grid[j].row, grid[j].col - 1);
			//Ux
				gridValue[4] = getRasterValue(temp, 11, 224 - grid[j].row, grid[j].col - 1);
			//ShortWave
				gridValue[5] = getRasterValue(temp, 5, 224 - grid[j].row, grid[j].col - 1);
			//Uy
				gridValue[6] = getRasterValue(temp, 3, 224 - grid[j].row, grid[j].col - 1);
			//P atm

				// PRECIPITATION TO THE FILE
				files[j][0] << setprecision(5) << ((double) i) / 24 << "\t" << setprecision(6) << gridValue[0] * 0.024 << "\t";
			files[j][0] << setprecision(5) << ((double) i) / 24 + (1.0 / 24) - 0.00001 << "\t" << setprecision(6) << gridValue[0] * 0.024 << "\n";
			//TEMPERATURE TO THE FILE
				files[j][1] << setprecision(5) << ((double) i) / 24 << "\t" << setprecision(2) << gridValue[1] << "\t";
			files[j][1] << setprecision(5) << ((double) i) / 24 + (1.0 / 24) - 0.00001 << "\t" << setprecision(2) << gridValue[1] << "\n";


			VP = gridValue[6] * (gridValue[2] / (gridValue[2] + 0.62198));
			VPsat = 611.2 * exp(17.62 * gridValue[1] / (243.12 + gridValue[1]));
			V = sqrt(gridValue[3] * gridValue[3] + gridValue[5] * gridValue[5]) * 86400.0;

			//RELATIVE HUMIDITY TO THE FILE
				files[j][2] << setprecision(5) << ((double) i) / 24 << "\t" << setprecision(2) << VP / VPsat << "\t";
			files[j][2] << setprecision(5) << ((double) i) / 24 + (1.0 / 24) - 0.00001 << "\t" << setprecision(2) << VP / VPsat << "\n";

			//WIND VELOCITY
				files[j][3] << setprecision(5) << ((double) i) / 24 << "\t" << setprecision(2) << V << "\t";
			files[j][3] << setprecision(5) << ((double) i) / 24 + (1.0 / 24) - 0.00001 << "\t" << setprecision(2) << V << "\n";

			//SOLAR RADIATION
				files[j][4] << setprecision(5) << ((double) i) / 24 << "\t" << setprecision(2) << 86400.0 * gridValue[4] << "\t";
			files[j][4] << setprecision(5) << ((double) i) / 24 + (1.0 / 24) - 0.00001 << "\t" << setprecision(2) << 86400.00 * gridValue[4] << "\n";

			//VAPOR PRESSURE
				files[j][5] << setprecision(5) << ((double) i) / 24 << "\t" << setprecision(2) << VP << "\t";
			files[j][5] << setprecision(5) << ((double) i) / 24 + (1.0 / 24) - 0.00001 << "\t" << setprecision(2) << VP << "\n";
		}
		GDALClose(temp);
		stepOneHour(&stime);
	}

	DEBUG           cout << "\n";

	for (int i = 0; i < NumGrid; i++) {
		for (int j = 0; j < 6; j++) {
			files[i][j].flush();
			files[i][j].close();
		}
	}
	if (est <= 0) {
		cout << "\n\nSUM TING RONG\n\n";
		return 0;
	}
	cout << "\n\nJob finished successfully!\n";
}
