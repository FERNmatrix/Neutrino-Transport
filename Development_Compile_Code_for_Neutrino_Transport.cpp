// g++ -o NT Development_Compile_Code_for_Neutrino_Transport.cpp
// Development Compile Code for Neutrino Transport Network
// Raghav Chari, Adam Cole
// #include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <ctime>

#include <fstream>
using namespace std;
 using std::string;

// Function signatures in main:

// Global Variables

double Model = '001';
double N_g = 40;
double stoptime = 1.0e2;
double t = 1.0e-19;
double t_W0 = 1.0e-10;
double dt_min = 1.0e-16;
double t_max = 1.0e-01;
double dt = 1.0;
double dt_grw = 1.003;
double dt_dec = 0.90;

double dt_FE  = dt;
double dt_EA  = dt;
double dt_PE  = dt;

bool reStart  = 0;

double G_A = 7.5d-01;    
double G_B = 1.0d+02;    
double G_C = sqrt(50.0);

int cycleM = 10e9;    
int cycleD = 100;     
int cycleW = 10;      

double tolPE = 1.0e-01;  
double tolC  = 1.0e-09;

int FE    = 0;
int EA    = 1;
int FE_PE = 2;
int QSS1  = 4;
int QSS2  = 5;
int BE    = 6;

char scheme[18] = "ExplicitAsymptotic";

char Comment = '';
char PlotFileDir[6] = "Output";
char PlotFileName[8] = "PlotFile";
int PlotFileNumber = 0;
int nPlotFiles = 100; 



 Class defenitions 

 class Utilities{
	private:
	public:
		double eC[N_G];
		double de[N_G];

 		void initalizeNES(model, N_G){
 			double *eC;
 			eC = ReadData1D('',N_G);
 			double *de;
 			de = ReadData1D('',N_G)
 			double dV[N_G];
 			for (int i = 0; i <= (N_G-1); ++i){

 			dV[i]= ( (eC[i]+0.5d0*de[i]).^3 - (eC[i]-0.5d0*de[i]).^3 )/3.0d0;

 			}
 			double *R_In;
 			R_In = ReadData2D('',N_G,N_G);

 			R_Out = conj(R_In);

 			switch(model){
 				case 1:
 					int mu = 145.254;
 					int kT = 20.5399;
 					N_Eq[i] = 1.0 ./  exp( (eC[i]-mu[i])./kT ) + 1.0;
 					for (int j = 1; j <= N_G; ++j){
 						for (int i = 1; i <= N_G; ++i){
 							if( j < i )
 								R_In(i,j) = R_In(j,i) * exp( ( eC(j) - eC(i) ) / kT );
 						}
 					}
 				case 2:	
 					int mu = 045.835;
 					int kT = 15.9751;
 					N_Eq[i] = 1.0 ./  exp( (eC[i]-mu[i])./kT ) + 1.0;
 					for (int j = 1; j <= N_G; ++j){
 						for (int i = 1; i <= N_G; ++i){
 							if( j < i )
 								R_In(i,j) = R_In(j,i) * exp( ( eC(j) - eC(i) ) / kT );
 						}
 					}
 				case 3:
					int mu = 020.183;
 					int kT = 07.7141;
 					N_Eq[i] = 1.0 ./  exp( (eC[i]-mu[i])./kT ) + 1.0;
 					for (int j = 1; j <= N_G; ++j){
 						for (int i = 1; i <= N_G; ++i){
 							if( j < i )
 								R_In(i,j) = R_In(j,i) * exp( ( eC(j) - eC(i) ) / kT );
 						}
 					}
 				case 4:
 					int mu = 009.118;
 					int kT = 07.5830;
 					N_Eq[i] = 1.0 ./  exp( (eC[i]-mu[i])./kT ) + 1.0;
 					for (int j = 1; j <= N_G; ++j){
 						for (int i = 1; i <= N_G; ++i){
 							if( j < i )
 								R_In(i,j) = R_In(j,i) * exp( ( eC(j) - eC(i) ) / kT );
 						}
 					}
 				case 5:
 					int mu = 003.886;
 					int kT = 03.1448;
 					N_Eq[i] = 1.0 ./  exp( (eC[i]-mu[i])./kT ) + 1.0;
 					for (int j = 1; j <= N_G; ++j){
 						for (int i = 1; i <= N_G; ++i){
 							if( j < i )
 								R_In(i,j) = R_In(j,i) * exp( ( eC(j) - eC(i) ) / kT );

// 							R_Out = R_In'

 							
 						}
 					}
 				default: N_Eq = 1.0
 	}



// Main

int main() {

	FILE *pFile;
	pFile = fopen("output/NT.data","w");
	int totalTimeSteps = 0;
	Logfile(stop_time, t, t_W0, dt, G_A, G_B, G_C, tolC, Scheme, Comment);

	initializeNES( Model, N_g);

	eC = get_eC();
	dV = get_dV();
	R_In = get_R_In();
	R_Out = get_R_Out();
	N_Eq = get_N_Eq();

	cout << " t       dt       totalTimeSteps       ";

	fprintf(pFile, " t       dt       totalTimeSteps       "); 

	while(t < stoptime){

cout<< t<< "\n";
//		cout << "/n%7.4d %7.4f", t, dt;
fprintf(pFile, "/n%7.4f %7.4f", t, dt); 
// 
	t += dt;
	totalTimeSteps ++;
	}

}




// Functions Used in Main

double get_eC(){

double AAA = Utilities::eC;

return AAA;

}

double get_dV(){

double AAA = Utilities::dV;

return AAA;

}

double get_R_In(){

double AAA = Utilities::R_In;

return AAA;

}

double get_R_Out(){

double AAA = Utilities::R_Out;

return AAA;

}

double get_N_Eq(){

double AAA = Utilities::N_Eq;

return AAA;

}

double * ReadData1D(FileName, N) {
 FILE *fr;
 fr = fopen(FileName, "r")
	double Data1D[40];

	while (!feof(fr)){

		fread(Data1D, sizeof(Data1D), 1, fr);
		cout << Data1D
	}
	fclose (fr);
	return Data1D;
}

 void ReadData2D(FileName, N) {
	File *fr;
	fr = fopen(FileName, "r")
	double Data2D[40];
	while (!feof(fr)){

		fread(Data2D, sizeof(Data2D), 1, fr);
		cout << Data2D

		int flat[100];
    for (int i=0; i<100; ++i)
       flat[i]=i;
   		   int (&square)[10][10] = *reinterpret_cast<int(*)[10][10]>(flat);
   for (int i=0; i<10; ++i)
        cout << square[5][i] << ' ';
   cout << '\n';
   fclose (fr);
   return Data2D

}

	void logfile(t_end, t, t_W0, dt, G_A, G_B, G_C, tolC, Scheme, Comment) 
	FILE *Log
	Log = fopen("Log.txt", "a");
	if(log!=NULL){
		fprintf(Log,' t_end =%5d\n', t_end );
    	fprintf(Log,' t     =%5d\n', t );
    	fprintf(Log,' t_W0  =%5d\n', t_W0 );
   		fprintf(Log,' dt    =%5d\n', dt );
    	fprintf(Log,' G_A   =%5d\n', G_A );
    	fprintf(Log,' G_B   =%5d\n', G_B );
    	fprintf(Log,' G_C   =%5d\n', G_C );
    	fprintf(Log,' tolC  =%5d\n', tolC );
    	fprintf(Log,' Scheme=%s\n', Scheme );
    	fprintf(Log,' Comment: %s\n', Comment );
    	fprintf(Log, '\n' );
    	fclose(Log);  
	}






