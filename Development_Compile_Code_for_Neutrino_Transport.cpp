/* 
*g++ Development_Compile_Code_for_Neutrino_Transport.cpp -o -lgsl -lgslcblas -lm -lstdc++
*g++ -Wall -I/opt/homebrew/Cellar/gsl/2.7.1/include -c Development_Compile_Code_for_Neutrino_Transport.cpp

*Development Compile Code for Neutrino Transport Network

AUTHORS:
* ---------------
* Raghav Chari
* Adam Cole
*----------------

*/
#include "stdio.h"
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <ctime>
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <string.h>
static FILE* pFile; 


using namespace std;
 using std::string;

#define N_G 40
clock_t startCPU, stopCPU;
#define START_CPU if ((startCPU=clock())==-1) {printf("Error calling clock"); exit(1);}
#define STOP_CPU if ((stopCPU=clock())==-1) {printf("Error calling clock"); exit(1);}
#define PRINT_CPU (printf("Timer: %7.4e ms used", 1000*(double)(stopCPU-startCPU)/CLOCKS_PER_SEC));
#define FPRINTF_CPU (printf(pFile, "in %7.4e seconds\n", (double)(stopCPU-startCPU)/CLOCKS_PER_SEC));
#define FPRINTF_CPU2 (printf(pFile2, "in %7.4e seconds\n", (double)(stopCPU-startCPU)/CLOCKS_PER_SEC));
#define FPRINTF_CPUD (printf(pFileD, "in %g seconds\n", (double)(stopCPU-startCPU)/CLOCKS_PER_SEC));
#define PRINT_CPU_TEST (printf("\nTimer Test: %g ms used by CPU\n", 1000*(double)(stopCPU-startCPU)/CLOCKS_PER_SEC));

// Global Variables

int model = 1;
double stoptime = 1.0e2;
double t = 1.0e-15;
double t_W0 = 1.0e-10;
double dt_min = 1.0e-16;
double t_max = 1.0e-01;
double dt = 1.0e-15;
double dt_grw = 1.003;
double dt_dec = 0.90;
double t_end = 0.0100;



double dt_FE  = dt;
double dt_EA  = dt;
double dt_PE  = dt;

bool reStart  = 0;

double G_A = 7.5e-01;    
double G_B = 1.0e+02;    
double G_C = sqrt(50.0);

double cycleM = 10e9;    
int cycleD = 100;     
int cycleW = 10;      

double tolPE = 1.0e-01;  
double tolC  = 1.0e-09;
double tolN = 0.0100;
double Tol;


int FE    = 0;
int EA    = 1;
int FE_PE = 2;
int QSS1  = 4;
int QSS2  = 5;
int BE    = 6;

int Scheme = 1;

char Comment = ' ';
char PlotFileDir[] = "Output";
char PlotFileName[] = "PlotFile";
int PlotFileNumber = 0;
int nPlotFiles = 100; 

// Global Arrays used in Main

double ** tempA[N_G][N_G];
double * tempKvec[N_G];
double ** tempF[N_G][N_G];
double ** tempKmatrix[N_G][N_G];
double ** tempFp[N_G][N_G];
double ** tempkp[N_G][N_G];
double *N_0[N_G];
double *Nold[N_G];
double *Nnew[N_G];
double *eC[N_G];
double *de[N_G];
double *dV[N_G];
double **R_In[40][40];
double **R_Out[40][40];
double *N_Eq[N_G];
double *N[N_G];

// Function Definitions
double max_element(double __lcpp_x);
double min_element(double __lcpp_x);
double sum(double, double);

// Write Plotfile

int  wrtCount;
int  wrtTimes;
int logspace();
int Branch;
int FileNumber = 0;
 
//myFunctions

void logfile() {
	FILE* Log = fopen("logFile.txt", "a");
  if (Log!=NULL){
	std::fprintf(Log," t_end =%5d\n", t_end );
    std::fprintf(Log," t     =%5d\n", t );
      std::fprintf(Log," t_W0  =%5d\n", t_W0 );
      std::fprintf(Log," dt    =%5d\n", dt );
      std::fprintf(Log," G_A   =%5d\n", G_A );
      std::fprintf(Log," G_B   =%5d\n", G_B );
      std::fprintf(Log," G_C   =%5d\n", G_C );
      std::fprintf(Log," tolC  =%5d\n", tolC );
      std::fprintf(Log," Scheme=%s\n", Scheme );
      std::fprintf(Log," Comment: %s\n", Comment );
      std::fprintf(Log, "\n" );
      fclose(Log);   
  }
}


// Function signatures in main:

// Functions Used in Main

// Get Functions

double ** get_cMat(){
  return ** tempA;
}
double * get_kVec(){
  return * tempKvec;
}

double ** get_F0(){
  return ** tempF;
}

double ** get_k0(){
  return ** tempKmatrix;
}
double ** get_kp(){
  return ** tempkp;
}
double ** get_Fp(){
  return ** tempFp;
}
// Write_Plotfile

int Write_Plotfile() {

std::fprintf(pFile, "%7.4e %7.4e", t, dt);
FileNumber = FileNumber + 1; 

return  FileNumber; 

}



//Utilities initalizeNES 
double * get_eC(){

double * AAA = *eC;

return AAA;
}

//Utilities initalizeNES;
double * get_dV(){

double * AAA = *dV;

return AAA;

}
//Utilities initalizeNES;
double ** get_R_In(){

double ** AAA = **R_In;

return AAA;

}
//Utilities initalizeNES;
double ** get_R_Out(){

double ** AAA =  **R_Out;

return AAA;

}
//Utilities initalizeNES;
double * get_N_Eq(){

double * AAA = *N_Eq;

return AAA;

}

// BuildCollisionMatrix 

void BuildCollisionMatrix(double **A, double **K, double N_Eq_i, double N_Eq_j, double diff_i, double diff_j) {

for (int i = 0; i < N_G; ++i){
    for (int j = 0; j < N_G; ++j){
       A[i][j] = 0;
  }
}

for (int i = 0; i < N_G; ++i){
  for (int j = 0; j < N_G; ++j){
   K[i][j] = 0;
  }
}

for (int i = 0; i < N_G; ++i){
    for (int j = 0; j < N_G; ++j){
      if( i != j){
        double B = ( **R_In[i][j] * *dV[i] + **R_Out[i][j] * *dV[j] );
            double C = *dV[i] * *N[i] + *dV[j] * *N[j];
        
           double a = ( R_In[i][j] - R_Out[i][j] ) * *dV[i];
            double b = B + ( R_In[i][j] - R_Out[i][j] ) * C;
            double c = **R_In[i][j] * C;
            double d = (b*b) - (4.0 * a * c);
        
              N_Eq_i = 0.5 * ( b - sqrt( d ) ) / a;
              N_Eq_j = ( C - N_Eq_i * *dV[i] ) / *dV[j];
      }
      else {
          N_Eq_i = *N[i];
          N_Eq_j = *N[j];
        }
      
    }
}

for (int i = 0; i < N_G; ++i){
	for (int j = 0; j < N_G; ++j){
		if( i != j){
	diff_i = abs( N_Eq_i - *N[i] ) / max( 1.0e-16, N_Eq_i );
    diff_j = abs( N_Eq_j - *N[j] ) / max( 1.0e-16, N_Eq_j );
    if( diff_i > Tol || diff_j > Tol )  {
	


      A[i][j] = A[i][j]+ (1.0 - *N[i]) * **R_In[i][j] * *dV[j];
      A[i][j] = A[i][j] - **R_Out[i][j] * *dV[j] * (1.0 -  *N[j]);
      *K[i] = *K[i] + **R_Out[i][j] * *dV[j] + ( (**R_In[i][j] * *dV[j] - **R_Out[i][j] * *dV[j] * *N[j]));
	}
		}
	}
	
}
}
//ComputeRates


void ComputeRates( double F[40][40], double k[40][40]) {
  **F = **R_In[40][40] * ( *N[40] * *dV[40] );
    **k = **F + **R_Out[40][40] * ( ( 1.0 - *N[40] ) * *dV[40] );

  }


 // class defenitions 

 class Utilities{
  private:
  public:
    void initalizeNES(int model, const int [N_G]);
    double static * ReadData1D(char FileName[], double N);
    double static ** ReadData2D(char FileName[], double M, double N);
    static void startTimer(){
      START_CPU // Start a timer for rate calculation
    }

    static void stopTimer(){

      STOP_CPU;
      printf("\n");
      PRINT_CPU;
      printf("\n");
    }
  };



      void Utilities::initalizeNES( int model, const int [N_G]){
      *eC = ReadData1D("Data/NES_RATES_EnergyBinCenter.dat",N_G);
      *de = ReadData1D("Data/NES_RATES_EnergyBinWidth.dat",N_G);
      

      
    for (int i = 0; i <= (N_G-1); ++i){
    
    *dV[i]= (pow(*eC[i]+0.5e0**de[i],3) - pow(*eC[i]-0.5e0**de[i],3))/3.0e0;


      }

    **R_In = ReadData2D("Data/NES_RATES_R_In___001.dat",N_G,N_G);


    //int gsl_matrix_complex_conjtrans_memcpy(gsl_matrix *R_Out, const gsl_matrix *R_In); // This function makes R_Out the conjugate transpose of R_In 

      switch(model) {
        case 1: {
          double mu = 145.254;
          double kt = 20.5399;
          for (int i = 1; i <= N_G; ++i){
             *N_Eq[N_G] = 1 /  (exp( (*eC[i]-mu)/kt ) + 1);
          }
          for (int j = 1; j <= N_G; ++j){
            for (int i = 1; i <= N_G; ++i){
              if( j < i )
                **R_In[j][i] = **R_In[j][i] * exp( ( eC[j] - eC[i] ) / kt );
            }
          }

        break;
        }
        //int gsl_matrix_complex_conjtrans_memcpy(gsl_matrix *R_Out, const gsl_matrix *R_In);

        case 2: {
           double mu = 045.835;
           double kt = 15.9751;
          for (int i = 1; i <= N_G; ++i){
             *N_Eq[N_G] = 1 /  (exp( (*eC[i]-mu)/kt ) + 1);
          }
          for (int j = 1; j <= N_G; ++j){
            for (int i = 1; i <= N_G; ++i){
              if( j < i )
                **R_In[j][i] = **R_In[j][i] * exp( ( eC[j] - eC[i] ) / kt );
            }
          }

        break;
        }
        //int gsl_matrix_complex_conjtrans_memcpy(gsl_matrix *R_Out, const gsl_matrix *R_In);

        case 3: {
           double mu = 020.183;
           double kt = 07.7141;
          for (int i = 1; i <= N_G; ++i){
             *N_Eq[N_G] = 1 /  (exp( (*eC[i]-mu)/kt ) + 1);
          }
          for (int j = 1; j <= N_G; ++j){
            for (int i = 1; i <= N_G; ++i){
              if( j < i )
                **R_In[j][i] = **R_In[j][i] * exp( ( eC[j] - eC[i] ) / kt );
            }
          }

        break;
        }
        //int gsl_matrix_complex_conjtrans_memcpy(gsl_matrix *R_Out, const gsl_matrix *R_In);

        case 4: {
           double mu = 009.118;
           double kt = 07.5830;
          for (int i = 1; i <= N_G; ++i){
             *N_Eq[N_G] = 1 /  (exp( (*eC[i]-mu)/kt ) + 1);
          }
          for (int j = 1; j <= N_G; ++j){
            for (int i = 1; i <= N_G; ++i){
              if( j < i )
                **R_In[j][i] = **R_In[j][i]  * exp( ( eC[j] - eC[i] ) / kt );
            }
          }

        break;
        }
        //int gsl_matrix_complex_conjtrans_memcpy(gsl_matrix *R_Out, const gsl_matrix *R_In);
        case 5: {
           double mu = 003.886;
           double kt = 03.1448;
          for (int i = 1; i <= N_G; ++i){
             *N_Eq[N_G] = 1 /  (exp( (*eC[i]-mu)/kt ) + 1);
          }
          for (int j = 1; j <= N_G; ++j){
            for (int i = 1; i <= N_G; ++i){
              if( j < i )
                **R_In[j][i] = **R_In[j][i]  * exp( ( eC[j] - eC[i] ) / kt );
              
            }
          }

        break;
        } 
        // int gsl_matrix_complex_conjtrans_memcpy(gsl_matrix *R_Out, const gsl_matrix *R_In);
        default: *N_Eq[N_G] = 1.0;
      
      }
}
  

double * Utilities::ReadData1D(char FileName[], double N) {
 FILE *fr;
 fr = fopen(FileName, "r");
  double Data1D[40];

  while (!feof(fr)){

    fread(Data1D, sizeof(Data1D), 1, fr);
    std::cout << Data1D;
  }
  fclose (fr);
  return Data1D;
}

double ** Utilities::ReadData2D(char FileName[], double M, double N) {
  FILE *fr;
  fr = fopen(FileName, "r");
  double ** Data2D[40][40];
  while (!feof(fr)){

    fread(Data2D, sizeof(Data2D), 1, fr);
    std::cout << ReadData2D;

    int flat[100];
    for (int i=0; i<100; ++i)
       flat[i]=i;
         int (&square)[10][10] = *reinterpret_cast<int(*)[10][10]>(flat);
   for (int i=0; i<10; ++i)
        std::cout << square[5][i] << ' ';
   std::cout << '\n';
   fclose (fr);
   return ** Data2D;
  }
}

// Main

int main() {
Utilities::startTimer();

  FILE *pFile;
  pFile = fopen("output/NT.data","w");
  int totalTimeSteps = 0;
  logfile();

std::fprintf(pFile, "# data_output \n");
std::fprintf(pFile, "#       t        dt  %\n");

Utilities initalizeNES; 

  double multi[40][40];
  double **a;

  for (int i = 0; i < N_G; ++i){
    for (int j = 0; j < N_G; ++j){
      if (i = j) a[i][j] = *dV[i];
      else if (i != j) a[i][j] = 0;
    }
  }

  for (int i = 0; i < N_G; ++i){
    for (int j = 0; j < N_G; ++j){
      multi[i][j] = 0;
    }
  }

  for (int i = 0; i < N_G; ++i){
    for (int j = 0; j < N_G; ++j){
      for(int k = 0; k < N_G; ++k){
        multi[i][j] += **R_In[i][k] * a[k][j]; 

      }
    }
  }

double R_In_H = multi[40][40];

  for (int i = 0; i < N_G; ++i){
    for (int j = 0; j < N_G; ++j){
      multi[i][j] = 0;
    }
  }

  for (int i = 0; i < N_G; ++i){
    for (int j = 0; j < N_G; ++j){
      for(int k = 0; k < N_G; ++k){
        multi[i][j] += **R_Out[i][k] * a[k][j]; 

      }
    }
  }

double R_Out_H = multi[40][40];

if (reStart == true) {

}

else {

  for (int i = 0; i < N_G; ++i) {
*N_0[N_G] = G_A * exp( - 0.5 * pow(( ( *eC[i] - G_B ) / G_C ),2));
  
  }
}


for (int i = 0; i < N_G; ++i){
  Nold[N_G] = N_0[N_G];
}


// BuildCollisionMatrix(R_In, R_Out, Nold, dV, const int [N_G], 0.0);
BuildCollisionMatrix;

double ** cMat[40][40];
cMat[40][40] = get_cMat();
double * kVec[40];
kVec[40] = get_kVec();

bool done   = false;
bool reStep = false;
int cycle  = 0;
int true_cycle = 0;
int nIterations = 0;
int nTrueIterations = 0;
int maxFPIterations = 10000; 
int mAA = 3; 


//PlotFileNumber = Write_Plotfile ( t, dt, Nold, eC, dV, kVec, cMat, 0, FE, PlotFileDir, PlotFileName, PlotFileNumber );
PlotFileNumber = Write_Plotfile();
//wrtTimes  =  logspace( log10(t_W0), log10(t_end), nPlotFiles );
wrtTimes  =  logspace();
wrtCount = 1;

while (done != true) {

 true_cycle = true_cycle + 1; 

  if (reStep != true) {

  cycle = cycle +1;

  }

 // Nold = ApplyPerturbation(Nold, amp dV, N_G); 

  switch (Scheme) {

  case 1: {

  Branch = EA;

  //void BuildCollisionMatrix(R_In, R_Out, Nold, dV, N_G, 0.0);
  BuildCollisionMatrix;

    cMat[40][40] = get_cMat();
    kVec[40] = get_kVec();
  	dt_FE = min_element(1.0 / *kVec[40]);

  if (reStep) {

    dt = (dt_dec * dt);
    reStep = false;
  }
    else {
      dt = dt_grw * dt;
	}	
  

  if (dt <= dt_FE) {

    *Nnew[N_G] = *Nold[N_G] + dt * (**cMat[40][40] * *Nold[N_G]);
}
  else {

    *Nnew[N_G] = *Nold[N_G] + (dt / (1.0 + *kVec[40] *dt)) * (**cMat[40][40] * *Nold[N_G]);
}

if( abs( sum(*Nnew[N_G], *dV[40]) - sum( *Nold[N_G], *dV[40] ) )  / sum ( *Nold[N_G] ,*dV[40] ) > tolC ) {

  reStep = true;
}


if ( max_element( abs( Nnew - Nold ) / max( *Nold[N_G], 1.0e-8 ) ) > tolN ) {

  reStep = true;

}

  }
  case 2: {

  Branch = FE_PE;

  //void BuildCollisionMatrix(R_In, R_Out, Nold, dV, N_G, tolPE);
  BuildCollisionMatrix;

  cMat[40][40] = get_cMat();
  kVec[40] = get_kVec();


  dt = min( min_element( 1.0 / *kVec[40] ), dt_grw * dt);

  *Nnew[N_G] = *Nold[N_G] + dt * ( **cMat[40][40] * *Nold[N_G]);
  }

  case 4: {

  Branch = QSS1;

  if (reStep) {
    dt = dt_dec * dt;
    reStep = false;
  }

    else
      dt = dt_grw * dt;
  

  ComputeRates;
 
  double ** F0[40][40];
  F0[40][40] = get_F0();
  double ** k0[40][40];
  k0[40][40] = get_k0();
  
  double exp_kdt = exp( - dt * **k0[40][40] );
      
    dt_FE = min_element( 1.0 / **k0[40][40] );
      
    if( dt <= dt_FE ) {

    //void BuildCollisionMatrix(R_In, R_Out, Nold, dV, N_G, 0.0); 
	BuildCollisionMatrix;

    cMat[40][40] = get_cMat();
    kVec[40] = get_kVec();


    *Nnew[N_G] = *Nold[N_G] + dt * ( **cMat[40][40] * *Nold[N_G] );

    nTrueIterations = nTrueIterations + 1 ;

  }

    else {
      //void BuildCollisionMatrix(R_In, R_Out, Nold, dV, N_G, 0.0);
	  BuildCollisionMatrix;

      cMat[40][40] = get_cMat();
      kVec[40] = get_kVec();


      *Nnew[N_G] = *Nold[N_G] + ( **cMat[40][40] * *Nold[N_G] ) * ( 1.0 - exp_kdt ) / **k0[40][40];
    }

    if ( abs( sum(*Nnew[N_G], *dV[40]) - sum(*Nold[N_G], *dV[40]) ) / sum(*Nold[N_G], *dV[40]) > tolC ) {
      reStep = true;

      }
  }
	case 5: {


 Branch = QSS2;

    if (reStep) {
      	dt = dt_dec * dt; 
        reStep = false;
	}
      else {
        dt = dt_grw * dt;
	  }
  

    //void ComputeRates( R_In, R_Out, Nold, dV);
	ComputeRates;
    double ** F0[40][40];
  	F0[40][40] = get_F0();
  	double ** k0[40][40];
  	k0[40][40] = get_k0();
      
  dt_FE = min_element( 1.0 / **k0[40][40] );

  if (dt <= dt_FE) {
    //BuildCollisionMatrix(R_In, R_Out, Nold, dV, N_G, 0.0);
	BuildCollisionMatrix;
    cMat[40][40] = get_cMat();
    kVec[40] = get_kVec();


  *Nnew[N_G] = *Nold[N_G] + dt * ( **cMat[40][40] * *Nold[N_G] );

  nTrueIterations = nTrueIterations + 1;

  }

  else {
    double r0 = 1.0 / (**k0[40][40] * dt);      
         
         double Alpha0 = (pow(160 *r0, 3) + pow(60 * r0,2) + (11 * r0) + 1) /(pow(360 * r0,3) + pow(60 * r0,2) + (12 * r0) + 1);
     
         double Np = *Nold[N_G] + dt * ( **F0[40][40] - **k0[40][40] * *Nold[N_G] ) / (1.0 + (Alpha0 * **k0[40][40] * dt));
         
         ComputeRates;
          
		double ** kp[40][40];
  		kp[40][40] = get_kp();
  		double ** Fp[40][40];
  		Fp[40][40] = get_Fp();

         double kBAR = 0.5 * ( ** kp[40][40] + **k0[40][40] );
      
         double rBAR = 1.0 / (kBAR * dt);
      
         double AlphaBAR = (pow(160 * rBAR,3) + pow(60 * rBAR,2) + (11 * rBAR) + 1) /(pow(360 * rBAR,3) + pow(60 * rBAR,2) + (12 * rBAR) + 1);
            
         double Ft = AlphaBAR * ** Fp[40][40] + ( 1.0 - AlphaBAR ) * **F0[40][40];
    
      
         *Nnew[N_G] = *Nold[N_G] + dt * (Ft - kBAR * *Nold[N_G]) / (1.0 + (Alpha0 * kBAR * dt));
          
         nTrueIterations = nTrueIterations + 2;
         
      }
      
      if( abs( sum(*Nnew[N_G], *dV[40]) - sum(*Nold[N_G], *dV[40]) ) / sum(*Nold[N_G], *dV[40]) > tolC ) {
        reStep = true;
      }
  
	}

 /*case 6: {

  Branch = BE;

  //void BuildCollisionMatrix(R_In, R_Out, Nold, dV, N_G, 0.0);
  BuildCollisionMatrix;
  cMat[40][40] = get_cMat();
  kVec[40] = get_kVec();


   dt_FE = min( 1.0 / *kVec[40] );
        
    if( reStep ) {
       dt = dt_dec * dt;
               reStep = false;
   }
    else {
        dt = dt_grw * dt;
       }

    if( dt <= dt_FE ) {
         
          *Nnew[N_G] = *Nold[N_G] + dt * ( **cMat[40][4] * *Nold[N_G] );
        
          nTrueIterations = nTrueIterations + 1;
      }

    else 
  //  [Nnew, nIterations] = NewtonRaphson( Nold, dt, R_In_H, R_Out_H, N_g, tolBE );
	NewtonRaphson;
  		


   nTrueIterations = nTrueIterations + nIterations;

    if ( max( abs( *Nnew[N_G] - *Nold[N_G] ) / max( *Nold[N_G], 1.0e-8 ) ) > tolN ) {
            reStep = true; 

    }

  }
  */
  default: 

   Branch = FE;
   //void BuildCollisionMatrix(R_In, R_Out, Nold, dV, N_G, tolPE);
	BuildCollisionMatrix;
    cMat[40][40] = get_cMat();
    kVec[40] = get_kVec();

   dt = min( min_element( 1.0 / *kVec[40] ), dt_grw * dt );
   *Nnew[N_G] = *Nold[N_G] + dt * ( **cMat[40][40] * *Nold[N_G] );
  



  }


  if  ( reStep != true ) {

  Nold[N_G] = Nnew[N_G];
  t    = t + dt;
  }

  if (t >= t_end || cycle >= cycleM) {

  done = true;
  }

/*
  if ( mod(cycle, cycleD) == 1) {
  disp* std::fprintf( '  Cycle = %d, t = %d, dt = %d ', cycle, t, dt ) );
  }
  */

  //if  ( t >= wrtTimes(wrtCount) ) {
	if  ( t >= wrtTimes ) {
// Write data file

  //PlotFileNumber = Write_Plotfile( t, dt);
  PlotFileNumber = Write_Plotfile();

  wrtCount = wrtCount  + 1;

  


  Utilities::stopTimer();
  	}

	}

}

//********************************************************


// Initializing elements of matrix mult to 0.
//    for(i = 0; i < r1; ++i)
//        for(j = 0; j < c2; ++j)
//       {
//           mult[i][j]=0;
//        }
//
    // Multiplying matrix a and b and storing in array mult.
 //   for(i = 0; i < r1; ++i)
 //       for(j = 0; j < c2; ++j)
 //           for(k = 0; k < c1; ++k)
 //           {
 //               mult[i][j] += a[i][k] * b[k][j];
 //          }






// NewtonRaphson

/*
void NewtonRaphson() {

 *Nnew = *Nold;

bool converged   = false;
  int nIterations = 0;
  while (( converged != true ) )

     dN = ( eye( N_G ) - dt * JAC_FD( *Nnew, R_In, R_Out, N_G ) ) / ( ( Nold - Nnew ) + dt * RHS_FD( *Nnew, R_In, R_Out, N_G ) );
    
    *Nnew = Nnew + dN;
    
    nIterations = nIterations + 1;
    
    if( norm( dN / Nnew ) < Tol )
      converged = true;
}

// JAC_FD

void JAC_FD() {
**R_In[40][40]( eye( size( **R_In[40][40]  ) ) ) = 0.0;
**R_Out[40][40]( eye( size( **R_Out[40][40] ) ) ) = 0.0;


//R_In(logical( eye( size( R_In  ) ) )) = 0.0;
//R_Out(logical( eye( size( R_Out ) ) )) = 0.0;



}


// RHS_FD
void RHS_FD() {
	double *RHS = ( ones( N_G, 1 ) - N ) * ( **R_In[40][40] * *N[40] ) - N * ( R_Out * ( ones( N_G, 1 ) - N ) );
}

void eye( int Size){
	int i = 0;
	int m = Size;
	int j = 0;
	int n = Size;
	int IdentityMatrix[Size][Size]; 
	for (i = 0; i < m; i++) {
		 for (j = 0; j < n; j++) {
			if (i=j){
				IdentityMatrix[i][j] = 1;
			}
			else if (i != j){
				IdentityMatrix[i][j] = 0;
			}
		 } 
		 }
}
    
*/