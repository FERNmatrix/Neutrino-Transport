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
#include <stdio.h>
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
using namespace std;
 using std::string;

#define N_G 40
clock_t startCPU, stopCPU;
#define START_CPU if ((startCPU=clock())==-1) {printf("Error calling clock"); exit(1);}
#define STOP_CPU if ((stopCPU=clock())==-1) {printf("Error calling clock"); exit(1);}
#define PRINT_CPU (printf("Timer: %7.4e ms used", 1000*(double)(stopCPU-startCPU)/CLOCKS_PER_SEC));
#define FPRINTF_CPU (fprintf(pFile, "in %7.4e seconds\n", (double)(stopCPU-startCPU)/CLOCKS_PER_SEC));
#define FPRINTF_CPU2 (fprintf(pFile2, "in %7.4e seconds\n", (double)(stopCPU-startCPU)/CLOCKS_PER_SEC));
#define FPRINTF_CPUD (fprintf(pFileD, "in %g seconds\n", (double)(stopCPU-startCPU)/CLOCKS_PER_SEC));
#define PRINT_CPU_TEST (printf("\nTimer Test: %g ms used by CPU\n", 1000*(double)(stopCPU-startCPU)/CLOCKS_PER_SEC));
// Global Variables

int model = 1;
// int N_G = 40;
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

double G_A = 7.5e-01;    
double G_B = 1.0e+02;    
double G_C = sqrt(50.0);

double cycleM = 10e9;    
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

char scheme[] = "ExplicitAsymptotic";

char Comment = ' ';
char PlotFileDir[] = "Output";
char PlotFileName[] = "PlotFile";
int PlotFileNumber = 0;
int nPlotFiles = 100; 
double ** tempA[N_G][N_G];
double * tempKvec[N_G];
double ** tempF[N_G][N_G];
double ** tempKmatrix[N_G][N_G];


// Function signatures in main:


 // class defenitions 

 class Utilities{
  private:
  public:
    double eC[N_G];
    double de[N_G];

    void initalizeNES();
    double * ReadData1D(char FileName[], double N);
    double ** ReadData2D(char FileName[], double M, double N);
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



void Utilities::initalizeNES(){
      double * eC = ReadData1D("Data/NES_RATES_EnergyBinCenter.dat",N_G);
      double * de = ReadData1D("Data/NES_RATES_EnergyBinWidth.dat",N_G);
      double dV[N_G];

      double N_Eq[N_G];
    for (int i = 0; i <= (N_G-1); ++i){
    
    dV[i]= (pow(eC[i]+0.5e0*de[i],3) - pow(eC[i]-0.5e0*de[i],3))/3.0e0;


      }
 //     double *R_In;
      double **R_In = ReadData2D("Data/NES_RATES_R_In___001.dat",N_G,N_G);

//      double *R_Out;

//      R_Out = conj(R_In);
//      N_Eq =1.0/ (exp(int gsl_vector_div(gsl_vector *(eC-mu), const gsl_vector *kt))+1.0);


    int gsl_matrix_complex_conjtrans_memcpy(gsl_matrix *R_Out, const gsl_matrix *R_In); // This function makes R_Out the conjugate transpose of R_In 
//      int k = 0;
      switch(model) {
        case 1:
          double mu = 145.254;
          double kt = 20.5399;
          for (int i = 1; i <= N_G; ++i){
             N_Eq[N_G] = 1 /  (exp( (eC[i]-mu)/kt ) + 1);
          }
          for (int j = 1; j <= N_G; ++j){
            for (int i = 1; i <= N_G; ++i){
              if( j < i )
                R_In[j][i] = R_In[j][i] * exp( ( eC[j] - eC[i] ) / kt );
            }
          }

        break;

        int gsl_matrix_complex_conjtrans_memcpy(gsl_matrix *R_Out, const gsl_matrix *R_In);

        case 2: 
           mu = 045.835;
           kt = 15.9751;
          for (int i = 1; i <= N_G; ++i){
             N_Eq[N_G] = 1 /  (exp( (eC[i]-mu)/kt ) + 1);
          }
          for (int j = 1; j <= N_G; ++j){
            for (int i = 1; i <= N_G; ++i){
              if( j < i )
                R_In[j][i] = R_In[j][i] * exp( ( eC[j] - eC[i] ) / kt );
            }
          }

        break;

        int gsl_matrix_complex_conjtrans_memcpy(gsl_matrix *R_Out, const gsl_matrix *R_In);

        case 3:
           mu = 020.183;
           kt = 07.7141;
          for (int i = 1; i <= N_G; ++i){
             N_Eq[N_G] = 1 /  (exp( (eC[i]-mu)/kt ) + 1);
          }
          for (int j = 1; j <= N_G; ++j){
            for (int i = 1; i <= N_G; ++i){
              if( j < i )
                R_In[j][i] = R_In[j][i] * exp( ( eC[j] - eC[i] ) / kt );
            }
          }

        break;

        int gsl_matrix_complex_conjtrans_memcpy(gsl_matrix *R_Out, const gsl_matrix *R_In);

        case 4:
           mu = 009.118;
           kt = 07.5830;
          for (int i = 1; i <= N_G; ++i){
             N_Eq[N_G] = 1 /  (exp( (eC[i]-mu)/kt ) + 1);
          }
          for (int j = 1; j <= N_G; ++j){
            for (int i = 1; i <= N_G; ++i){
              if( j < i )
                R_In[j][i] = R_In[j][i]  * exp( ( eC[j] - eC[i] ) / kt );
            }
          }

        break;

        int gsl_matrix_complex_conjtrans_memcpy(gsl_matrix *R_Out, const gsl_matrix *R_In);
        case 5:
           mu = 003.886;
           kt = 03.1448;
          for (int i = 1; i <= N_G; ++i){
             N_Eq[N_G] = 1 /  (exp( (eC[i]-mu)/kt ) + 1);
          }
          for (int j = 1; j <= N_G; ++j){
            for (int i = 1; i <= N_G; ++i){
              if( j < i )
                R_In[j][i] = R_In[j][i]  * exp( ( eC[j] - eC[i] ) / kt );
              
            }
          }

        break;

        int gsl_matrix_complex_conjtrans_memcpy(gsl_matrix *R_Out, const gsl_matrix *R_In);
        default: N_Eq[N_G] = 1.0;
      
      }
}
  
  


double * Utilities::ReadData1D(char FileName[], double N) {
 FILE *fr;
 fr = fopen(FileName, "r");
  double Data1D[40];

  while (!feof(fr)){

    fread(Data1D, sizeof(Data1D), 1, fr);
    cout << Data1D;
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
    cout << ReadData2D;

    int flat[100];
    for (int i=0; i<100; ++i)
       flat[i]=i;
         int (&square)[10][10] = *reinterpret_cast<int(*)[10][10]>(flat);
   for (int i=0; i<10; ++i)
        cout << square[5][i] << ' ';
   cout << '\n';
   fclose (fr);
   return ** Data2D;

}

// Main

int main() {
Utilities::startTimer();

  FILE *pFile;
  pFile = fopen("output/NT.data","w");
  int totalTimeSteps = 0;
  Logfile(stop_time, t, t_W0, dt, G_A, G_B, G_C, tolC, scheme, Comment);

fprintf(pFile, "# data_output \n");
fprintf(pFile, "#       t        dt  \n");

  initializeNES( Model, N_g);

  eC = get_eC();
  dV = get_dV();
  R_In[][] = get_R_In();
  R_Out[][] = get_R_Out();
  N_Eq = get_N_Eq();

  double multi[40][40];

  for (int i = 0; i < N_G; ++i){
    for (int j = 0; j < N_G; ++j){
      if (i = j) a[i][j] = dV[i];
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
        multi[i][j] += R_In[i][k] * a[k][j] 

      }
    }
  }

R_In_H = multi

  for (int i = 0; i < N_G; ++i){
    for (int j = 0; j < N_G; ++j){
      multi[i][j] = 0;
    }
  }

  for (int i = 0; i < N_G; ++i){
    for (int j = 0; j < N_G; ++j){
      for(int k = 0; k < N_G; ++k){
        multi[i][j] += R_Out[i][k] * a[k][j] 

      }
    }
  }

R_Out_H = multi

if (reStart == true) {
  end;
}

else {

  for (int i = 0; i < N_G; ++i) {

 N_0[i] = G_A * exp( - 0.5 * pow(( ( eC[i] - G_B ) / G_C ),2));
  
  }
}





for (int i = 0; i < N_G; ++i){
 Nold = N_0[i];
}



void BuildCollisionMatrix(R_In, R_Out, Nold, dV, N_G, 0.0);

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


PlotFileNumber = Write_Plotfile ( t, dt, Nold, eC, dV, kVec, cMat, 0, 0, 0, 0, FE, dt, dt, dt, PlotFileDir, PlotFileName, PlotFileNumber );

wrtTimes = logspace( log10(t_W0), log10(t_end), nPlotFiles );
wrtCount = 1;

while (done != true) {

 true_cycle = true_cycle + 1; 

  if not( reStep ) {

  cycle = cycle +1;

  }



  Nold = ApplyPerturbation(Nold, amp dV, N_G);

  switch Scheme {

  case 'ExplicitAsymptotic'

  Branch = EA;

void BuildCollisionMatrix(R_In, R_Out, Nold, dV, N_G, 0.0);

cMat[40][40] = get_cMat();
kVec[40] = get_kVec();

  dt_FE = min(1.0 / kVec);

  if reStep {

    dt = (dt_dec * dt);
    reStep = false;

    else 
      dt = dt_grw * dt;

  }

  if dt <= dt_FE {

    Nnew = Nold + dt * (cMat * Nold);
}
  else {

    Nnew = Nold + (dt / (1.0 + kVec *dt)) * (cMat * Nold);
}

if( abs( sum(Nnew * dV) - sum(Nold * dV) ) / sum(Nold * dV) > tolC ) {

  reStep = true;
}

if ( max( abs( Nnew - Nold ) / max( Nold, 1.0d-8 ) ) > tolN ) {

  reStep = true;

}


  case 'PartialEquilibrium'

  Branch = FE_PE;

  void BuildCollisionMatrix(R_In, R_Out, Nold, dV, N_G, tolPE);

  cMat[40][40] = get_cMat();
  kVec[40] = get_kVec();


  dt = min( min( 1.0 / kVec ), dt_grw * dt);

  Nnew = Nold + dt * (cMat * Nold);

  case 'QuasiSteadyState'

  Branch = QSS1;

  if reStep {
    dt = dt_dec * dt;
    reStep = false;
  }

    else
      dt = dt_grw * dt;
  

  void ComputeRates( R_In, R_Out, Nold, dV);

  double ** F0[40][40];
  F0[40][40] = get_F0();
  double ** k0[40][40];
  k0[40][40] = get_k0();
  
  exp_kdt = exp( - dt * k0 );
      
    dt_FE = min( 1.0 / k0 );
      
    if( dt <= dt_FE ) {

    void BuildCollisionMatrix(R_In, R_Out, Nold, dV, N_G, 0.0); 

    cMat[40][40] = get_cMat();
    kVec[40] = get_kVec();


    Nnew = Nold + dt * ( cMat * Nold );

    nTrueIterations = nTrueIterations + 1 ;

  }

    else {
      void BuildCollisionMatrix(R_In, R_Out, Nold, dV, N_G, 0.0);

      cMat[40][40] = get_cMat();
      kVec[40] = get_kVec();


      Nnew = Nold + ( cMat * Nold ) * ( 1.0 - exp_kdt ) / k0;
    }

    if ( abs( sum(Nnew * dV) - sum(Nold * dV) ) / sum(Nold * dV) > tolC ) {
      reStep = true;

      }

    case 'QSScorrection'

    Branch = QSS2;

    if reStep {
      dt = dt_dec * dt; 
        reStep = false;
      else
        dt = dt_grw * dt;
  }

    void ComputeRates( R_In, R_Out, Nold, dV );
    F0
    k0
      
  dt_FE = min( 1.0 ./ k0 );

  if dt <= dt_FE {
    void BuildCollisionMatrix(R_In, R_Out, Nold, dV, N_G, 0.0);

    cMat[40][40] = get_cMat();
    kVec[40] = get_kVec();


  Nnew = Nold + dt * ( cMat * Nold );

  nTrueIterations = nTrueIterations + 1;

  }

  else {
    r0 = 1.0 / (k0 * dt);      
         
         Alpha0 = ((160 * r0^3) + (60 * r0^2) + (11 * r0) + 1) /((360 * r0^3) + (60 * r0^2) + (12 * r0) + 1);
     
         Np = Nold + dt * ( F0 - k0 * Nold ) / (1.0 + (Alpha0 * k0 * dt));
         
         [Fp, kp] = ComputeRates( R_In, R_Out, Np, dV );
      
         kBAR = 0.5 * ( kp + k0 );
      
         rBAR = 1.0 / (kBAR * dt);
      
         AlphaBAR = ((160 * rBAR^3) + (60 * rBAR^2) + (11 * rBAR) + 1) /((360 * rBAR^3) + (60 * rBAR^2) + (12 * rBAR) + 1);
            
         Ft = AlphaBAR * Fp + ( 1.0 - AlphaBAR ) * F0;
    
      
         Nnew = Nold + dt * (Ft - kBAR * Nold) / (1.0 + (Alpha0 * kBAR * dt));
          
         nTrueIterations = nTrueIterations + 2 ;
         
      }
      
      if( abs( sum(Nnew * dV) - sum(Nold * dV) ) / sum(Nold * dV) > tolC ) {
        reStep = true;
      }
  }


  case 'BackwardEuler'

  Branch = BE;

  void BuildCollisionMatrix(R_In, R_Out, Nold, dV, N_G, 0.0);

  cMat[40][40] = get_cMat();
  kVec[40] = get_kVec();


   dt_FE = min( 1.0 / kVec );
        
    if( reStep ) {
       dt = dt_dec * dt;
               reStep = false;
   }
    else {
        dt = dt_grw * dt;
       }

    if( dt <= dt_FE ) {
         
          Nnew = Nold + dt * ( cMat * Nold );
        
          nTrueIterations = nTrueIterations + 1;
      }

    else 
     [Nnew, nIterations] = NewtonRaphson( Nold, dt, R_In_H, R_Out_H, N_g, tolBE );

   nTrueIterations = nTrueIterations + nIterations;

    if ( max( abs( Nnew - Nold ) / max( Nold, 1.0d-8 ) ) > tolN ) {
            reStep = true; 

    }


  default: 

   Branch = FE;
   void BuildCollisionMatrix(R_In, R_Out, Nold, dV, N_G, tolPE);
    cMat[40][40] = get_cMat();
    kVec[40] = get_kVec();

   dt = min( min( 1.0 ./ kVec ), dt_grw * dt );
   Nnew = Nold + dt * ( cMat * Nold );
  

}



if  not( reStep ) {

  Nold = Nnew;
  t    = t + dt;
}

if t >= t_end || cycle >= cycleM {

  done = true;
}

if  if mod(cycle, cycleD) == 1 {
  disp( fprintf( '  Cycle = %d, t = %d, dt = %d ', cycle, t, dt ) );
}

if  t >= wrtTimes(wrtCount) {

// Write data file

  PlotFileNumber = Write_Plotfile( t, dt, Nold, eC, dV, kVec, cMat, cycle, true_cycle, nIterations, nTrueIterations, Branch, dt_FE, dt_EA, dt_PE, PlotFileDir, PlotFileName, PlotFileNumber);

wrtCount = wrtCount  + 1;

  }

}
Utilities::stopTimer();

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
// Write_Plotfile

int Write_Plotfile( double t, double dt )


fprintf(pFile, "%7.4f %7.4f", t, dt)
FileNumber = FileNumber + 1;
return FileNumber;


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

// Logfile

  void logfile(t_end, t, t_W0, dt, G_A, G_B, G_C, tolC, Scheme, Comment) 
  FILE *Log;
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

// BuildCollisionMatrix 

void BuildCollisionMatrix( double R_In[], double R_Out[], double N, double dV, int N_G, double tol) {

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
        B = ( R_In[i][j] * dV[i] + R_Out[i][j] * dV[j] );
            C = dV[i] * N[i] + dV[j] * N[j];
        
            a = ( R_In[i][j] - R_Out[i][j] ) * dV[i];
            b = B + ( R_In[i][j] - R_Out[i][j] ) * C;
            c = R_In[i][j] * C;
            d = (b*b) - (4.0 * a * c);
        
            N_Eq_i = 0.5 * ( b - sqrt( d ) ) / a;
            N_Eq_j = ( C - N_Eq_i * dV[i] ) / dV[j];
      }
      else {
        N_Eq_i = N[i];
            N_Eq_j = N[j];
        }
      
    }
}


diff_i = abs( N_Eq_i - N[i] ) / max( [ 1.d-16 N_Eq_i ] );
    diff_j = abs( N_Eq_j - N[j] ) / max( [ 1.d-16 N_Eq_j ] );
    if( or( diff_i > tol, diff_j > tol ) )

      A[i][j] = A[i][j]+ (1.0 - N[i]) * R_In[i][j] * dV[j];
      A[i][j] = A[i][j] - R_Out[i][j] * dV[j] * (1.0 -  N[j]);
      k[i] = k[i] + R_Out[i][j] * dV[j] + ( (R_In[i][j] * dV[j] - R_Out[i][j] * dV[j] * N[j]));

}
// ComputeRates

void ComputeRates( double R_In[], double R_Out[], double N, double dV ) {

  F = R_In * ( N * dV );
    k = F + R_Out * ( ( 1.0 - N ) * dV );

  }

// NewtonRaphson

// void NewtonRaphson( double Nold, double dt, double R_In[], double R_Out[], int N_G, double Tol )


