/* Xibelly Eliseth Mosquera Escobar
  Implementation of the Radon Transform (RT) in the ray aproximation.

  The ray equation is:

   T(p)  = int_{S} dS / c(X)      (1)

   X -> (x,y) in 2D

 where T is the travel time, S is the ray path and c(X) is the velocity of the wave propagation.

 The aim is to obtain the medium's parameter -c- from the travel time measured on the surface and
 knowing the ray path. To this end the RT will be used, due to the form that this integral transform is.
 Also the Central Slice Theorem -CST-will be used, since this present a relation between the RT and the FT.
 So the IFT will be used to invert the equation (1).
 */

/*
ANALISIS AND DESIGN

To solve our problem we have to do:

-Receive the travel times, the model (1/c -slowness-) and the ray trajectories. 
Verify that these are apload in the correct way.

-Receive the iterations number (# of data lines).

-Do a fit to the trajectories -Nonlinear adjustment- to ensure the geometry of the trajectories

NOTE: AS A PRACTICAL EXERCISE, IN FIRST APPROXIMATION WE WILL CONSIDER THAT THE TRAJECTORIES 
 DESCRIBING THE CORRECT GEOMETRY, IS TO SAY, THESE ARE PARABOLAS. WITH THIS ASSUMPTION WE DON'T
 HAVE TO DO ANY ADJUSTMENT.

-Perform the RT to the trajectories.

-Use the CST to obtain c. This is done in the next way, 
 apply a 1D FT on the RT that was performed previously,
 and then apply a 2D IFT on the last set of data.

-Write in disk the results of the before step.

CONTROL TESTS

-the clock() function will be used to check the ejecution time

*/


/*
SET OF TESTS

In oder to understand the CST, we will compare the 2D FT of the data -ds/c- 
with the 1D FT of the  PRT of the data, is to say, the FT of the model -travel times- 
obtained trough the PRT.


 */


/* IMPLEMENTATION -> PARABOLIC RADON TRANSFORM */

/*Note: We are using the code Mmyradon2.c to implement the PRT in the ray approach*/


#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<fftw3.h>
#include<stdbool.h>
#include<complex.h>

#include"myradon2.c"


/*
PATHS: /home/xibelly/Madagascar/rsfsrc/user/pyang myradon2.c
PATHS: /home/xibelly/Madagascar/rsfsrc/system/seismic/Mradon.c radon.c
PATHS: /home/xibelly/Madagascar/include
 */

/*NOTA: el indicador sf_...  es para indicar el formato por defecto de MADAGASCAR -> rsf */

///////////////////////GLOBAL VARIABLES////////////////////

int gridx;
int gridz;
int NRAYS;
double GRADIENTE;
int MAX_ITERATIONS;
double Xini;
double Zini;
double Xfin; 
double Zfin;  
int r;

int num_ray;
int Ray_id, NSLOWNESS, Nds;
char *in_file_ray, *in_file_ttime, *in_file_slowness, *out_file;
clock_t tini, tend, tacum;
double cpu_time_used;

FILE *out =NULL; 
FILE *read =NULL; 

#include "Read_params.c"

//////////////////////////STRUCTURES//////////////////////////
struct read {
  
  double *dd; /* input data -> model data (ds/c) - */
  double *mm; /* input data -> travel time  */
  
  double *slowness; /*input data -> the initial model */
  double *ds;
  double *ds_x; /* input data -> ray trajectory -data size - */
  double *ds_z;

 
};

struct read data;



//////////////////////////SUB-STRUCTURES//////////////////////////

#include"input.c"



///////////////////////FUNCTIONS/////////////////////////////

int fft_next_fast_size(int n)
{
    while(1) {
        int m=n;
        while ( (m%2) == 0 ) m/=2;
        while ( (m%3) == 0 ) m/=3;
        while ( (m%5) == 0 ) m/=5;
        if (m<=1)
            break; /* n is completely factorable by twos, threes, and fives */
        n++;
    }
    return n;
}


void matrix_transpose(complex *matrix, int nx, int nz)
{
  int ix, iz;
  complex *tmp=(complex*)malloc(nx*nz*sizeof(complex));
  if (tmp==NULL) {printf("out of memory!\n"); exit(1);}
  for(iz=0; iz<nz; iz++)
    for(ix=0; ix<nx; ix++)
      tmp[iz+nz*ix]=matrix[ix+nx*iz];
  
  memcpy(matrix, tmp, nx*nz*sizeof(complex));
  free(tmp);
}

/*Computes the PRT to the load data -rays- */

int radon(char *out_file)  //Nds : # of ray points
{
  
  bool adj, inv, par;
  char *invmode;
  int iw, ip, ix, np, nt, nfft, nw, nx, niter;
  double dp, p0, dt, t0, dx, ox, x0, w, eps;
  double *p, *xx, *tmpr;             //**dd y **mm 
  complex *cdd, *cmm; 
  fftw_complex *tmpc;
  fftw_plan fft1, ifft1;

    
  adj=true;
  /* if y, perform adjoint operation */
  inv=true;
  /* if y, perform inverse operation */
  
  /* parameters */
  nt = Nds;
  /* number of samples in time axis, //# of ponist that constitues the ray */
 
  
  if (adj||inv)
    { /* m(tau,p)=sum_{i=0}^{nx} d(t=tau+p*x_i,x_i) */
      
      nx = Nds;  
      /* number of offset if the input in the data domain, //# of ponist that constitues the ray */
    
      /* specify slope axis */
      np = 1; //A ray
      
      /* number of p values (if adj=y) */

      dt = 1.0; 

      /*interval of axis time*/

      dp = 0.001;
    
      /* p sampling (if adj=y) */
      p0 = -0.2;
      
      /* p origin (if adj=y) */
      
      if(inv)
	{			
	  invmode="toeplitz";
	  /* inverse method: 'ls' if least-squares; 'toeplitz' -> 't' if use FFT */			
	  eps=0.01;  
	  /* regularization parameter */
	}
      else
	{
	  invmode=NULL;
	}
      
    }
  
    
  nfft=2*fft_next_fast_size(nt);
  nw=nfft/2+1;
  p=(double *)malloc(np*sizeof(double));
  xx=(double *)malloc(nx*sizeof(double));
  cdd=(complex*)malloc(nw*nx*sizeof(complex));
  cmm=(complex*)malloc(nw*np*sizeof(complex));
  tmpr=(double*)fftw_malloc(nfft*sizeof(double));
  tmpc=(fftw_complex*)fftw_malloc(nw*sizeof(fftw_complex));
  fft1=fftw_plan_dft_r2c_1d(nfft,tmpr,tmpc,FFTW_MEASURE);	
  ifft1=fftw_plan_dft_c2r_1d(nfft,tmpc,tmpr,FFTW_MEASURE);

  printf("STATE OF LOADING PARAMTERS FOR PRT IS: SUCESS\n");
  

  for(ip=0; ip<np; ip++) 
    p[ip]=p0+ip*dp;	
    
  ox = 0.0;
  /* data origin in x */
  dx = 0.025;
  
  for(ix=0; ix<nx; ix++) 
    xx[ix]=ox+ix*dx;
  
  par=true;
  /* if y, parabolic Radon transform */
  x0=1.0;   
  /* reference offset */
  
  for (ix=0; ix < nx; ix++)/* normalize offsets */
    {
      if (par) xx[ix] *= xx[ix]/(x0*x0);
      
    }
  
  if(adj||inv)
    {/* m(tau,p)=sum_{i=0}^{nx} d(t=tau+p*x_i,x_i) */
      for(ix=0; ix<nx; ix++) /* loop over offsets */
	{
	  memset(tmpr, 0, nfft*sizeof(double));
	  tmpr[ix] = data.dd[ix];	
	  fftw_execute(fft1);/* FFT: dd-->cdd */
	  memcpy(&cdd[ix*nw], tmpc, nw*sizeof(complex));
	}
      matrix_transpose(cdd, nw, nx);
      
      for(ip=0; ip<np; ip++) /* loop over slopes */
	{
	  memset(tmpr, 0, nfft*sizeof(double));
	  tmpr[ip] = data.mm[num_ray];               //num_ray
	  fftw_execute(fft1);/* FFT: mm-->cmm */
	  memcpy(&cmm[ip*nw], tmpc, nw*sizeof(double));
	  
	}
      matrix_transpose(cmm, nw, np);
    }
  
  
  myradon2_init(np, nx, dp, p, xx);
  for(iw=0; iw<nw; iw++) 
    {
      w=2.*M_PI*iw/(nfft*dt);
      myradon2_set(w);
      myradon2_lop(adj, true, np, nx, &cmm[iw*np], &cdd[iw*nx]);//Linear operator RT
      if(adj&&inv)
	{
	  if (invmode[0]=='t')
	    myradon2_inv(&cmm[iw*np], &cmm[iw*np], eps);//Inverts the model
	  
	}
    }
  
  for(ip=0; ip<np; ip++) // loop over slopes // 
    {			
      memcpy(tmpc, &cmm[ip*nw], nw*sizeof(complex));
      fftw_execute(ifft1); // IFFT: cmm-->mm //
      data.mm[ip]=tmpr[ip]/nfft;
      
    }
  
  //wrting in disk the output ->the PRT

  
  out = fopen(out_file,"w");

  if(out==NULL)
    printf("THE FILE CAN NOT BE CREATED\n");
  

  for(ip=0; ip<np; ip++) 
    {			
      fprintf(out,"%lf\n", data.mm[ip]);
    }
  
   
  
  free(p);
  free(xx);
  free(cdd); 
  free(cmm); 
  fftw_free(tmpr);
  fftw_free(tmpc);
  fftw_destroy_plan(fft1);
  fftw_destroy_plan(ifft1);

  printf("THE STATE OF THE PARABOLIC RADON TRANSFORM IS: SUCESS\n");
  
  return 0;
} 

/*Reads output  -the RT data- and computes the 1D FT-. To obtain the wave velocity a 2D IFT is applied to before step*/

int fourier()
{
  FILE *out_wave1D = NULL;
  FILE *out_wave1 = NULL;
  FILE *out_wave2D = NULL;

  char buff1[200], buff2[200], buff3[200];
  
  int i, j, nread;
 
  int N = 1; //A ray
  
  fftw_plan my_plan1, my_plan2, my_plan3;
  fftw_complex *in_radon, *in,*out_origin, *in_fourier,*out_fourier1D, *out_fourier2D, *FT_RT;

 
  in_radon = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  out_fourier1D = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  

  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nds*Nds);
  out_origin = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nds*Nds);
   
  FT_RT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  out_fourier2D = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N);

  ///////////////////////////////////READING THE RADON DATA/////////////////////////
  
       
  for(i=0; i<N; i++)
    {
  
      in_radon[i][0] = data.mm[i];  //REAL PART -travel time -> model or Radon domain-

      in_radon[i][1] = 0.0;          //IMAGINARY PART 
      
      printf("%lf\n",in_radon[i][0]);
    }

  ////////////////////////////////Calculating 1D FFT OF RT-DATA /////////////////////////
  
  my_plan1 = fftw_plan_dft_1d(N, in_radon, out_fourier1D, FFTW_FORWARD, FFTW_MEASURE);
   
  fftw_execute(my_plan1);
    
  sprintf(buff1,"%s_%d","1DFT_of_RT",Ray_id);
  
  out_wave1D = fopen(buff1,"w");

  if(out_wave1D==NULL)
    printf("THE FILE CAN NOT BE CREATED\n");

  for(i=0; i<N; i++)
    {
      fprintf(out_wave1D,"%lf %lf\n", out_fourier1D[i][0], out_fourier1D[i][1]);
    }

  printf("THE STATE OF 1D FT OVER RT IS: SUCESS\n");

 ////////////////////////////////Calculating 2D FFT OF THE ORIGINAL DATA ///////////////////
 
  	
  for(i=0; i<Nds; i++)
    {
      for(j=0; j<Nds; j++)
	{
	  
	  in[Nds * i + j][0] = data.dd[Nds * i + j];   //REAL PART -ray path/slowness -> data domain - 
	  
	  in[Nds * i + j][1] = 0.0;          //IMAGINARY PART 
	}
    }
	
	
  my_plan2 = fftw_plan_dft_2d(Nds, Nds, in, out_origin, FFTW_FORWARD, FFTW_MEASURE);
   
  fftw_execute(my_plan2);
  
  sprintf(buff2,"%s_%d","2DFT_ORIGINAL_DATA",Ray_id);

  out_wave1 = fopen(buff2,"w");

  if(out_wave1==NULL)
    printf("THE FILE CAN NOT BE CREATED\n");

  for(i=0; i<Nds; i++)
    {
      for(j=0; j<Nds; j++)
	{
	  fprintf(out_wave1,"%lf %lf\n", out_origin[Nds * i + j][0], out_origin[Nds * i + j][1]);
	}

    }
  printf("THE STATE OF 2D FT OVER THE DATA ds/c IS: SUCESS\n");
  
////////////////////////////////Calculating 2D IFFT OF 1DFT_RT-DATA /////////////////////////
  
//HERE WE OBTAIN THE SLOWNESS RELATED TO THE TRAVEL TIMES AND RAY PATHS THAT ARE LOAD//

  
  for(i=0; i<N; i++)
    {
      for(j=0; j<N; j++)
	{
	  
	  FT_RT[N* i + j][0] =  out_fourier1D[N* i + j][0];
	  
	  FT_RT[N* i + j][1] =  out_fourier1D[N* i + j][1];
	}
    }
  
  
  my_plan3 = fftw_plan_dft_2d(N, N, FT_RT, out_fourier2D, FFTW_BACKWARD, FFTW_MEASURE);

   
  fftw_execute(my_plan3);


  sprintf(buff3,"%s_%d","2DIFT_of_1DFT_RT",Ray_id);
  
  out_wave2D = fopen(buff3,"w");

  if(out_wave2D==NULL)
    printf("THE FILE CAN NOT BE CREATED\n");

  for(i=0; i<N; i++)
    {
      for(j=0; j<N; j++)
	{
	  fprintf(out_wave2D,"%lf %lf\n", out_fourier2D[N* i + j][0]/N*N, out_fourier2D[N* i + j][1]/N*N);
	}
    }

  printf("THE STATE OF 2D IFT OVER 1D FT OF RT IS: SUCESS\n");

  /////////////////////////////////////////////////////////////////////////////////////////
  
   	
  fclose(out_wave1);
  fclose(out_wave1D);
  fclose(out_wave2D);


  fftw_destroy_plan(my_plan1);
  fftw_destroy_plan(my_plan2);
  fftw_destroy_plan(my_plan3);

  fftw_free(in);
  fftw_free(in_radon);
  fftw_free(out_origin);
  fftw_free(out_fourier1D);
  fftw_free(out_fourier2D);
  fftw_free(FT_RT);

  return 0;
}


//////////////////////////////////MAIN PROGRAM/////////////////////
int main(int argc, char **argv){

  int i,j, ch;
    
  printf("%d\n",argc);
  
  if(argc != 7)
    {
      printf("ERROR--> use as:\n");
      printf("%s param_file input_file_ray in_file_ttime in_file_model output_file Ray_id\n",argv[0]);
      exit(0);  
    }
  
  in_file_ray  = argv[2];
  in_file_ttime  = argv[3];
  in_file_slowness  = argv[4];
  out_file  = argv[5];
  Ray_id = atoi(argv[6]); //Ray ID 
  
  printf("%s %s %s %s %s %s %d\n",argv[0], argv[1], in_file_ray, in_file_ttime, in_file_slowness, out_file, Ray_id);
 

  /*Loading the parameters data*/

  tini = clock();

  Param_SPM(argv[1]);

  NSLOWNESS = gridx*gridz;

  printf("%d %d %d %lf %d %lf %lf %lf %lf %d\n", gridx, gridz, NRAYS, GRADIENTE, MAX_ITERATIONS, Xini, Zini, Xfin, Zfin, r);

  //---------------------------------------------------------Number of lines of in_file_ray

  read = fopen(in_file_ray,"r"); 

  if(read==NULL)
    printf("THE FILE CAN NOT BE OPENED\n");

  Nds = 0;

  while((ch = fgetc(read)) != EOF)
    if(ch == '\n')
      Nds++;
  
  printf("THE INPUT FILE %s HAS %d LINES\n",in_file_ray, Nds);
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  data.dd = (double *) malloc(Nds* sizeof(double));              //data -> ds/c(x)
  data.ds_x = (double *) malloc(Nds*sizeof(double));              //ds
  data.ds_z = (double *) malloc(Nds*sizeof(double));
  data.ds = (double *) malloc(Nds*sizeof(double));
  data.slowness = (double *) malloc(NSLOWNESS *sizeof(double));  //slowness 1/c(x)
  data.mm = (double *) malloc(NRAYS* sizeof(double));         //travel time
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  
  /*Reading the file data - loading data*/
  
  read_file1(in_file_ray, Nds); //READ RAYS' COORDINATES

  read_file2(in_file_ttime, NRAYS);  //READ FILE TRAVEL TIME

  read_file3(in_file_slowness, NSLOWNESS);  //READ FILE SLOWNESS

  for(i=0; i<Nds; i++)
    {
      for(j=0; j<NSLOWNESS; j++)
	{
	  data.dd[i] = (data.ds[i]) / (data.slowness[j]); //DUDA!!!
	  
	}
      
    }
  
  /*Calculate the Radon Transform to rays */

  num_ray = Ray_id -1;

  printf("THE SELECTED RAY IS :%d\n", Ray_id);
  
  printf("THE TRAVEL TIME OF THIS RAY IS:\n");
  printf("%lf\n",data.mm[num_ray]);
  
  radon(out_file);
  
  /*Calculate the 1D FT of RT and the 2D IFT of FT_RT -> COMPUTES CST: Central Slice Theorem*/
  
  fourier();

  tend = clock();
  
  cpu_time_used = ((double) (tend - tini)) / CLOCKS_PER_SEC;
  
  printf("CPU TIME USED: %16.8lf\n",cpu_time_used);
  
  return 0;
  
}
 









