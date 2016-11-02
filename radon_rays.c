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

int N;
char *in_file_ray, *in_file_ttime, *in_file_slowness, *out_file;
clock_t tini, tend, tacum;
double cpu_time_used;

FILE *out =NULL; 


//////////////////////////STRUCTURES//////////////////////////
struct read {
  
  double *dd; /* input data -> model data -size nt*nx- */
  double *mm; /* input data -> travel time -data size nt*np- */
  
  double *slowness; /*input data -> the initial model */
  double *ds; /* input data -> ray trajectory -data size nt*nx- */
  

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

int radon(char *out_file, int N)
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
  nt = N;
  /* number of samples in time axis */
 
  
  if (adj||inv)
    { /* m(tau,p)=sum_{i=0}^{nx} d(t=tau+p*x_i,x_i) */
      
      nx = N;
      /* number of offset if the input in the data domain */
    
      /* specify slope axis */
      np = N;
      
      /* number of p values (if adj=y) */
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
  xx==(double *)malloc(nx*sizeof(double));
  cdd=(complex*)malloc(nw*nx*sizeof(complex));
  cmm=(complex*)malloc(nw*np*sizeof(complex));
  tmpr=(double*)fftw_malloc(nfft*sizeof(double));
  tmpc=(fftw_complex*)fftw_malloc(nw*sizeof(fftw_complex));
  fft1=fftw_plan_dft_r2c_1d(nfft,tmpr,tmpc,FFTW_MEASURE);	
  ifft1=fftw_plan_dft_c2r_1d(nfft,tmpc,tmpr,FFTW_MEASURE);
  

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
	  tmpr[ip] = data.mm[ip];
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
  
  exit(0);
} 

/*Reads output  -the RT data- and computes the 1D FT-. To obtain the wave velocity a 2D IFT is applied to before step*/

int fourier(int N)
{
  FILE *out_wave1D = NULL;
  FILE *out_wave1 = NULL;
  FILE *out_wave2D = NULL;
  
  int i, j, nread;
 
  
  fftw_plan my_plan1, my_plan2, my_plan3;
  fftw_complex *in_radon, *in,*out_origin, *in_fourier,*out_fourier1D, *out_fourier2D, *FT_RT;

  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  in_radon = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  out_origin = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  out_fourier1D = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  out_fourier2D = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N);
  FT_RT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);

  


  ///////////////////////////////////READING THE RADON DATA/////////////////////////
  
  tini = clock();
      
  for(i=0; i<=N; i++)
    {
  
      in_radon[i][0] = data.mm[i];  //REAL PART -travel time -> model or Radon domain-

      in_radon[i][1] = 0.0;          //IMAGINARY PART 
      
    }

  ////////////////////////////////Calculating 1D FFT OF RT-DATA /////////////////////////
  
  my_plan1 = fftw_plan_dft_1d(N, in_radon, out_fourier1D, FFTW_FORWARD, FFTW_MEASURE);
   
  fftw_execute(my_plan1);

  out_wave1D = fopen("1DFT_of_RT.dat","w");

  for(i=0; i<=N; i++)
    {
      fprintf(out_wave1D,"%lf %lf\n", out_fourier1D[i][0], out_fourier1D[i][1]);
    }


 ////////////////////////////////Calculating 2D FFT OF THE ORIGINAL DATA ///////////////////

  	
  for(i=0; i<=N; i++)
    {
  
      in[i][0] = data.dd[i];   //REAL PART -ray path/slowness -> data domain - 

      in[i][1] = 0.0;          //IMAGINARY PART 
      
    }
	
	
  my_plan2 = fftw_plan_dft_2d(N, N, in, out_origin, FFTW_FORWARD, FFTW_MEASURE);
   
  fftw_execute(my_plan2);

  out_wave1 = fopen("2DFT_ORIGINAL_DATA.dat","w");

  for(i=0; i<=N; i++)
    {
      fprintf(out_wave1,"%lf %lf\n", out_origin[i][0], out_origin[i][1]);
    }


  
////////////////////////////////Calculating 2D IFFT OF 1DFT_RT-DATA /////////////////////////
//HERE WE OBTAIN THE SLOWNESS//

  
  for(i=0; i<=N; i++)
    {
      FT_RT[i][0] =  out_fourier1D[i][0];
      
      FT_RT[i][1] =  out_fourier1D[i][1];
    }

  
  
  my_plan3 = fftw_plan_dft_2d(N, N, FT_RT, out_fourier2D, FFTW_BACKWARD, FFTW_MEASURE);

   
  fftw_execute(my_plan3);

  out_wave2D = fopen("2DIFT_of_1DFT_RT.dat","w");

  for(i=0; i<=N; i++)
    {
      fprintf(out_wave2D,"%lf %lf\n", out_fourier2D[i][0]/N*N, out_fourier2D[i][1]/N*N);
      
    }


   
  tend = clock();
  
  cpu_time_used = ((double) (tend - tini)) / CLOCKS_PER_SEC;
  
  printf("CPU TIME USED: %16.8lf\n",cpu_time_used);
  
  
 	
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

  int i;
  
  printf("%d\n",argc);
  
  if(argc != 6)
    {
      printf("ERROR--> use as:\n");
      printf("%s #iterations input_file_ray in_file_ttime in_file_model output_file\n",argv[0]);
      exit(0);  
    }
  
  N   = atoi(argv[1]);
  in_file_ray  = argv[2];
  in_file_ttime  = argv[3];
  in_file_slowness  = argv[4];
  out_file  = argv[5];
  
  printf("%s %d %s %s %s %s\n",argv[0], N, in_file_ray, in_file_ttime, in_file_slowness, out_file);
  
  
  data.dd = (double *) malloc(N* sizeof(double));              //data -> ds/c(x)
  data.mm = (double *) malloc(N* sizeof(double));              //travel time
  data.ds = (double *) malloc(N *sizeof(double));              //ds
  data.slowness = (double *) malloc(N *sizeof(double));        //slowness 1/c(x)
  
  /*Reading the file data - loading data*/
  
  read_file1(in_file_ray, N);

  read_file2(in_file_ttime, N);

  read_file3(in_file_slowness, N);

  for(i=0; i<N; i++)
    {
      data.dd[i] = (data.ds[i]) / (data.slowness[i]);
      
    }
  
  /*Calculate the Radon Transform to rays */
  
  radon(out_file, N);
  
  /*Calculate the 1D FT of RT and the 2D IFT of FT_RT -> COMPUTES CST: Central Slice Theorem*/
  
  fourier(N);
  
  return 0;
  
}
 









