/* Xibelly Eliseth Mosquera Escobar
Implementation of the Radon Transform (RT) in the ray aproximation.

The ray equation is:

   T(p)  = int_{S} dS / c(X)      (1)

X -> (x,y) in 2D

where T is the travel time, S is the ray path and c(X) is the velocity of the wave propagation.

The aim is to obtain the medium's parameter -c- from the travel time measured on the surface and knowing the ray path. To this end the RT will be used, due to the form that this integral transform is. Also the Central Slice Theorem -CST-will be used, since this present a relation between the RT and the FT. So the IFT will be used to invert the equation (1).
 */

/*
ANALISIS AND DESIGN

To solve our problem we have to do:

-Receive the travel times and the ray trajectories, and verify that these are apload in the correct way

-Receive the iterations number (# of data lines)

-Do a fit to the trajectories -Nonlinear adjustment- to ensure the geometry of the trajectories

NOTE: AS A PRACTICAL EXERCISE, IN FIRST APPROXIMATION WE WILL CONSIDER THAT THE TRAJECTORIES DESCRIBING THE CORRECT GEOMETRY, IS TO SAY, THESE ARE PARABOLAS. WITH THIS ASSUMPTION WE DON'T HAVE TO DO ANY ADJUSTMENT.

-Perform the RT to the trajectories

-Use the CST to obtain c. This is done in the next way, apply a 1D FT on the RT that was performed previously, and then apply a 2D IFT on the last set of data.

-Write in disk the results of the before step.

CONTROL TESTS

-the clock() function will be used to check the ejecution time

*/

/* IMPLEMENTATION -> PARABOLIC RADON TRANSFORM */

/*Note: The idea is to use the code Mmyradon2.c to implement the PRT in the ray approach*/



#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
//#include<complex.h>

#include<fftw3.h>
#include"rsf.h"
#include"myradon2.h"

/*PATHS: /home/xibelly/Madagascar/rsfsrc/user/pyang myradon2.c
PATHS: /home/xibelly/Madagascar/rsfsrc/system/seismic/Mradon.c radon.c
PATHS: /home/xibelly/Madagascar/include */


/*LIBRARY PATHS

  /home/xibelly/Madagascar/rsf/include -> source rsf.h

  /home/xibelly/Madagascar/rsf/lib    -> file.o of rsf.h

  /home/xibelly/Madagascar/rsfsrc/build/api/c  -> file.o of rsf.h

*/

/*NOTA: el indicador sf_...  es para indicar el formato por defecto de MADAGASCAR -> rsf */

///////////////////////GLOBAL VARIABLES////////////////////

int N;
char *in_file_ray, *in_file_ttime, *out_file;
clock_t tini, tend, tacum;
double cpu_time_used;




//////////////////////////STRUCTURES//////////////////////////

struct ray {
  
  double *ttime;
  double *path;

};

struct ray RAY;



///////////////////////FUNCTIONS/////////////////////////////
void matrix_transpose(sf_complex *matrix, int nx, int nz)
{
  int ix, iz;
  sf_complex *tmp=(sf_complex*)malloc(nx*nz*sizeof(sf_complex));
  if (tmp==NULL) {printf("out of memory!\n"); exit(1);}
  for(iz=0; iz<nz; iz++)
    for(ix=0; ix<nx; ix++)
      tmp[iz+nz*ix]=matrix[ix+nx*iz];
  
  memcpy(matrix, tmp, nx*nz*sizeof(sf_complex));
  free(tmp);
}

/*Computes the LRT or PRT to the load data -rays- */

int radon(char *in_file_ray, char *in_file_ttime, char *out_file)
{
  sf_file in, in2, out;
  
  in = sf_input("in_file_ray");	    /* input data -> ray trajectory  */
  in2 = sf_input("in_file_ttime");  /* input data -> travel time */
  out = sf_output("out_file");      /* output radon */
  
  
  bool adj, inv, par;
  char *invmode;
  int iw, ip, ix, np, nt, nfft, nw, nx, niter;
  float dp, p0, dt, t0, dx, ox, x0, w, eps;
  float *p, *xx, **dd, **mm, *tmpr;
  sf_complex *cdd, *cmm;
  fftwf_complex *tmpc;
  fftwf_plan fft1, ifft1;
  sf_file offset=NULL;
  
    
  if (!sf_getbool("adj",&adj)) adj=true;
  /* if y, perform adjoint operation */
  if (!sf_getbool("inv",&inv)) inv=false; 
  /* if y, perform inverse operation */
  
  /* read input file parameters */
  nt = N;
  /* number of samples in time axis */
 
  
  if (adj||inv){ /* m(tau,p)=sum_{i=0}^{nx} d(t=tau+p*x_i,x_i) */
    
    nx = N;
    /* number of offset if the input in the data domain */
    
    /* specify slope axis */
    np = N;
    
    /* number of p values (if adj=y) */
    dp = 0.001;
    
    /* p sampling (if adj=y) */
    p0 = -0.2;
      
    /* p origin (if adj=y) */
    
    if(inv){			
      if ( !(invmode=sf_getstring("invmode")) ) invmode="toeplitz";
      /* inverse method: 'ls' if least-squares; 'toeplitz' -> 't' if use FFT */			
      if ((invmode[0]=='l') && (niter=100)); //!sf_getint("niter",&niter))
      /* number of CGLS iterations */
      if (eps=0.01) ;  //!sf_getfloat("eps",&eps)
      /* regularization parameter */
    } else {
      invmode=NULL;
    }
    
   
  }
  
    
  nfft=2*kiss_fft_next_fast_size(nt);
  nw=nfft/2+1;
  p=sf_floatalloc(np);
  xx=sf_floatalloc(nx);
  dd=sf_floatalloc2(nt, nx);
  mm=sf_floatalloc2(nt, np);
  cdd=(sf_complex*)malloc(nw*nx*sizeof(sf_complex));
  cmm=(sf_complex*)malloc(nw*np*sizeof(sf_complex));
  tmpr=(float*)fftwf_malloc(nfft*sizeof(float));
  tmpc=(fftwf_complex*)fftwf_malloc(nw*sizeof(fftwf_complex));
  fft1=fftwf_plan_dft_r2c_1d(nfft,tmpr,tmpc,FFTW_MEASURE);	
  ifft1=fftwf_plan_dft_c2r_1d(nfft,tmpc,tmpr,FFTW_MEASURE);
  
  for(ip=0; ip<np; ip++) p[ip]=p0+ip*dp;	
  if (adj||inv) {/* m(tau,p)=sum_{i=0}^{nx} d(t=tau+p*x_i,x_i) */
    sf_floatread(dd[0], nt*nx, in); 

    ox = 0.0;
    /* data origin in x */
    dx = 0.025;

    sf_floatread(mm[0], nt*np, in2);
  
    
  }
  
  if (NULL != offset) {
    sf_floatread(xx,nx,offset);
    sf_fileclose(offset);
  } else {
    for(ix=0; ix<nx; ix++) xx[ix]=ox+ix*dx;
  }
  
  if (!sf_getbool("parab",&par)) par=true; //par=false;
  /* if y, parabolic Radon transform */
  if (x0=1.) ;   
  /* reference offset */
	
	for (ix=0; ix < nx; ix++)/* normalize offsets */
	  {
	    if (par) xx[ix] *= xx[ix]/(x0*x0);
	    //else if (x0!=1.) xx[ix] /= x0;
	  }
	
	if(adj||inv){/* m(tau,p)=sum_{i=0}^{nx} d(t=tau+p*x_i,x_i) */
	  for(ix=0; ix<nx; ix++) /* loop over offsets */
	    {
	      memset(tmpr, 0, nfft*sizeof(float));
	      memcpy(tmpr, dd[ix], nt*sizeof(float));
	      fftwf_execute(fft1);/* FFT: dd-->cdd */
	      memcpy(&cdd[ix*nw], tmpc, nw*sizeof(sf_complex));
	    }
	  matrix_transpose(cdd, nw, nx);
	  
	  for(ip=0; ip<np; ip++) /* loop over slopes */
	    {
	      memset(tmpr, 0, nfft*sizeof(float));
	      memcpy(tmpr, mm[ip], nt*sizeof(float));
	      fftwf_execute(fft1);/* FFT: mm-->cmm */
	      memcpy(&cmm[ip*nw], tmpc, nw*sizeof(float));
	      			
	    }
	  matrix_transpose(cmm, nw, np);
	}
	

	myradon2_init(np, nx, dp, p, xx);
	for(iw=0; iw<nw; iw++) 
	  {
	    w=2.*SF_PI*iw/(nfft*dt);
	    myradon2_set(w);
	    myradon2_lop(adj, false, np, nx, &cmm[iw*np], &cdd[iw*nx]);//Linear operator RT
	    if(adj&&inv){
			if (invmode[0]=='t' )
			  myradon2_inv(&cmm[iw*np], &cmm[iw*np], eps);
			else
			  sf_csolver_reg(myradon2_lop, sf_ccgstep, sf_ccopy_lop, np, 
					 np, nx, &cmm[iw*np], &cdd[iw*nx], niter,eps,"end");
	    }
	  }
	
	
	if(adj||inv){// m(tau,p)=sum_{i=0}^{nx} d(t=tau+p*x_i,x_i) 
	  matrix_transpose(cmm, np, nw);
	  for(ip=0; ip<np; ip++) // loop over slopes // 
	    {			
	      memcpy(tmpc, &cmm[ip*nw], nw*sizeof(sf_complex));
	      fftwf_execute(ifft1); // IFFT: cmm-->mm //
	      for(iw=0; iw<nt; iw++) mm[ip][iw]=tmpr[iw]/nfft; //mm[ip][iw]=tmpr[iw]/nfft;
	    }
	  
	  sf_floatwrite(mm[0], nt*np, out);
	}

	
	
	free(p);
	free(xx);
	free(*dd); free(dd);
	free(*mm); free(mm);
	free(cdd); 
	free(cmm); 
	fftwf_free(tmpr);
	fftwf_free(tmpc);
	fftwf_destroy_plan(fft1);
    	fftwf_destroy_plan(ifft1);
	
    	exit(0);
}

/*Reads out_file  -the RT data- and computes the 1D FT-. To obtain the wave velocity a 2D IFT is applied to before step*/

int fourier(char *out_file)
{
  FILE *read = NULL;
  FILE *read2 = NULL;	
  FILE *out_wave = NULL;
  FILE *out_wave1 = NULL;
  FILE *out_wave2D = NULL;
  
  int i, j, nread;
 
  
  static sf_complex *dc;

  fftw_plan my_plan1, my_plan2, my_plan3;
  fftw_complex *in_radon, *in,*out_origin, *in_fourier,*out_fourier1D, *out_fourier2D, *FT_RT;

  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  in_radon = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  out_origin = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  out_fourier1D = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  out_fourier2D = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N);
  FT_RT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);

  //FT_RT = (double *) malloc(N* sizeof(double *));


  ///////////////////////////////////READING THE RADON DATA/////////////////////////
  
  read = fopen(out_file,"r");

  tini = clock();
      
  for(i=0; i<=N; i++)
    {
      nread = fscanf(read,"%lf",&RAY.path[i]);

      in_radon[i][0] = RAY.ttime[i];  //REAL PART

      in_radon[i][1] = 0.0;          //IMAGINARY PART 
      
    }

  ////////////////////////////////Calculating 1D FFT OF RT-DATA /////////////////////////
  
  my_plan1 = fftw_plan_dft_1d(N, in_radon, out_fourier1D, FFTW_FORWARD, FFTW_ESTIMATE);
   
  fftw_execute(my_plan1);

  out_wave = fopen("1DFT_of_RT.dat","w");

  for(i=0; i<=N; i++)
    {
      fprintf(out_wave,"%lf %lf\n", out_fourier1D[i][0], out_fourier1D[i][1]);
    }


 ////////////////////////////////Calculating 2D FFT OF THE ORIGINAL DATA ///////////////////

  read2 = fopen(in_file_ray,"r");
	
  for(i=0; i<=N; i++)
    {
      nread = fscanf(read2,"%lf",&RAY.path[i]);

      in[i][0] = RAY.path[i];  //REAL PART

      in[i][1] = 0.0;          //IMAGINARY PART 
      
    }
	
	
  my_plan2 = fftw_plan_dft_2d(N, N, in, out_origin, FFTW_FORWARD, FFTW_ESTIMATE);
   
  fftw_execute(my_plan2);

  out_wave1 = fopen("2DFT_ORIGINAL_DATA.dat","w");

  for(i=0; i<=N; i++)
    {
      
      fprintf(out_wave1,"%lf %lf\n", out_origin[i][0]/N*N, out_origin[i][1]/N*N);
    }


  
////////////////////////////////Calculating 2D IFFT OF 1DFT_RT-DATA /////////////////////////
//HERE WE OBTAIN THE SLOWNESS//

  
  for(i=0; i<=N; i++)
    {
      FT_RT[i][0] =  out_fourier1D[i][0];
      
    }

  
  
  my_plan3 = fftw_plan_dft_2d(N, N, FT_RT, out_fourier2D, FFTW_BACKWARD, FFTW_ESTIMATE);

   
  fftw_execute(my_plan3);

  out_wave2D = fopen("2DIFT_of_1DFT_RT.dat","w");

  for(i=0; i<=N; i++)
    {
      fprintf(out_wave2D,"%lf %lf\n", out_fourier2D[i][0]/N*N, out_fourier2D[i][1]/N*N);
      
    }


   
  tend = clock();
  
  cpu_time_used = ((double) (tend - tini)) / CLOCKS_PER_SEC;
  
  printf("CPU TIME USED: %16.8lf\n",cpu_time_used);
  
  
  fclose(read);
  fclose(read2);	
  fclose(out_wave);
  fclose(out_wave2D);


  fftw_destroy_plan(my_plan1);
  fftw_destroy_plan(my_plan2);
  fftw_destroy_plan(my_plan3);

  fftw_free(in_radon);
  fftw_free(out_wave);
  fftw_free(out_wave1);
  fftw_free(out_wave2D);

  return 0;
}


//////////////////////////////////MAIN PROGRAM/////////////////////
int main(int argc, char **argv){

 
  printf("%d\n",argc);

  if(argc != 5)
    {
      printf("ERROR--> use as:\n");
      printf("%s #iterations input_file_ray input_file_ttime  output_file\n",argv[0]);
      exit(0);  
    }
  
  N   = atoi(argv[1]);
  in_file_ray  = argv[2];
  in_file_ttime  = argv[3];
  out_file  = argv[2];

  printf("%s %d %s %s %s\n",argv[0], N, in_file_ray, in_file_ttime, out_file);

  RAY.path = (double *) malloc(N* sizeof(double *));
  RAY.ttime = (double *) malloc(N* sizeof(double *));

  /*Calculate the Radon Transform to rays */

  radon(in_file_ray, in_file_ttime, out_file);

  /*Calculate the 1D FT of RT and the 2D IFT of FT_RT -> COMPUTES CST: Central Slice Theorem*/
  
  fourier(out_file);
         
  return 0;
  
}










