/* Linear/parabolic radon operator in frequency domain
   Note: I borrowed a lot from /system/seismic/radon.c+Mradon.c. 
   The distinction:
   1) I am using FFTW because I am inexperienced in invoking kiss_fft. 
   2) I am using FFT-based Toeplitz inversion, while radon.c uses Levinson 
   recursion.
*/
/*
  Copyright (C) 2014  Xi'an Jiaotong University, UT Austin (Pengliang Yang)

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

  References: 
  1) Kostov C., 1990. "Toeplitz structure in Slant-Stack inversion": 
  SEG Extended Abstracts, 1647-1650.
  2) Sacchi, Mauricio D., and Milton Porsani. "Fast high resolution 
  parabolic Radon transform." Society of Exploration Geophysicists 
  69th Annual International Meeting, SPRO P. Vol. 1. No. 1. 1999.
  3) Vogel, Curtis R. Computational methods for inverse problems. 
  Vol. 23. Siam, 2002.	[Chapter 5.2.]
*/

#include <stdlib.h>
#include <time.h>
#include <complex.h>
#include <stdbool.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "myradon2.h"
#include "ctoeplitz_reg.c"

static int np, nx;
static double w, dp, *p, *xx;

void myradon2_init(int np_, int nx_, double dp_, double *p_, double *xx_)
/*< initialization >*/
{
    np=np_;
    nx=nx_;
    dp=dp_;
    p=p_;
    xx=xx_;
}

void myradon2_set(double w_)
/*< set up frequency w >*/
{
    w=w_;
}

void myradon2_lop(bool adj, bool add, int nm, int nd, complex *mm, complex *dd)
/*< radon linear operator >*/
{
    int ix, ip;
    complex sumc;

    if (nm != np || nd != nx)
      printf("mismatched data sizes\n");
	
    
    if(adj){/* mm(p,w)=sum_{ix=0}^{nx} dd(xx[ix],w)*exp(i*w*p*xx[ix]) */
	for(ip=0; ip<np; ip++) /* loop over slopes */
	{
	  sumc= 0.0 + 0.0 *I;             
	    for(ix=0; ix<nx; ix++) {

	      sumc+=expf(w*p[ip]*xx[ix])*dd[ix];
	     

	    }
	    mm[ip]=sumc;
	}
    }

    
}


static bool allocated=false;
static complex *c;

void myradon2_inv(complex *mm, complex *adj_dd, double eps)
/*< fast Toeplitz matrix inversion for radon transform 
  mm: model to be inverted
  adj_dd: adjoint radon of data
  eps: regularization parameter
  >*/
{
    int ip, ix;
    complex sumc;
    if (!allocated){
      c = (complex*)malloc(np*sizeof(complex));
      allocated=true;
    }
	
    for(ip=0; ip<np; ip++) 
    {
      sumc= 0.0 + 0.0 *I;
      for(ix=0; ix<nx; ix++) {

	    sumc+=expf(w*ip*dp*xx[ix]);

	}
	c[ip]=sumc;
    }
	
    ctoeplitz_inv(np, eps, c, mm, adj_dd);
}
