/* This file is automatically generated. DO NOT EDIT! */

#ifndef _myradon2_h
#define _myradon2_h


void myradon2_init(int np_, int nx_, double dp_, double *p_, double *xx_);
/*< initialization >*/


void myradon2_set(double w_);
/*< set up frequency w >*/


void myradon2_lop(bool adj, bool add, int nm, int nd, complex *mm, complex *dd);
/*< radon linear operator >*/


void myradon2_inv(complex *mm, complex *adj_dd, double eps);
/*< fast Toeplitz matrix inversion for radon transform 
  mm: model to be inverted
  adj_dd: adjoint radon of data
  eps: regularization parameter
  >*/

#endif
