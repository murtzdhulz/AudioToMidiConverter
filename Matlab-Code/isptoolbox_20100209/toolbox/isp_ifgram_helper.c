#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void ifgram(double *X, int N, int W, int H, double SR, int maxbin, int nhops, double *Freal, double *Fimag, double *D);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  double *X, SR, *Dreal, *Dimag, *F;
  int N, xlen, W, H, maxbin, nhops;

  if (nrhs != 6)
    mexErrMsgTxt("Exactly six input arguments expected.");

  if (nlhs != 2)
    mexErrMsgTxt("Exactly two output arguments expected.");

  X = mxGetPr(prhs[0]);
  xlen = mxGetM(prhs[0])*mxGetN(prhs[0]);

  /* Decided to choose defaults in Matlab
     N = (nlhs>1) ? mxGetScalar(prhs[1]):256;
     W = (nlhs>2) ? mxGetScalar(prhs[2]):N;
     H = (nlhs>3) ? mxGetScalar(prhs[3]):W/2;
     SR = (nlhs>4) ? mxGetScalar(prhs[4]):SR;
     maxbin = (nlhs>5) ? mxGetScalar(prhs[5]):(1+N/2);
  */

  N = mxGetScalar(prhs[1]);
  W = mxGetScalar(prhs[2]);
  H = mxGetScalar(prhs[3]);
  SR = mxGetScalar(prhs[4]);
  maxbin = mxGetScalar(prhs[5]);

  nhops = 1 + floor((xlen-W)/H);

  plhs[0] = mxCreateDoubleMatrix(maxbin, nhops, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(maxbin, nhops, mxCOMPLEX);

  if ((plhs[0]==NULL) || (plhs[1]==NULL))
    mexErrMsgTxt("Unable to allocate space for output.");

  F = mxGetPr(plhs[0]);
  Dreal = mxGetPr(plhs[1]);
  Dimag = mxGetPi(plhs[1]);

  ifgram(X, N, W, H, SR, maxbin, nhops, F, Dreal, Dimag);
  
}

void ifgram(double *X, int N, int W, int H, double SR, int maxbin, int nhops, double *F, double *Dreal, double *Dimag) {
  double *win, *dwin, norm, *ww;
  int i, tmp, nmw1;

  mxArray *mxwu, *mxdu, *mxwufft, *mxdufft;
  double *du, *wu;
  double *t1real, *t1imag, *t2real, *t2imag;
  double T;  
  int split;
  int h;
  int errdu, errwu;
  
  win = malloc(W*sizeof(win[0]));
  dwin = malloc(W*sizeof(dwin[0]));
  ww = malloc(maxbin*sizeof(ww[0]));

  T = ((double) W)/SR;
  norm = 0;
  for (i=0; i<W; i++) {
    win[i] = 0.5*(1-cos(2.0*M_PI*((double) i)/((double) W)));
    dwin[i] = - M_PI/T * sin((double) i/((double) W) * 2.0*M_PI);
    norm += win[i];
  }
  norm = 2/norm;

  for (i=0; i<maxbin; i++) {
    ww[i]=2.0*M_PI*i*SR/N;
  }

  nmw1 = floor( .5*(N-W) );
  split=ceil(.5*N);

  mxdu = mxCreateDoubleMatrix(N, 1, mxREAL);
  mxwu = mxCreateDoubleMatrix(N, 1, mxREAL);

  du = mxGetPr(mxdu);
  wu = mxGetPr(mxwu);

  for (h=0; h<nhops; h++) {
    double a,b, da, db, sign;

    /* Pad or truncate samples if N != W, and apply windows now, while
       the length is right */
    if (N>=W) {
      for (i=0; i<W; i++) {
	wu[(nmw1+split+i)%N] = X[i]*win[i];
	du[(nmw1+split+i)%N] = X[i]*dwin[i];
      }
    } else {
      for (i=-nmw1; i<N-nmw1; i++) {
	wu[(N+i+nmw1-split)%N] = X[i]*win[i];
	du[(N+i+nmw1-split)%N] = X[i]*dwin[i];
      }
    }

    errdu = mexCallMATLAB(1, &mxdufft, 1, &mxdu, "fft");
    errwu = mexCallMATLAB(1, &mxwufft, 1, &mxwu, "fft");

    if (errdu || errwu) {
      mexErrMsgTxt("Error calling fft command.");
    }

    t1real = mxGetPr(mxdufft);
    t1imag = mxGetPi(mxdufft);
    t2real = mxGetPr(mxwufft);
    t2imag = mxGetPi(mxwufft);

    if (t1real==NULL || t2real==NULL)
      mexErrMsgTxt("Error calling fft command.");
      
    for (i=0; i<maxbin; i++) {
      a = t2real[i];
      b = (t2imag==NULL) ? 0:t2imag[i];
      da = t1real[i] - ww[i]*b;
      db = ((t1imag==NULL) ? 0:t1imag[i]) + ww[i]*a;


      if (a || b) {
	*(F++) = (a*db-b*da)/((a*a+b*b)*(2*M_PI));
	*(Dreal++) = a*norm;
	*(Dimag++) = -b*norm;
      } else {
	*(F++) = 0;
	*(Dreal++) = 0;
	*(Dimag++) = 0;
      }

    }
    
    mxDestroyArray(mxdufft);
    mxDestroyArray(mxwufft);

    X += H;
  }

  free(win);
  free(dwin);
  free(ww);
}


