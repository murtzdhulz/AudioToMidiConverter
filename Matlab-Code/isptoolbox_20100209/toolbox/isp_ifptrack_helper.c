#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void ifptrack(double *p, double *m, double *I, double *Sreal, double *Simag, int xlen, int ylen, double dgoodThrfminl, double fminl, double fminu, double fmaxl, double fmaxu);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, 
		 const mxArray *prhs[]) {

  /*
    [p,m] = ifptrack_helper(I, S, .75*sr/w, fminl,fminu,fmaxl,fmaxu);
  */

  double *I, *Sreal, *Simag, *p, *m;
  int xlen, ylen;
  double dgoodThr, fminl, fminu,  fmaxl,  fmaxu;
  if (nrhs != 7)
    mexErrMsgTxt("Wrong number of inputs.");

  if (nlhs != 2)
    mexErrMsgTxt("Exactly two outputs required.");

  I=mxGetPr(prhs[0]);
  xlen=mxGetM(prhs[0]);
  ylen=mxGetN(prhs[0]);
  Sreal=mxGetPr(prhs[1]);
  Simag=mxGetPi(prhs[1]);
  if (xlen != mxGetM(prhs[1]) || ylen != mxGetN(prhs[1]))
    mexErrMsgTxt("Sizes of I and S do not match.");
    
  dgoodThr=mxGetScalar(prhs[2]);
  fminl=mxGetScalar(prhs[3]);
  fminu=mxGetScalar(prhs[4]);
  fmaxl=mxGetScalar(prhs[5]);
  fmaxu=mxGetScalar(prhs[6]);

  plhs[0] = mxCreateDoubleMatrix(xlen, ylen, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(xlen, ylen, mxREAL);

  p=mxGetPr(plhs[0]);
  m=mxGetPr(plhs[1]);

  if (plhs[0]==NULL || plhs[1]==NULL)
    mexErrMsgTxt("Unable to allocate space for output matrices.");

  ifptrack(p, m, I, Sreal, Simag, xlen, ylen, dgoodThr, fminl, fminu, fmaxl, fmaxu);
  
}

void ifptrack(double *p, double *m, double *I, double *Sreal, double *Simag, int xlen, int ylen, double dgoodThr, double fminl, double fminu, double fmaxl, double fmaxu) {

  char *dgood = calloc(xlen*ylen, 1);
  double* ptrI=I;

  int y, x;
  char* ds=dgood;
  double bump, frq, mag;
  int t, bin;

  for (y=0; y<ylen; y++) {
    if (fabs(ptrI[1] - ptrI[0]) < dgoodThr)
      ds[0]=1;
    ds++;
    for (x=1; x<xlen-1; x++) {
      if (fabs(ptrI[x+1]-ptrI[x-1]) < dgoodThr)
	*ds = 1;
      ds++;
    }
    if (fabs(ptrI[xlen-1] - ptrI[xlen-2]) < dgoodThr)
      *ds=1;
    ds++;
    ptrI += xlen;
  }
  

  ds=dgood;
  for (y=0; y<ylen; y++) {
    ds++;
    for (x=xlen-2; x--;) {
      ds[0] &= (ds[-1] | ds[1]);
      ds++;
    }
    ds++;
  }

  ds=dgood;
  for (t=0; t<ylen; t++) {
    x=0;
   
    while (x<xlen) {
      while (!ds[x] && ++x<xlen);
      
      if (x==xlen)
	continue;
      
      bin=x;
      
      frq=0;
      mag=0;
      while (x<xlen && ds[x]) {
	bump = sqrt(Sreal[x]*Sreal[x] + Simag[x]*Simag[x]);
	mag += bump;
	frq += bump*I[x++];
      }
      if (mag)
	frq /= mag;


      if (frq<fminu) {
	if (frq<fminl) {
	  frq=0;
	  mag=0;
	} else {
	  mag = mag*(frq-fminl)/(fminu-fminl);
	}
      }

      if (frq>fmaxl){
	if (frq>fmaxu) {
	  frq=0;
	  mag=0;
	} else {
	  mag = mag*(fmaxu-frq)/(fmaxu-fmaxl);
	}
      }
      
      bin = floor(.5*(bin+x));
      
      p[bin] = frq;
      m[bin] = mag;
    }
    p += xlen;
    m += xlen;
    I += xlen;
    Sreal += xlen;
    Simag += xlen;
    ds += xlen;
  }

  free(dgood);
  
}
