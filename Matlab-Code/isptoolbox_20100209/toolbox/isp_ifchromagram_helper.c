#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void hz2octs(double *p, int len);
void quantize(double *p, double *p_org, int len, int nchr);
void chromagram_IF(double *oct, double *p, double *m, int xlen, int ylen, int nchr);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, 
		 const mxArray *prhs[]) {

  /*
    Y = chromagram_IF_helper(p, m, nchr);
  */

  double *oct, *p, *m;
  int xlen, ylen, nchr;

  if (nrhs != 3)
    mexErrMsgTxt("Wrong number of inputs.");

  if (nlhs != 1)
    mexErrMsgTxt("Exactly one output required.");

  p=mxGetPr(prhs[0]);
  xlen=mxGetM(prhs[0]);
  ylen=mxGetN(prhs[0]);
  m=mxGetPr(prhs[1]);
  if (xlen != mxGetM(prhs[1]) || ylen != mxGetN(prhs[1]))
    mexErrMsgTxt("Sizes of argument one and two do not match.");
    
  nchr=mxGetScalar(prhs[2]);

  plhs[0] = mxCreateDoubleMatrix(nchr, ylen, mxREAL);
  if (plhs[0]==NULL)
    mexErrMsgTxt("Unable to allocate space for output matrix.");

  if (xlen*ylen == 0)
    return;

  oct = mxGetPr(plhs[0]);
  chromagram_IF(oct, p, m, xlen, ylen, nchr);
}

void chromagram_IF(double *oct, double *p_org, double *m, int xlen, int ylen, int nchr) {

  int x,y;

  /* The next four lines of code were added 'cause it isn't nice to
     alter the input arguments */
  double *p;
  double *pp;
  double *mp;
  double *octp;

  p = malloc(xlen*ylen*sizeof(*p));
  for (x=0; x<xlen*ylen; x++)
    p[x] = p_org[x];

  pp = p;
  mp = m;
  octp = oct;
  
  hz2octs(p, xlen*ylen);
  quantize(p, p_org, xlen*ylen, nchr);

  for (y=0; y<ylen; y++) {
    for (x=0; x<xlen; x++) {
      if (*(p_org++)) {
	octp[((int) pp[x]) % nchr] += mp[x];
      }
    }
    mp += xlen;
    pp += xlen;
    octp += nchr;
  }
  free(p);
}


void quantize(double *p, double *p_org, int len, int nchr) {
  int hist[100];
  int i;
  int best;
  double correction;
  
  for (i=0; i<100; i++) 
    hist[i]=0;

  for (i=0; i<len; i++) {
    if (p_org[i]) {
      p[i] *= nchr;
      hist[(int) floor(100*(0.5+p[i]-floor(0.5+p[i])))]++;
    }
  }

  best = 0;
  for (i=1; i<100; i++) {
    if (hist[i]>hist[best])
      best=i;
  }

  correction=0.5-((.5+(double) best)/100 - 0.5);
  
  for (i=len-1; i--;)
    if (p[i]) 
      p[i] = floor(p[i]+correction);
}

void hz2octs(double *p, int len) {
  int i;
  double c = 1.0/(440.0/16.0);
  for (i=len-1; i--;) {
    if (p[i])
      p[i]=log(c*p[i])/M_LN2;
  }

}
