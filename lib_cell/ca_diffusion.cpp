#include "subcell.hpp"

void CSubcell::computeIci(void)
{

#ifdef ___PERIODIC
  Ici[0+0*nnx+0*nnxnny]=(ci[0+0*nnx+0*nnxnny+1]+ci[(nnx-1)+0*nnx+0*nnxnny]-2*ci[0+0*nnx+0*nnxnny])/tauiL+
    (ci[0+(0+1)*nnx+0*nnxnny]+ci[0+(nny-1)*nnx+0*nnxnny]-2*ci[0+0*nnx+0*nnxnny])/tauiT+
    (ci[0+0*nnx+(0+1)*nnxnny]+ci[0+0*nnx+(nnz-1)*nnxnny]-2*ci[0+0*nnx+0*nnxnny])/tauiT;
  
  Ici[0+(nny-1)*nnx+(nnz-1)*nnxnny]=(ci[0+(nny-1)*nnx+(nnz-1)*nnxnny+1]+ci[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]-2*ci[0+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiL+
    (ci[0+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+ci[0+(0)*nnx+(nnz-1)*nnxnny]-2*ci[0+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiT+
    (ci[0+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+ci[0+(nny-1)*nnx+(0)*nnxnny]-2*ci[0+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiT;
  
  Ici[0+(nny-1)*nnx+0*nnxnny]=(ci[0+(nny-1)*nnx+0*nnxnny+1]+ci[(nnx-1)+(nny-1)*nnx+0*nnxnny]-2*ci[0+(nny-1)*nnx+0*nnxnny])/tauiL+
    (ci[0+((nny-1)-1)*nnx+0*nnxnny]+ci[0+(0)*nnx+0*nnxnny]-2*ci[0+(nny-1)*nnx+0*nnxnny])/tauiT+
    (ci[0+(nny-1)*nnx+(0+1)*nnxnny]+ci[0+(nny-1)*nnx+(nnz-1)*nnxnny]-2*ci[0+(nny-1)*nnx+0*nnxnny])/tauiT;
  
  Ici[0+0*nnx+(nnz-1)*nnxnny]=(ci[0+0*nnx+(nnz-1)*nnxnny+1]+ci[(nnx-1)+0*nnx+(nnz-1)*nnxnny]-2*ci[0+0*nnx+(nnz-1)*nnxnny])/tauiL+
    (ci[0+(0+1)*nnx+(nnz-1)*nnxnny]+ci[0+(nny-1)*nnx+(nnz-1)*nnxnny]-2*ci[0+0*nnx+(nnz-1)*nnxnny])/tauiT+
    (ci[0+0*nnx+((nnz-1)-1)*nnxnny]+ci[0+0*nnx+(0)*nnxnny]-2*ci[0+0*nnx+(nnz-1)*nnxnny])/tauiT;
  
  Ici[(nnx-1)+(nny-1)*nnx+0*nnxnny]=(ci[(nnx-1)+(nny-1)*nnx+0*nnxnny-1]+ci[(0)+(nny-1)*nnx+0*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+0*nnxnny])/tauiL+
    (ci[(nnx-1)+((nny-1)-1)*nnx+0*nnxnny]+ci[(nnx-1)+(0)*nnx+0*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+0*nnxnny])/tauiT+
    (ci[(nnx-1)+(nny-1)*nnx+(0+1)*nnxnny]+ci[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+0*nnxnny])/tauiT;
  
  Ici[(nnx-1)+0*nnx+(nnz-1)*nnxnny]=(ci[(nnx-1)+0*nnx+(nnz-1)*nnxnny-1]+ci[(0)+0*nnx+(nnz-1)*nnxnny]-2*ci[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/tauiL+
    (ci[(nnx-1)+(0+1)*nnx+(nnz-1)*nnxnny]+ci[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]-2*ci[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/tauiT+
    (ci[(nnx-1)+0*nnx+((nnz-1)-1)*nnxnny]+ci[(nnx-1)+0*nnx+(0)*nnxnny]-2*ci[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/tauiT;
  
  Ici[(nnx-1)+0*nnx+0*nnxnny]=(ci[(nnx-1)+0*nnx+0*nnxnny-1]+ci[(0)+0*nnx+0*nnxnny]-2*ci[(nnx-1)+0*nnx+0*nnxnny])/tauiL+
    (ci[(nnx-1)+(0+1)*nnx+0*nnxnny]+ci[(nnx-1)+(nny-1)*nnx+0*nnxnny]-2*ci[(nnx-1)+0*nnx+0*nnxnny])/tauiT+
    (ci[(nnx-1)+0*nnx+(0+1)*nnxnny]+ci[(nnx-1)+0*nnx+(nnz-1)*nnxnny]-2*ci[(nnx-1)+0*nnx+0*nnxnny])/tauiT;
  
  Ici[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]=(ci[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny-1]+ci[(0)+(nny-1)*nnx+(nnz-1)*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiL+
    (ci[(nnx-1)+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+ci[(nnx-1)+(0)*nnx+(nnz-1)*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiT+
    (ci[(nnx-1)+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+ci[(nnx-1)+(nny-1)*nnx+(0)*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiT;
  //x fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int j=1;j<nny-1;j++)
  {
    Ici[0+j*nnx+0*nnxnny]=(ci[0+j*nnx+0*nnxnny+1]+ci[(nnx-1)+j*nnx+0*nnxnny]-2*ci[0+j*nnx+0*nnxnny])/tauiL+
      (ci[0+(j+1)*nnx+0*nnxnny]+ci[0+(j-1)*nnx+0*nnxnny]-2*ci[0+j*nnx+0*nnxnny])/tauiT+
      (ci[0+j*nnx+(0+1)*nnxnny]+ci[0+j*nnx+(nnz-1)*nnxnny]-2*ci[0+j*nnx+0*nnxnny])/tauiT;
    Ici[(nnx-1)+j*nnx+0*nnxnny]=(ci[(nnx-1)+j*nnx+0*nnxnny-1]+ci[(0)+j*nnx+0*nnxnny]-2*ci[(nnx-1)+j*nnx+0*nnxnny])/tauiL+
      (ci[(nnx-1)+(j+1)*nnx+0*nnxnny]+ci[(nnx-1)+(j-1)*nnx+0*nnxnny]-2*ci[(nnx-1)+j*nnx+0*nnxnny])/tauiT+
      (ci[(nnx-1)+j*nnx+(0+1)*nnxnny]+ci[(nnx-1)+j*nnx+(nnz-1)*nnxnny]-2*ci[(nnx-1)+j*nnx+0*nnxnny])/tauiT;
    Ici[0+j*nnx+(nnz-1)*nnxnny]=(ci[0+j*nnx+(nnz-1)*nnxnny+1]+ci[(nnx-1)+j*nnx+(nnz-1)*nnxnny]-2*ci[0+j*nnx+(nnz-1)*nnxnny])/tauiL+
      (ci[0+(j+1)*nnx+(nnz-1)*nnxnny]+ci[0+(j-1)*nnx+(nnz-1)*nnxnny]-2*ci[0+j*nnx+(nnz-1)*nnxnny])/tauiT+
      (ci[0+j*nnx+((nnz-1)-1)*nnxnny]+ci[0+j*nnx+(0)*nnxnny]-2*ci[0+j*nnx+(nnz-1)*nnxnny])/tauiT;
    Ici[(nnx-1)+j*nnx+(nnz-1)*nnxnny]=(ci[(nnx-1)+j*nnx+(nnz-1)*nnxnny-1]+ci[(0)+j*nnx+(nnz-1)*nnxnny]-2*ci[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/tauiL+
      (ci[(nnx-1)+(j+1)*nnx+(nnz-1)*nnxnny]+ci[(nnx-1)+(j-1)*nnx+(nnz-1)*nnxnny]-2*ci[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/tauiT+
      (ci[(nnx-1)+j*nnx+((nnz-1)-1)*nnxnny]+ci[(nnx-1)+j*nnx+(0)*nnxnny]-2*ci[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/tauiT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      Ici[0+j*nnx+k*nnxnny]=(ci[0+j*nnx+k*nnxnny+1]+ci[(nnx-1)+j*nnx+k*nnxnny]-2*ci[0+j*nnx+k*nnxnny])/tauiL+
        (ci[0+(j+1)*nnx+k*nnxnny]+ci[0+(j-1)*nnx+k*nnxnny]-2*ci[0+j*nnx+k*nnxnny])/tauiT+
        (ci[0+j*nnx+(k+1)*nnxnny]+ci[0+j*nnx+(k-1)*nnxnny]-2*ci[0+j*nnx+k*nnxnny])/tauiT;
      Ici[(nnx-1)+j*nnx+k*nnxnny]=(ci[(nnx-1)+j*nnx+k*nnxnny-1]+ci[(0)+j*nnx+k*nnxnny]-2*ci[(nnx-1)+j*nnx+k*nnxnny])/tauiL+
        (ci[(nnx-1)+(j+1)*nnx+k*nnxnny]+ci[(nnx-1)+(j-1)*nnx+k*nnxnny]-2*ci[(nnx-1)+j*nnx+k*nnxnny])/tauiT+
        (ci[(nnx-1)+j*nnx+(k+1)*nnxnny]+ci[(nnx-1)+j*nnx+(k-1)*nnxnny]-2*ci[(nnx-1)+j*nnx+k*nnxnny])/tauiT;
    }
  }
  //y fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=1;i<(nnx-1);i++)
  {
    Ici[i+0*nnx+0*nnxnny]=(ci[i+0*nnx+0*nnxnny+1]+ci[i+0*nnx+0*nnxnny-1]-2*ci[i+0*nnx+0*nnxnny])/tauiL+
      (ci[i+(0+1)*nnx+0*nnxnny]+ci[i+(nny-1)*nnx+0*nnxnny]-2*ci[i+0*nnx+0*nnxnny])/tauiT+
      (ci[i+0*nnx+(0+1)*nnxnny]+ci[i+0*nnx+(nnz-1)*nnxnny]-2*ci[i+0*nnx+0*nnxnny])/tauiT;
    Ici[i+(nny-1)*nnx+0*nnxnny]=(ci[i+(nny-1)*nnx+0*nnxnny+1]+ci[i+(nny-1)*nnx+0*nnxnny-1]-2*ci[i+(nny-1)*nnx+0*nnxnny])/tauiL+
      (ci[i+((nny-1)-1)*nnx+0*nnxnny]+ci[i+(0)*nnx+0*nnxnny]-2*ci[i+(nny-1)*nnx+0*nnxnny])/tauiT+
      (ci[i+(nny-1)*nnx+(0+1)*nnxnny]+ci[i+(nny-1)*nnx+(nnz-1)*nnxnny]-2*ci[i+(nny-1)*nnx+0*nnxnny])/tauiT;
    Ici[i+0*nnx+(nnz-1)*nnxnny]=(ci[i+0*nnx+(nnz-1)*nnxnny+1]+ci[i+0*nnx+(nnz-1)*nnxnny-1]-2*ci[i+0*nnx+(nnz-1)*nnxnny])/tauiL+
      (ci[i+(0+1)*nnx+(nnz-1)*nnxnny]+ci[i+(nny-1)*nnx+(nnz-1)*nnxnny]-2*ci[i+0*nnx+(nnz-1)*nnxnny])/tauiT+
      (ci[i+0*nnx+((nnz-1)-1)*nnxnny]+ci[i+0*nnx+(0)*nnxnny]-2*ci[i+0*nnx+(nnz-1)*nnxnny])/tauiT;
    Ici[i+(nny-1)*nnx+(nnz-1)*nnxnny]=(ci[i+(nny-1)*nnx+(nnz-1)*nnxnny+1]+ci[i+(nny-1)*nnx+(nnz-1)*nnxnny-1]-2*ci[i+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiL+
      (ci[i+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+ci[i+(0)*nnx+(nnz-1)*nnxnny]-2*ci[i+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiT+
      (ci[i+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+ci[i+(nny-1)*nnx+(0)*nnxnny]-2*ci[i+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      Ici[i+0*nnx+k*nnxnny]=(ci[i+0*nnx+k*nnxnny+1]+ci[i+0*nnx+k*nnxnny-1]-2*ci[i+0*nnx+k*nnxnny])/tauiL+
        (ci[i+(0+1)*nnx+k*nnxnny]+ci[i+(nny-1)*nnx+k*nnxnny]-2*ci[i+0*nnx+k*nnxnny])/tauiT+
        (ci[i+0*nnx+(k+1)*nnxnny]+ci[i+0*nnx+(k-1)*nnxnny]-2*ci[i+0*nnx+k*nnxnny])/tauiT;
      Ici[i+(nny-1)*nnx+k*nnxnny]=(ci[i+(nny-1)*nnx+k*nnxnny+1]+ci[i+(nny-1)*nnx+k*nnxnny-1]-2*ci[i+(nny-1)*nnx+k*nnxnny])/tauiL+
        (ci[i+((nny-1)-1)*nnx+k*nnxnny]+ci[i+(0)*nnx+k*nnxnny]-2*ci[i+(nny-1)*nnx+k*nnxnny])/tauiT+
        (ci[i+(nny-1)*nnx+(k+1)*nnxnny]+ci[i+(nny-1)*nnx+(k-1)*nnxnny]-2*ci[i+(nny-1)*nnx+k*nnxnny])/tauiT;
    }
    //z fixed
#pragma ivdep
#pragma vector always
    for (int j=1;j<nny-1;j++)
    {
      Ici[i+j*nnx+0*nnxnny]=(ci[i+j*nnx+0*nnxnny+1]+ci[i+j*nnx+0*nnxnny-1]-2*ci[i+j*nnx+0*nnxnny])/tauiL+
        (ci[i+(j+1)*nnx+0*nnxnny]+ci[i+(j-1)*nnx+0*nnxnny]-2*ci[i+j*nnx+0*nnxnny])/tauiT+
        (ci[i+j*nnx+(0+1)*nnxnny]+ci[i+j*nnx+(nnz-1)*nnxnny]-2*ci[i+j*nnx+0*nnxnny])/tauiT;
      Ici[i+j*nnx+(nnz-1)*nnxnny]=(ci[i+j*nnx+(nnz-1)*nnxnny+1]+ci[i+j*nnx+(nnz-1)*nnxnny-1]-2*ci[i+j*nnx+(nnz-1)*nnxnny])/tauiL+
        (ci[i+(j+1)*nnx+(nnz-1)*nnxnny]+ci[i+(j-1)*nnx+(nnz-1)*nnxnny]-2*ci[i+j*nnx+(nnz-1)*nnxnny])/tauiT+
        (ci[i+j*nnx+((nnz-1)-1)*nnxnny]+ci[i+j*nnx+(0)*nnxnny]-2*ci[i+j*nnx+(nnz-1)*nnxnny])/tauiT;
    }
  }
#pragma ivdep
#pragma vector always
  for (int k=1;k<nnz-1;k++)
  {
    Ici[0+0*nnx+k*nnxnny]=(ci[0+0*nnx+k*nnxnny+1]+ci[(nnx-1)+0*nnx+k*nnxnny]-2*ci[0+0*nnx+k*nnxnny])/tauiL+
      (ci[0+(0+1)*nnx+k*nnxnny]+ci[0+(nny-1)*nnx+k*nnxnny]-2*ci[0+0*nnx+k*nnxnny])/tauiT+
      (ci[0+0*nnx+(k+1)*nnxnny]+ci[0+0*nnx+(k-1)*nnxnny]-2*ci[0+0*nnx+k*nnxnny])/tauiT;
    Ici[0+(nny-1)*nnx+k*nnxnny]=(ci[0+(nny-1)*nnx+k*nnxnny+1]+ci[(nnx-1)+(nny-1)*nnx+k*nnxnny]-2*ci[0+(nny-1)*nnx+k*nnxnny])/tauiL+
      (ci[0+((nny-1)-1)*nnx+k*nnxnny]+ci[0+(0)*nnx+k*nnxnny]-2*ci[0+(nny-1)*nnx+k*nnxnny])/tauiT+
      (ci[0+(nny-1)*nnx+(k+1)*nnxnny]+ci[0+(nny-1)*nnx+(k-1)*nnxnny]-2*ci[0+(nny-1)*nnx+k*nnxnny])/tauiT;
    Ici[(nnx-1)+0*nnx+k*nnxnny]=(ci[(nnx-1)+0*nnx+k*nnxnny-1]+ci[(0)+0*nnx+k*nnxnny]-2*ci[(nnx-1)+0*nnx+k*nnxnny])/tauiL+
      (ci[(nnx-1)+(0+1)*nnx+k*nnxnny]+ci[(nnx-1)+(nny-1)*nnx+k*nnxnny]-2*ci[(nnx-1)+0*nnx+k*nnxnny])/tauiT+
      (ci[(nnx-1)+0*nnx+(k+1)*nnxnny]+ci[(nnx-1)+0*nnx+(k-1)*nnxnny]-2*ci[(nnx-1)+0*nnx+k*nnxnny])/tauiT;
    Ici[(nnx-1)+(nny-1)*nnx+k*nnxnny]=(ci[(nnx-1)+(nny-1)*nnx+k*nnxnny-1]+ci[(0)+(nny-1)*nnx+k*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+k*nnxnny])/tauiL+
      (ci[(nnx-1)+((nny-1)-1)*nnx+k*nnxnny]+ci[(nnx-1)+(0)*nnx+k*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+k*nnxnny])/tauiT+
      (ci[(nnx-1)+(nny-1)*nnx+(k+1)*nnxnny]+ci[(nnx-1)+(nny-1)*nnx+(k-1)*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+k*nnxnny])/tauiT;
  }
#else
  Ici[0+0*nnx+0*nnxnny]=(ci[0+0*nnx+0*nnxnny+1]+ci[0+0*nnx+0*nnxnny+1]-2*ci[0+0*nnx+0*nnxnny])/tauiL+
    (ci[0+(0+1)*nnx+0*nnxnny]+ci[0+(0+1)*nnx+0*nnxnny]-2*ci[0+0*nnx+0*nnxnny])/tauiT+
    (ci[0+0*nnx+(0+1)*nnxnny]+ci[0+0*nnx+(0+1)*nnxnny]-2*ci[0+0*nnx+0*nnxnny])/tauiT;
  Ici[0+(nny-1)*nnx+(nnz-1)*nnxnny]=(ci[0+(nny-1)*nnx+(nnz-1)*nnxnny+1]+ci[0+(nny-1)*nnx+(nnz-1)*nnxnny+1]-2*ci[0+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiL+
    (ci[0+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+ci[0+((nny-1)-1)*nnx+(nnz-1)*nnxnny]-2*ci[0+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiT+
    (ci[0+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+ci[0+(nny-1)*nnx+((nnz-1)-1)*nnxnny]-2*ci[0+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiT;
  Ici[0+(nny-1)*nnx+0*nnxnny]=(ci[0+(nny-1)*nnx+0*nnxnny+1]+ci[0+(nny-1)*nnx+0*nnxnny+1]-2*ci[0+(nny-1)*nnx+0*nnxnny])/tauiL+
    (ci[0+((nny-1)-1)*nnx+0*nnxnny]+ci[0+((nny-1)-1)*nnx+0*nnxnny]-2*ci[0+(nny-1)*nnx+0*nnxnny])/tauiT+
    (ci[0+(nny-1)*nnx+(0+1)*nnxnny]+ci[0+(nny-1)*nnx+(0+1)*nnxnny]-2*ci[0+(nny-1)*nnx+0*nnxnny])/tauiT;
  Ici[0+0*nnx+(nnz-1)*nnxnny]=(ci[0+0*nnx+(nnz-1)*nnxnny+1]+ci[0+0*nnx+(nnz-1)*nnxnny+1]-2*ci[0+0*nnx+(nnz-1)*nnxnny])/tauiL+
    (ci[0+(0+1)*nnx+(nnz-1)*nnxnny]+ci[0+(0+1)*nnx+(nnz-1)*nnxnny]-2*ci[0+0*nnx+(nnz-1)*nnxnny])/tauiT+
    (ci[0+0*nnx+((nnz-1)-1)*nnxnny]+ci[0+0*nnx+((nnz-1)-1)*nnxnny]-2*ci[0+0*nnx+(nnz-1)*nnxnny])/tauiT;
  Ici[(nnx-1)+(nny-1)*nnx+0*nnxnny]=(ci[(nnx-1)+(nny-1)*nnx+0*nnxnny-1]+ci[(nnx-1)+(nny-1)*nnx+0*nnxnny-1]-2*ci[(nnx-1)+(nny-1)*nnx+0*nnxnny])/tauiL+
    (ci[(nnx-1)+((nny-1)-1)*nnx+0*nnxnny]+ci[(nnx-1)+((nny-1)-1)*nnx+0*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+0*nnxnny])/tauiT+
    (ci[(nnx-1)+(nny-1)*nnx+(0+1)*nnxnny]+ci[(nnx-1)+(nny-1)*nnx+(0+1)*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+0*nnxnny])/tauiT;
  Ici[(nnx-1)+0*nnx+(nnz-1)*nnxnny]=(ci[(nnx-1)+0*nnx+(nnz-1)*nnxnny-1]+ci[(nnx-1)+0*nnx+(nnz-1)*nnxnny-1]-2*ci[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/tauiL+
    (ci[(nnx-1)+(0+1)*nnx+(nnz-1)*nnxnny]+ci[(nnx-1)+(0+1)*nnx+(nnz-1)*nnxnny]-2*ci[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/tauiT+
    (ci[(nnx-1)+0*nnx+((nnz-1)-1)*nnxnny]+ci[(nnx-1)+0*nnx+((nnz-1)-1)*nnxnny]-2*ci[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/tauiT;
  Ici[(nnx-1)+0*nnx+0*nnxnny]=(ci[(nnx-1)+0*nnx+0*nnxnny-1]+ci[(nnx-1)+0*nnx+0*nnxnny-1]-2*ci[(nnx-1)+0*nnx+0*nnxnny])/tauiL+
    (ci[(nnx-1)+(0+1)*nnx+0*nnxnny]+ci[(nnx-1)+(0+1)*nnx+0*nnxnny]-2*ci[(nnx-1)+0*nnx+0*nnxnny])/tauiT+
    (ci[(nnx-1)+0*nnx+(0+1)*nnxnny]+ci[(nnx-1)+0*nnx+(0+1)*nnxnny]-2*ci[(nnx-1)+0*nnx+0*nnxnny])/tauiT;
  Ici[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]=(ci[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny-1]+ci[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny-1]-2*ci[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiL+
    (ci[(nnx-1)+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+ci[(nnx-1)+((nny-1)-1)*nnx+(nnz-1)*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiT+
    (ci[(nnx-1)+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+ci[(nnx-1)+(nny-1)*nnx+((nnz-1)-1)*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiT;
  //x fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int j=1;j<nny-1;j++)
  {
    Ici[0+j*nnx+0*nnxnny]=(ci[0+j*nnx+0*nnxnny+1]+ci[0+j*nnx+0*nnxnny+1]-2*ci[0+j*nnx+0*nnxnny])/tauiL+
      (ci[0+(j+1)*nnx+0*nnxnny]+ci[0+(j-1)*nnx+0*nnxnny]-2*ci[0+j*nnx+0*nnxnny])/tauiT+
      (ci[0+j*nnx+(0+1)*nnxnny]+ci[0+j*nnx+(0+1)*nnxnny]-2*ci[0+j*nnx+0*nnxnny])/tauiT;
    Ici[(nnx-1)+j*nnx+0*nnxnny]=(ci[(nnx-1)+j*nnx+0*nnxnny-1]+ci[(nnx-1)+j*nnx+0*nnxnny-1]-2*ci[(nnx-1)+j*nnx+0*nnxnny])/tauiL+
      (ci[(nnx-1)+(j+1)*nnx+0*nnxnny]+ci[(nnx-1)+(j-1)*nnx+0*nnxnny]-2*ci[(nnx-1)+j*nnx+0*nnxnny])/tauiT+
      (ci[(nnx-1)+j*nnx+(0+1)*nnxnny]+ci[(nnx-1)+j*nnx+(0+1)*nnxnny]-2*ci[(nnx-1)+j*nnx+0*nnxnny])/tauiT;
    Ici[0+j*nnx+(nnz-1)*nnxnny]=(ci[0+j*nnx+(nnz-1)*nnxnny+1]+ci[0+j*nnx+(nnz-1)*nnxnny+1]-2*ci[0+j*nnx+(nnz-1)*nnxnny])/tauiL+
      (ci[0+(j+1)*nnx+(nnz-1)*nnxnny]+ci[0+(j-1)*nnx+(nnz-1)*nnxnny]-2*ci[0+j*nnx+(nnz-1)*nnxnny])/tauiT+
      (ci[0+j*nnx+((nnz-1)-1)*nnxnny]+ci[0+j*nnx+((nnz-1)-1)*nnxnny]-2*ci[0+j*nnx+(nnz-1)*nnxnny])/tauiT;
    Ici[(nnx-1)+j*nnx+(nnz-1)*nnxnny]=(ci[(nnx-1)+j*nnx+(nnz-1)*nnxnny-1]+ci[(nnx-1)+j*nnx+(nnz-1)*nnxnny-1]-2*ci[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/tauiL+
      (ci[(nnx-1)+(j+1)*nnx+(nnz-1)*nnxnny]+ci[(nnx-1)+(j-1)*nnx+(nnz-1)*nnxnny]-2*ci[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/tauiT+
      (ci[(nnx-1)+j*nnx+((nnz-1)-1)*nnxnny]+ci[(nnx-1)+j*nnx+((nnz-1)-1)*nnxnny]-2*ci[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/tauiT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      Ici[0+j*nnx+k*nnxnny]=(ci[0+j*nnx+k*nnxnny+1]+ci[0+j*nnx+k*nnxnny+1]-2*ci[0+j*nnx+k*nnxnny])/tauiL+
        (ci[0+(j+1)*nnx+k*nnxnny]+ci[0+(j-1)*nnx+k*nnxnny]-2*ci[0+j*nnx+k*nnxnny])/tauiT+
        (ci[0+j*nnx+(k+1)*nnxnny]+ci[0+j*nnx+(k-1)*nnxnny]-2*ci[0+j*nnx+k*nnxnny])/tauiT;
      Ici[(nnx-1)+j*nnx+k*nnxnny]=(ci[(nnx-1)+j*nnx+k*nnxnny-1]+ci[(nnx-1)+j*nnx+k*nnxnny-1]-2*ci[(nnx-1)+j*nnx+k*nnxnny])/tauiL+
        (ci[(nnx-1)+(j+1)*nnx+k*nnxnny]+ci[(nnx-1)+(j-1)*nnx+k*nnxnny]-2*ci[(nnx-1)+j*nnx+k*nnxnny])/tauiT+
        (ci[(nnx-1)+j*nnx+(k+1)*nnxnny]+ci[(nnx-1)+j*nnx+(k-1)*nnxnny]-2*ci[(nnx-1)+j*nnx+k*nnxnny])/tauiT;
    }
  }
  //y fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=1;i<(nnx-1);i++)
  {
    Ici[i+0*nnx+0*nnxnny]=(ci[i+0*nnx+0*nnxnny+1]+ci[i+0*nnx+0*nnxnny-1]-2*ci[i+0*nnx+0*nnxnny])/tauiL+
      (ci[i+(0+1)*nnx+0*nnxnny]+ci[i+(0+1)*nnx+0*nnxnny]-2*ci[i+0*nnx+0*nnxnny])/tauiT+
      (ci[i+0*nnx+(0+1)*nnxnny]+ci[i+0*nnx+(0+1)*nnxnny]-2*ci[i+0*nnx+0*nnxnny])/tauiT;
    Ici[i+(nny-1)*nnx+0*nnxnny]=(ci[i+(nny-1)*nnx+0*nnxnny+1]+ci[i+(nny-1)*nnx+0*nnxnny-1]-2*ci[i+(nny-1)*nnx+0*nnxnny])/tauiL+
      (ci[i+((nny-1)-1)*nnx+0*nnxnny]+ci[i+((nny-1)-1)*nnx+0*nnxnny]-2*ci[i+(nny-1)*nnx+0*nnxnny])/tauiT+
      (ci[i+(nny-1)*nnx+(0+1)*nnxnny]+ci[i+(nny-1)*nnx+(0+1)*nnxnny]-2*ci[i+(nny-1)*nnx+0*nnxnny])/tauiT;
    Ici[i+0*nnx+(nnz-1)*nnxnny]=(ci[i+0*nnx+(nnz-1)*nnxnny+1]+ci[i+0*nnx+(nnz-1)*nnxnny-1]-2*ci[i+0*nnx+(nnz-1)*nnxnny])/tauiL+
      (ci[i+(0+1)*nnx+(nnz-1)*nnxnny]+ci[i+(0+1)*nnx+(nnz-1)*nnxnny]-2*ci[i+0*nnx+(nnz-1)*nnxnny])/tauiT+
      (ci[i+0*nnx+((nnz-1)-1)*nnxnny]+ci[i+0*nnx+((nnz-1)-1)*nnxnny]-2*ci[i+0*nnx+(nnz-1)*nnxnny])/tauiT;
    Ici[i+(nny-1)*nnx+(nnz-1)*nnxnny]=(ci[i+(nny-1)*nnx+(nnz-1)*nnxnny+1]+ci[i+(nny-1)*nnx+(nnz-1)*nnxnny-1]-2*ci[i+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiL+
      (ci[i+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+ci[i+((nny-1)-1)*nnx+(nnz-1)*nnxnny]-2*ci[i+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiT+
      (ci[i+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+ci[i+(nny-1)*nnx+((nnz-1)-1)*nnxnny]-2*ci[i+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      Ici[i+0*nnx+k*nnxnny]=(ci[i+0*nnx+k*nnxnny+1]+ci[i+0*nnx+k*nnxnny-1]-2*ci[i+0*nnx+k*nnxnny])/tauiL+
        (ci[i+(0+1)*nnx+k*nnxnny]+ci[i+(0+1)*nnx+k*nnxnny]-2*ci[i+0*nnx+k*nnxnny])/tauiT+
        (ci[i+0*nnx+(k+1)*nnxnny]+ci[i+0*nnx+(k-1)*nnxnny]-2*ci[i+0*nnx+k*nnxnny])/tauiT;
      Ici[i+(nny-1)*nnx+k*nnxnny]=(ci[i+(nny-1)*nnx+k*nnxnny+1]+ci[i+(nny-1)*nnx+k*nnxnny-1]-2*ci[i+(nny-1)*nnx+k*nnxnny])/tauiL+
        (ci[i+((nny-1)-1)*nnx+k*nnxnny]+ci[i+((nny-1)-1)*nnx+k*nnxnny]-2*ci[i+(nny-1)*nnx+k*nnxnny])/tauiT+
        (ci[i+(nny-1)*nnx+(k+1)*nnxnny]+ci[i+(nny-1)*nnx+(k-1)*nnxnny]-2*ci[i+(nny-1)*nnx+k*nnxnny])/tauiT;
    }
    //z fixed
#pragma ivdep
#pragma vector always
    for (int j=1;j<nny-1;j++)
    {
      Ici[i+j*nnx+0*nnxnny]=(ci[i+j*nnx+0*nnxnny+1]+ci[i+j*nnx+0*nnxnny-1]-2*ci[i+j*nnx+0*nnxnny])/tauiL+
        (ci[i+(j+1)*nnx+0*nnxnny]+ci[i+(j-1)*nnx+0*nnxnny]-2*ci[i+j*nnx+0*nnxnny])/tauiT+
        (ci[i+j*nnx+(0+1)*nnxnny]+ci[i+j*nnx+(0+1)*nnxnny]-2*ci[i+j*nnx+0*nnxnny])/tauiT;
      Ici[i+j*nnx+(nnz-1)*nnxnny]=(ci[i+j*nnx+(nnz-1)*nnxnny+1]+ci[i+j*nnx+(nnz-1)*nnxnny-1]-2*ci[i+j*nnx+(nnz-1)*nnxnny])/tauiL+
        (ci[i+(j+1)*nnx+(nnz-1)*nnxnny]+ci[i+(j-1)*nnx+(nnz-1)*nnxnny]-2*ci[i+j*nnx+(nnz-1)*nnxnny])/tauiT+
        (ci[i+j*nnx+((nnz-1)-1)*nnxnny]+ci[i+j*nnx+((nnz-1)-1)*nnxnny]-2*ci[i+j*nnx+(nnz-1)*nnxnny])/tauiT;
    }
  }
#pragma ivdep
#pragma vector always
  for (int k=1;k<nnz-1;k++)
  {
    Ici[0+0*nnx+k*nnxnny]=(ci[0+0*nnx+k*nnxnny+1]+ci[0+0*nnx+k*nnxnny+1]-2*ci[0+0*nnx+k*nnxnny])/tauiL+
      (ci[0+(0+1)*nnx+k*nnxnny]+ci[0+(0+1)*nnx+k*nnxnny]-2*ci[0+0*nnx+k*nnxnny])/tauiT+
      (ci[0+0*nnx+(k+1)*nnxnny]+ci[0+0*nnx+(k-1)*nnxnny]-2*ci[0+0*nnx+k*nnxnny])/tauiT;
    Ici[0+(nny-1)*nnx+k*nnxnny]=(ci[0+(nny-1)*nnx+k*nnxnny+1]+ci[0+(nny-1)*nnx+k*nnxnny+1]-2*ci[0+(nny-1)*nnx+k*nnxnny])/tauiL+
      (ci[0+((nny-1)-1)*nnx+k*nnxnny]+ci[0+((nny-1)-1)*nnx+k*nnxnny]-2*ci[0+(nny-1)*nnx+k*nnxnny])/tauiT+
      (ci[0+(nny-1)*nnx+(k+1)*nnxnny]+ci[0+(nny-1)*nnx+(k-1)*nnxnny]-2*ci[0+(nny-1)*nnx+k*nnxnny])/tauiT;
    Ici[(nnx-1)+0*nnx+k*nnxnny]=(ci[(nnx-1)+0*nnx+k*nnxnny-1]+ci[(nnx-1)+0*nnx+k*nnxnny-1]-2*ci[(nnx-1)+0*nnx+k*nnxnny])/tauiL+
      (ci[(nnx-1)+(0+1)*nnx+k*nnxnny]+ci[(nnx-1)+(0+1)*nnx+k*nnxnny]-2*ci[(nnx-1)+0*nnx+k*nnxnny])/tauiT+
      (ci[(nnx-1)+0*nnx+(k+1)*nnxnny]+ci[(nnx-1)+0*nnx+(k-1)*nnxnny]-2*ci[(nnx-1)+0*nnx+k*nnxnny])/tauiT;
    Ici[(nnx-1)+(nny-1)*nnx+k*nnxnny]=(ci[(nnx-1)+(nny-1)*nnx+k*nnxnny-1]+ci[(nnx-1)+(nny-1)*nnx+k*nnxnny-1]-2*ci[(nnx-1)+(nny-1)*nnx+k*nnxnny])/tauiL+
      (ci[(nnx-1)+((nny-1)-1)*nnx+k*nnxnny]+ci[(nnx-1)+((nny-1)-1)*nnx+k*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+k*nnxnny])/tauiT+
      (ci[(nnx-1)+(nny-1)*nnx+(k+1)*nnxnny]+ci[(nnx-1)+(nny-1)*nnx+(k-1)*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+k*nnxnny])/tauiT;
  }

#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int k=1;k<nnz-1;k++)
  {
    for (int j=1;j<nny-1;j++)
    {
#pragma ivdep
#pragma vector always
      for (int i=1;i<nnx-1;i++)
      {
        Ici[i+j*nnx+k*nnxnny]=(ci[i+j*nnx+k*nnxnny+1]+ci[i+j*nnx+k*nnxnny-1]-2*ci[i+j*nnx+k*nnxnny])/tauiL+
          (ci[i+(j+1)*nnx+k*nnxnny]+ci[i+(j-1)*nnx+k*nnxnny]-2*ci[i+j*nnx+k*nnxnny])/tauiT+
          (ci[i+j*nnx+(k+1)*nnxnny]+ci[i+j*nnx+(k-1)*nnxnny]-2*ci[i+j*nnx+k*nnxnny])/tauiT;
      }
    }
  }
}
void CSubcell::computeIcnsr(void)
{
#ifdef ___PERIODIC
  //boundary
  Icnsr[0+0*nnx+0*nnxnny]=(cnsr[0+0*nnx+0*nnxnny+1]+cnsr[(nnx-1)+0*nnx+0*nnxnny]-2*cnsr[0+0*nnx+0*nnxnny])/taunsrL+
    (cnsr[0+(0+1)*nnx+0*nnxnny]+cnsr[0+(nny-1)*nnx+0*nnxnny]-2*cnsr[0+0*nnx+0*nnxnny])/taunsrT+
    (cnsr[0+0*nnx+(0+1)*nnxnny]+cnsr[0+0*nnx+(nnz-1)*nnxnny]-2*cnsr[0+0*nnx+0*nnxnny])/taunsrT;
  Icnsr[0+(nny-1)*nnx+(nnz-1)*nnxnny]=(cnsr[0+(nny-1)*nnx+(nnz-1)*nnxnny+1]+cnsr[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[0+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrL+
    (cnsr[0+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cnsr[0+(0)*nnx+(nnz-1)*nnxnny]-2*cnsr[0+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrT+
    (cnsr[0+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cnsr[0+(nny-1)*nnx+(0)*nnxnny]-2*cnsr[0+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrT;
  Icnsr[0+(nny-1)*nnx+0*nnxnny]=(cnsr[0+(nny-1)*nnx+0*nnxnny+1]+cnsr[(nnx-1)+(nny-1)*nnx+0*nnxnny]-2*cnsr[0+(nny-1)*nnx+0*nnxnny])/taunsrL+
    (cnsr[0+((nny-1)-1)*nnx+0*nnxnny]+cnsr[0+(0)*nnx+0*nnxnny]-2*cnsr[0+(nny-1)*nnx+0*nnxnny])/taunsrT+
    (cnsr[0+(nny-1)*nnx+(0+1)*nnxnny]+cnsr[0+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[0+(nny-1)*nnx+0*nnxnny])/taunsrT;
  Icnsr[0+0*nnx+(nnz-1)*nnxnny]=(cnsr[0+0*nnx+(nnz-1)*nnxnny+1]+cnsr[(nnx-1)+0*nnx+(nnz-1)*nnxnny]-2*cnsr[0+0*nnx+(nnz-1)*nnxnny])/taunsrL+
    (cnsr[0+(0+1)*nnx+(nnz-1)*nnxnny]+cnsr[0+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[0+0*nnx+(nnz-1)*nnxnny])/taunsrT+
    (cnsr[0+0*nnx+((nnz-1)-1)*nnxnny]+cnsr[0+0*nnx+(0)*nnxnny]-2*cnsr[0+0*nnx+(nnz-1)*nnxnny])/taunsrT;
  Icnsr[(nnx-1)+(nny-1)*nnx+0*nnxnny]=(cnsr[(nnx-1)+(nny-1)*nnx+0*nnxnny-1]+cnsr[(0)+(nny-1)*nnx+0*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+0*nnxnny])/taunsrL+
    (cnsr[(nnx-1)+((nny-1)-1)*nnx+0*nnxnny]+cnsr[(nnx-1)+(0)*nnx+0*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+0*nnxnny])/taunsrT+
    (cnsr[(nnx-1)+(nny-1)*nnx+(0+1)*nnxnny]+cnsr[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+0*nnxnny])/taunsrT;
  Icnsr[(nnx-1)+0*nnx+(nnz-1)*nnxnny]=(cnsr[(nnx-1)+0*nnx+(nnz-1)*nnxnny-1]+cnsr[(0)+0*nnx+(nnz-1)*nnxnny]-2*cnsr[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/taunsrL+
    (cnsr[(nnx-1)+(0+1)*nnx+(nnz-1)*nnxnny]+cnsr[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/taunsrT+
    (cnsr[(nnx-1)+0*nnx+((nnz-1)-1)*nnxnny]+cnsr[(nnx-1)+0*nnx+(0)*nnxnny]-2*cnsr[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/taunsrT;
  Icnsr[(nnx-1)+0*nnx+0*nnxnny]=(cnsr[(nnx-1)+0*nnx+0*nnxnny-1]+cnsr[(0)+0*nnx+0*nnxnny]-2*cnsr[(nnx-1)+0*nnx+0*nnxnny])/taunsrL+
    (cnsr[(nnx-1)+(0+1)*nnx+0*nnxnny]+cnsr[(nnx-1)+(nny-1)*nnx+0*nnxnny]-2*cnsr[(nnx-1)+0*nnx+0*nnxnny])/taunsrT+
    (cnsr[(nnx-1)+0*nnx+(0+1)*nnxnny]+cnsr[(nnx-1)+0*nnx+(nnz-1)*nnxnny]-2*cnsr[(nnx-1)+0*nnx+0*nnxnny])/taunsrT;
  Icnsr[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]=(cnsr[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny-1]+cnsr[(0)+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrL+
    (cnsr[(nnx-1)+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cnsr[(nnx-1)+(0)*nnx+(nnz-1)*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrT+
    (cnsr[(nnx-1)+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cnsr[(nnx-1)+(nny-1)*nnx+(0)*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrT;
  //x fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int j=1;j<nny-1;j++)
  {
    Icnsr[0+j*nnx+0*nnxnny]=(cnsr[0+j*nnx+0*nnxnny+1]+cnsr[(nnx-1)+j*nnx+0*nnxnny]-2*cnsr[0+j*nnx+0*nnxnny])/taunsrL+
      (cnsr[0+(j+1)*nnx+0*nnxnny]+cnsr[0+(j-1)*nnx+0*nnxnny]-2*cnsr[0+j*nnx+0*nnxnny])/taunsrT+
      (cnsr[0+j*nnx+(0+1)*nnxnny]+cnsr[0+j*nnx+(nnz-1)*nnxnny]-2*cnsr[0+j*nnx+0*nnxnny])/taunsrT;
    Icnsr[(nnx-1)+j*nnx+0*nnxnny]=(cnsr[(nnx-1)+j*nnx+0*nnxnny-1]+cnsr[(0)+j*nnx+0*nnxnny]-2*cnsr[(nnx-1)+j*nnx+0*nnxnny])/taunsrL+
      (cnsr[(nnx-1)+(j+1)*nnx+0*nnxnny]+cnsr[(nnx-1)+(j-1)*nnx+0*nnxnny]-2*cnsr[(nnx-1)+j*nnx+0*nnxnny])/taunsrT+
      (cnsr[(nnx-1)+j*nnx+(0+1)*nnxnny]+cnsr[(nnx-1)+j*nnx+(nnz-1)*nnxnny]-2*cnsr[(nnx-1)+j*nnx+0*nnxnny])/taunsrT;
    Icnsr[0+j*nnx+(nnz-1)*nnxnny]=(cnsr[0+j*nnx+(nnz-1)*nnxnny+1]+cnsr[(nnx-1)+j*nnx+(nnz-1)*nnxnny]-2*cnsr[0+j*nnx+(nnz-1)*nnxnny])/taunsrL+
      (cnsr[0+(j+1)*nnx+(nnz-1)*nnxnny]+cnsr[0+(j-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[0+j*nnx+(nnz-1)*nnxnny])/taunsrT+
      (cnsr[0+j*nnx+((nnz-1)-1)*nnxnny]+cnsr[0+j*nnx+(0)*nnxnny]-2*cnsr[0+j*nnx+(nnz-1)*nnxnny])/taunsrT;
    Icnsr[(nnx-1)+j*nnx+(nnz-1)*nnxnny]=(cnsr[(nnx-1)+j*nnx+(nnz-1)*nnxnny-1]+cnsr[(0)+j*nnx+(nnz-1)*nnxnny]-2*cnsr[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/taunsrL+
      (cnsr[(nnx-1)+(j+1)*nnx+(nnz-1)*nnxnny]+cnsr[(nnx-1)+(j-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/taunsrT+
      (cnsr[(nnx-1)+j*nnx+((nnz-1)-1)*nnxnny]+cnsr[(nnx-1)+j*nnx+(0)*nnxnny]-2*cnsr[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/taunsrT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      Icnsr[0+j*nnx+k*nnxnny]=(cnsr[0+j*nnx+k*nnxnny+1]+cnsr[(nnx-1)+j*nnx+k*nnxnny]-2*cnsr[0+j*nnx+k*nnxnny])/taunsrL+
        (cnsr[0+(j+1)*nnx+k*nnxnny]+cnsr[0+(j-1)*nnx+k*nnxnny]-2*cnsr[0+j*nnx+k*nnxnny])/taunsrT+
        (cnsr[0+j*nnx+(k+1)*nnxnny]+cnsr[0+j*nnx+(k-1)*nnxnny]-2*cnsr[0+j*nnx+k*nnxnny])/taunsrT;
      Icnsr[(nnx-1)+j*nnx+k*nnxnny]=(cnsr[(nnx-1)+j*nnx+k*nnxnny-1]+cnsr[(0)+j*nnx+k*nnxnny]-2*cnsr[(nnx-1)+j*nnx+k*nnxnny])/taunsrL+
        (cnsr[(nnx-1)+(j+1)*nnx+k*nnxnny]+cnsr[(nnx-1)+(j-1)*nnx+k*nnxnny]-2*cnsr[(nnx-1)+j*nnx+k*nnxnny])/taunsrT+
        (cnsr[(nnx-1)+j*nnx+(k+1)*nnxnny]+cnsr[(nnx-1)+j*nnx+(k-1)*nnxnny]-2*cnsr[(nnx-1)+j*nnx+k*nnxnny])/taunsrT;
    }
  }
  //y fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=1;i<(nnx-1);i++)
  {
    Icnsr[i+0*nnx+0*nnxnny]=(cnsr[i+0*nnx+0*nnxnny+1]+cnsr[i+0*nnx+0*nnxnny-1]-2*cnsr[i+0*nnx+0*nnxnny])/taunsrL+
      (cnsr[i+(0+1)*nnx+0*nnxnny]+cnsr[i+(nny-1)*nnx+0*nnxnny]-2*cnsr[i+0*nnx+0*nnxnny])/taunsrT+
      (cnsr[i+0*nnx+(0+1)*nnxnny]+cnsr[i+0*nnx+(nnz-1)*nnxnny]-2*cnsr[i+0*nnx+0*nnxnny])/taunsrT;
    Icnsr[i+(nny-1)*nnx+0*nnxnny]=(cnsr[i+(nny-1)*nnx+0*nnxnny+1]+cnsr[i+(nny-1)*nnx+0*nnxnny-1]-2*cnsr[i+(nny-1)*nnx+0*nnxnny])/taunsrL+
      (cnsr[i+((nny-1)-1)*nnx+0*nnxnny]+cnsr[i+(0)*nnx+0*nnxnny]-2*cnsr[i+(nny-1)*nnx+0*nnxnny])/taunsrT+
      (cnsr[i+(nny-1)*nnx+(0+1)*nnxnny]+cnsr[i+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[i+(nny-1)*nnx+0*nnxnny])/taunsrT;
    Icnsr[i+0*nnx+(nnz-1)*nnxnny]=(cnsr[i+0*nnx+(nnz-1)*nnxnny+1]+cnsr[i+0*nnx+(nnz-1)*nnxnny-1]-2*cnsr[i+0*nnx+(nnz-1)*nnxnny])/taunsrL+
      (cnsr[i+(0+1)*nnx+(nnz-1)*nnxnny]+cnsr[i+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[i+0*nnx+(nnz-1)*nnxnny])/taunsrT+
      (cnsr[i+0*nnx+((nnz-1)-1)*nnxnny]+cnsr[i+0*nnx+(0)*nnxnny]-2*cnsr[i+0*nnx+(nnz-1)*nnxnny])/taunsrT;
    Icnsr[i+(nny-1)*nnx+(nnz-1)*nnxnny]=(cnsr[i+(nny-1)*nnx+(nnz-1)*nnxnny+1]+cnsr[i+(nny-1)*nnx+(nnz-1)*nnxnny-1]-2*cnsr[i+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrL+
      (cnsr[i+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cnsr[i+(0)*nnx+(nnz-1)*nnxnny]-2*cnsr[i+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrT+
      (cnsr[i+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cnsr[i+(nny-1)*nnx+(0)*nnxnny]-2*cnsr[i+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      Icnsr[i+0*nnx+k*nnxnny]=(cnsr[i+0*nnx+k*nnxnny+1]+cnsr[i+0*nnx+k*nnxnny-1]-2*cnsr[i+0*nnx+k*nnxnny])/taunsrL+
        (cnsr[i+(0+1)*nnx+k*nnxnny]+cnsr[i+(nny-1)*nnx+k*nnxnny]-2*cnsr[i+0*nnx+k*nnxnny])/taunsrT+
        (cnsr[i+0*nnx+(k+1)*nnxnny]+cnsr[i+0*nnx+(k-1)*nnxnny]-2*cnsr[i+0*nnx+k*nnxnny])/taunsrT;
      Icnsr[i+(nny-1)*nnx+k*nnxnny]=(cnsr[i+(nny-1)*nnx+k*nnxnny+1]+cnsr[i+(nny-1)*nnx+k*nnxnny-1]-2*cnsr[i+(nny-1)*nnx+k*nnxnny])/taunsrL+
        (cnsr[i+((nny-1)-1)*nnx+k*nnxnny]+cnsr[i+(0)*nnx+k*nnxnny]-2*cnsr[i+(nny-1)*nnx+k*nnxnny])/taunsrT+
        (cnsr[i+(nny-1)*nnx+(k+1)*nnxnny]+cnsr[i+(nny-1)*nnx+(k-1)*nnxnny]-2*cnsr[i+(nny-1)*nnx+k*nnxnny])/taunsrT;
    }
    //z fixed
#pragma ivdep
#pragma vector always
    for (int j=1;j<nny-1;j++)
    {
      Icnsr[i+j*nnx+0*nnxnny]=(cnsr[i+j*nnx+0*nnxnny+1]+cnsr[i+j*nnx+0*nnxnny-1]-2*cnsr[i+j*nnx+0*nnxnny])/taunsrL+
        (cnsr[i+(j+1)*nnx+0*nnxnny]+cnsr[i+(j-1)*nnx+0*nnxnny]-2*cnsr[i+j*nnx+0*nnxnny])/taunsrT+
        (cnsr[i+j*nnx+(0+1)*nnxnny]+cnsr[i+j*nnx+(nnz-1)*nnxnny]-2*cnsr[i+j*nnx+0*nnxnny])/taunsrT;
      Icnsr[i+j*nnx+(nnz-1)*nnxnny]=(cnsr[i+j*nnx+(nnz-1)*nnxnny+1]+cnsr[i+j*nnx+(nnz-1)*nnxnny-1]-2*cnsr[i+j*nnx+(nnz-1)*nnxnny])/taunsrL+
        (cnsr[i+(j+1)*nnx+(nnz-1)*nnxnny]+cnsr[i+(j-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[i+j*nnx+(nnz-1)*nnxnny])/taunsrT+
        (cnsr[i+j*nnx+((nnz-1)-1)*nnxnny]+cnsr[i+j*nnx+(0)*nnxnny]-2*cnsr[i+j*nnx+(nnz-1)*nnxnny])/taunsrT;
    }
  }
#pragma ivdep
#pragma vector always
  for (int k=1;k<nnz-1;k++)
  {
    Icnsr[0+0*nnx+k*nnxnny]=(cnsr[0+0*nnx+k*nnxnny+1]+cnsr[(nnx-1)+0*nnx+k*nnxnny]-2*cnsr[0+0*nnx+k*nnxnny])/taunsrL+
      (cnsr[0+(0+1)*nnx+k*nnxnny]+cnsr[0+(nny-1)*nnx+k*nnxnny]-2*cnsr[0+0*nnx+k*nnxnny])/taunsrT+
      (cnsr[0+0*nnx+(k+1)*nnxnny]+cnsr[0+0*nnx+(k-1)*nnxnny]-2*cnsr[0+0*nnx+k*nnxnny])/taunsrT;
    Icnsr[0+(nny-1)*nnx+k*nnxnny]=(cnsr[0+(nny-1)*nnx+k*nnxnny+1]+cnsr[(nnx-1)+(nny-1)*nnx+k*nnxnny]-2*cnsr[0+(nny-1)*nnx+k*nnxnny])/taunsrL+
      (cnsr[0+((nny-1)-1)*nnx+k*nnxnny]+cnsr[0+(0)*nnx+k*nnxnny]-2*cnsr[0+(nny-1)*nnx+k*nnxnny])/taunsrT+
      (cnsr[0+(nny-1)*nnx+(k+1)*nnxnny]+cnsr[0+(nny-1)*nnx+(k-1)*nnxnny]-2*cnsr[0+(nny-1)*nnx+k*nnxnny])/taunsrT;
    Icnsr[(nnx-1)+0*nnx+k*nnxnny]=(cnsr[(nnx-1)+0*nnx+k*nnxnny-1]+cnsr[(0)+0*nnx+k*nnxnny]-2*cnsr[(nnx-1)+0*nnx+k*nnxnny])/taunsrL+
      (cnsr[(nnx-1)+(0+1)*nnx+k*nnxnny]+cnsr[(nnx-1)+(nny-1)*nnx+k*nnxnny]-2*cnsr[(nnx-1)+0*nnx+k*nnxnny])/taunsrT+
      (cnsr[(nnx-1)+0*nnx+(k+1)*nnxnny]+cnsr[(nnx-1)+0*nnx+(k-1)*nnxnny]-2*cnsr[(nnx-1)+0*nnx+k*nnxnny])/taunsrT;
    Icnsr[(nnx-1)+(nny-1)*nnx+k*nnxnny]=(cnsr[(nnx-1)+(nny-1)*nnx+k*nnxnny-1]+cnsr[(0)+(nny-1)*nnx+k*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+k*nnxnny])/taunsrL+
      (cnsr[(nnx-1)+((nny-1)-1)*nnx+k*nnxnny]+cnsr[(nnx-1)+(0)*nnx+k*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+k*nnxnny])/taunsrT+
      (cnsr[(nnx-1)+(nny-1)*nnx+(k+1)*nnxnny]+cnsr[(nnx-1)+(nny-1)*nnx+(k-1)*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+k*nnxnny])/taunsrT;
  }
#else
  //boundary

  Icnsr[0+0*nnx+0*nnxnny]=(cnsr[0+0*nnx+0*nnxnny+1]+cnsr[0+0*nnx+0*nnxnny+1]-2*cnsr[0+0*nnx+0*nnxnny])/taunsrL+
    (cnsr[0+(0+1)*nnx+0*nnxnny]+cnsr[0+(0+1)*nnx+0*nnxnny]-2*cnsr[0+0*nnx+0*nnxnny])/taunsrT+
    (cnsr[0+0*nnx+(0+1)*nnxnny]+cnsr[0+0*nnx+(0+1)*nnxnny]-2*cnsr[0+0*nnx+0*nnxnny])/taunsrT;
  Icnsr[0+(nny-1)*nnx+(nnz-1)*nnxnny]=(cnsr[0+(nny-1)*nnx+(nnz-1)*nnxnny+1]+cnsr[0+(nny-1)*nnx+(nnz-1)*nnxnny+1]-2*cnsr[0+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrL+
    (cnsr[0+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cnsr[0+((nny-1)-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[0+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrT+
    (cnsr[0+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cnsr[0+(nny-1)*nnx+((nnz-1)-1)*nnxnny]-2*cnsr[0+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrT;
  Icnsr[0+(nny-1)*nnx+0*nnxnny]=(cnsr[0+(nny-1)*nnx+0*nnxnny+1]+cnsr[0+(nny-1)*nnx+0*nnxnny+1]-2*cnsr[0+(nny-1)*nnx+0*nnxnny])/taunsrL+
    (cnsr[0+((nny-1)-1)*nnx+0*nnxnny]+cnsr[0+((nny-1)-1)*nnx+0*nnxnny]-2*cnsr[0+(nny-1)*nnx+0*nnxnny])/taunsrT+
    (cnsr[0+(nny-1)*nnx+(0+1)*nnxnny]+cnsr[0+(nny-1)*nnx+(0+1)*nnxnny]-2*cnsr[0+(nny-1)*nnx+0*nnxnny])/taunsrT;
  Icnsr[0+0*nnx+(nnz-1)*nnxnny]=(cnsr[0+0*nnx+(nnz-1)*nnxnny+1]+cnsr[0+0*nnx+(nnz-1)*nnxnny+1]-2*cnsr[0+0*nnx+(nnz-1)*nnxnny])/taunsrL+
    (cnsr[0+(0+1)*nnx+(nnz-1)*nnxnny]+cnsr[0+(0+1)*nnx+(nnz-1)*nnxnny]-2*cnsr[0+0*nnx+(nnz-1)*nnxnny])/taunsrT+
    (cnsr[0+0*nnx+((nnz-1)-1)*nnxnny]+cnsr[0+0*nnx+((nnz-1)-1)*nnxnny]-2*cnsr[0+0*nnx+(nnz-1)*nnxnny])/taunsrT;
  Icnsr[(nnx-1)+(nny-1)*nnx+0*nnxnny]=(cnsr[(nnx-1)+(nny-1)*nnx+0*nnxnny-1]+cnsr[(nnx-1)+(nny-1)*nnx+0*nnxnny-1]-2*cnsr[(nnx-1)+(nny-1)*nnx+0*nnxnny])/taunsrL+
    (cnsr[(nnx-1)+((nny-1)-1)*nnx+0*nnxnny]+cnsr[(nnx-1)+((nny-1)-1)*nnx+0*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+0*nnxnny])/taunsrT+
    (cnsr[(nnx-1)+(nny-1)*nnx+(0+1)*nnxnny]+cnsr[(nnx-1)+(nny-1)*nnx+(0+1)*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+0*nnxnny])/taunsrT;
  Icnsr[(nnx-1)+0*nnx+(nnz-1)*nnxnny]=(cnsr[(nnx-1)+0*nnx+(nnz-1)*nnxnny-1]+cnsr[(nnx-1)+0*nnx+(nnz-1)*nnxnny-1]-2*cnsr[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/taunsrL+
    (cnsr[(nnx-1)+(0+1)*nnx+(nnz-1)*nnxnny]+cnsr[(nnx-1)+(0+1)*nnx+(nnz-1)*nnxnny]-2*cnsr[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/taunsrT+
    (cnsr[(nnx-1)+0*nnx+((nnz-1)-1)*nnxnny]+cnsr[(nnx-1)+0*nnx+((nnz-1)-1)*nnxnny]-2*cnsr[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/taunsrT;
  Icnsr[(nnx-1)+0*nnx+0*nnxnny]=(cnsr[(nnx-1)+0*nnx+0*nnxnny-1]+cnsr[(nnx-1)+0*nnx+0*nnxnny-1]-2*cnsr[(nnx-1)+0*nnx+0*nnxnny])/taunsrL+
    (cnsr[(nnx-1)+(0+1)*nnx+0*nnxnny]+cnsr[(nnx-1)+(0+1)*nnx+0*nnxnny]-2*cnsr[(nnx-1)+0*nnx+0*nnxnny])/taunsrT+
    (cnsr[(nnx-1)+0*nnx+(0+1)*nnxnny]+cnsr[(nnx-1)+0*nnx+(0+1)*nnxnny]-2*cnsr[(nnx-1)+0*nnx+0*nnxnny])/taunsrT;
  Icnsr[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]=(cnsr[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny-1]+cnsr[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny-1]-2*cnsr[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrL+
    (cnsr[(nnx-1)+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cnsr[(nnx-1)+((nny-1)-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrT+
    (cnsr[(nnx-1)+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cnsr[(nnx-1)+(nny-1)*nnx+((nnz-1)-1)*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrT;
  //x fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int j=1;j<nny-1;j++)
  {
    Icnsr[0+j*nnx+0*nnxnny]=(cnsr[0+j*nnx+0*nnxnny+1]+cnsr[0+j*nnx+0*nnxnny+1]-2*cnsr[0+j*nnx+0*nnxnny])/taunsrL+
      (cnsr[0+(j+1)*nnx+0*nnxnny]+cnsr[0+(j-1)*nnx+0*nnxnny]-2*cnsr[0+j*nnx+0*nnxnny])/taunsrT+
      (cnsr[0+j*nnx+(0+1)*nnxnny]+cnsr[0+j*nnx+(0+1)*nnxnny]-2*cnsr[0+j*nnx+0*nnxnny])/taunsrT;
    Icnsr[(nnx-1)+j*nnx+0*nnxnny]=(cnsr[(nnx-1)+j*nnx+0*nnxnny-1]+cnsr[(nnx-1)+j*nnx+0*nnxnny-1]-2*cnsr[(nnx-1)+j*nnx+0*nnxnny])/taunsrL+
      (cnsr[(nnx-1)+(j+1)*nnx+0*nnxnny]+cnsr[(nnx-1)+(j-1)*nnx+0*nnxnny]-2*cnsr[(nnx-1)+j*nnx+0*nnxnny])/taunsrT+
      (cnsr[(nnx-1)+j*nnx+(0+1)*nnxnny]+cnsr[(nnx-1)+j*nnx+(0+1)*nnxnny]-2*cnsr[(nnx-1)+j*nnx+0*nnxnny])/taunsrT;
    Icnsr[0+j*nnx+(nnz-1)*nnxnny]=(cnsr[0+j*nnx+(nnz-1)*nnxnny+1]+cnsr[0+j*nnx+(nnz-1)*nnxnny+1]-2*cnsr[0+j*nnx+(nnz-1)*nnxnny])/taunsrL+
      (cnsr[0+(j+1)*nnx+(nnz-1)*nnxnny]+cnsr[0+(j-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[0+j*nnx+(nnz-1)*nnxnny])/taunsrT+
      (cnsr[0+j*nnx+((nnz-1)-1)*nnxnny]+cnsr[0+j*nnx+((nnz-1)-1)*nnxnny]-2*cnsr[0+j*nnx+(nnz-1)*nnxnny])/taunsrT;
    Icnsr[(nnx-1)+j*nnx+(nnz-1)*nnxnny]=(cnsr[(nnx-1)+j*nnx+(nnz-1)*nnxnny-1]+cnsr[(nnx-1)+j*nnx+(nnz-1)*nnxnny-1]-2*cnsr[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/taunsrL+
      (cnsr[(nnx-1)+(j+1)*nnx+(nnz-1)*nnxnny]+cnsr[(nnx-1)+(j-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/taunsrT+
      (cnsr[(nnx-1)+j*nnx+((nnz-1)-1)*nnxnny]+cnsr[(nnx-1)+j*nnx+((nnz-1)-1)*nnxnny]-2*cnsr[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/taunsrT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      Icnsr[0+j*nnx+k*nnxnny]=(cnsr[0+j*nnx+k*nnxnny+1]+cnsr[0+j*nnx+k*nnxnny+1]-2*cnsr[0+j*nnx+k*nnxnny])/taunsrL+
        (cnsr[0+(j+1)*nnx+k*nnxnny]+cnsr[0+(j-1)*nnx+k*nnxnny]-2*cnsr[0+j*nnx+k*nnxnny])/taunsrT+
        (cnsr[0+j*nnx+(k+1)*nnxnny]+cnsr[0+j*nnx+(k-1)*nnxnny]-2*cnsr[0+j*nnx+k*nnxnny])/taunsrT;
      Icnsr[(nnx-1)+j*nnx+k*nnxnny]=(cnsr[(nnx-1)+j*nnx+k*nnxnny-1]+cnsr[(nnx-1)+j*nnx+k*nnxnny-1]-2*cnsr[(nnx-1)+j*nnx+k*nnxnny])/taunsrL+
        (cnsr[(nnx-1)+(j+1)*nnx+k*nnxnny]+cnsr[(nnx-1)+(j-1)*nnx+k*nnxnny]-2*cnsr[(nnx-1)+j*nnx+k*nnxnny])/taunsrT+
        (cnsr[(nnx-1)+j*nnx+(k+1)*nnxnny]+cnsr[(nnx-1)+j*nnx+(k-1)*nnxnny]-2*cnsr[(nnx-1)+j*nnx+k*nnxnny])/taunsrT;
    }
  }
  //y fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=1;i<(nnx-1);i++)
  {
    Icnsr[i+0*nnx+0*nnxnny]=(cnsr[i+0*nnx+0*nnxnny+1]+cnsr[i+0*nnx+0*nnxnny-1]-2*cnsr[i+0*nnx+0*nnxnny])/taunsrL+
      (cnsr[i+(0+1)*nnx+0*nnxnny]+cnsr[i+(0+1)*nnx+0*nnxnny]-2*cnsr[i+0*nnx+0*nnxnny])/taunsrT+
      (cnsr[i+0*nnx+(0+1)*nnxnny]+cnsr[i+0*nnx+(0+1)*nnxnny]-2*cnsr[i+0*nnx+0*nnxnny])/taunsrT;
    Icnsr[i+(nny-1)*nnx+0*nnxnny]=(cnsr[i+(nny-1)*nnx+0*nnxnny+1]+cnsr[i+(nny-1)*nnx+0*nnxnny-1]-2*cnsr[i+(nny-1)*nnx+0*nnxnny])/taunsrL+
      (cnsr[i+((nny-1)-1)*nnx+0*nnxnny]+cnsr[i+((nny-1)-1)*nnx+0*nnxnny]-2*cnsr[i+(nny-1)*nnx+0*nnxnny])/taunsrT+
      (cnsr[i+(nny-1)*nnx+(0+1)*nnxnny]+cnsr[i+(nny-1)*nnx+(0+1)*nnxnny]-2*cnsr[i+(nny-1)*nnx+0*nnxnny])/taunsrT;
    Icnsr[i+0*nnx+(nnz-1)*nnxnny]=(cnsr[i+0*nnx+(nnz-1)*nnxnny+1]+cnsr[i+0*nnx+(nnz-1)*nnxnny-1]-2*cnsr[i+0*nnx+(nnz-1)*nnxnny])/taunsrL+
      (cnsr[i+(0+1)*nnx+(nnz-1)*nnxnny]+cnsr[i+(0+1)*nnx+(nnz-1)*nnxnny]-2*cnsr[i+0*nnx+(nnz-1)*nnxnny])/taunsrT+
      (cnsr[i+0*nnx+((nnz-1)-1)*nnxnny]+cnsr[i+0*nnx+((nnz-1)-1)*nnxnny]-2*cnsr[i+0*nnx+(nnz-1)*nnxnny])/taunsrT;
    Icnsr[i+(nny-1)*nnx+(nnz-1)*nnxnny]=(cnsr[i+(nny-1)*nnx+(nnz-1)*nnxnny+1]+cnsr[i+(nny-1)*nnx+(nnz-1)*nnxnny-1]-2*cnsr[i+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrL+
      (cnsr[i+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cnsr[i+((nny-1)-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[i+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrT+
      (cnsr[i+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cnsr[i+(nny-1)*nnx+((nnz-1)-1)*nnxnny]-2*cnsr[i+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      Icnsr[i+0*nnx+k*nnxnny]=(cnsr[i+0*nnx+k*nnxnny+1]+cnsr[i+0*nnx+k*nnxnny-1]-2*cnsr[i+0*nnx+k*nnxnny])/taunsrL+
        (cnsr[i+(0+1)*nnx+k*nnxnny]+cnsr[i+(0+1)*nnx+k*nnxnny]-2*cnsr[i+0*nnx+k*nnxnny])/taunsrT+
        (cnsr[i+0*nnx+(k+1)*nnxnny]+cnsr[i+0*nnx+(k-1)*nnxnny]-2*cnsr[i+0*nnx+k*nnxnny])/taunsrT;
      Icnsr[i+(nny-1)*nnx+k*nnxnny]=(cnsr[i+(nny-1)*nnx+k*nnxnny+1]+cnsr[i+(nny-1)*nnx+k*nnxnny-1]-2*cnsr[i+(nny-1)*nnx+k*nnxnny])/taunsrL+
        (cnsr[i+((nny-1)-1)*nnx+k*nnxnny]+cnsr[i+((nny-1)-1)*nnx+k*nnxnny]-2*cnsr[i+(nny-1)*nnx+k*nnxnny])/taunsrT+
        (cnsr[i+(nny-1)*nnx+(k+1)*nnxnny]+cnsr[i+(nny-1)*nnx+(k-1)*nnxnny]-2*cnsr[i+(nny-1)*nnx+k*nnxnny])/taunsrT;
    }
    //z fixed
#pragma ivdep
#pragma vector always
    for (int j=1;j<nny-1;j++)
    {
      Icnsr[i+j*nnx+0*nnxnny]=(cnsr[i+j*nnx+0*nnxnny+1]+cnsr[i+j*nnx+0*nnxnny-1]-2*cnsr[i+j*nnx+0*nnxnny])/taunsrL+
        (cnsr[i+(j+1)*nnx+0*nnxnny]+cnsr[i+(j-1)*nnx+0*nnxnny]-2*cnsr[i+j*nnx+0*nnxnny])/taunsrT+
        (cnsr[i+j*nnx+(0+1)*nnxnny]+cnsr[i+j*nnx+(0+1)*nnxnny]-2*cnsr[i+j*nnx+0*nnxnny])/taunsrT;
      Icnsr[i+j*nnx+(nnz-1)*nnxnny]=(cnsr[i+j*nnx+(nnz-1)*nnxnny+1]+cnsr[i+j*nnx+(nnz-1)*nnxnny-1]-2*cnsr[i+j*nnx+(nnz-1)*nnxnny])/taunsrL+
        (cnsr[i+(j+1)*nnx+(nnz-1)*nnxnny]+cnsr[i+(j-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[i+j*nnx+(nnz-1)*nnxnny])/taunsrT+
        (cnsr[i+j*nnx+((nnz-1)-1)*nnxnny]+cnsr[i+j*nnx+((nnz-1)-1)*nnxnny]-2*cnsr[i+j*nnx+(nnz-1)*nnxnny])/taunsrT;
    }
  }
#pragma ivdep
#pragma vector always
  for (int k=1;k<nnz-1;k++)
  {
    Icnsr[0+0*nnx+k*nnxnny]=(cnsr[0+0*nnx+k*nnxnny+1]+cnsr[0+0*nnx+k*nnxnny+1]-2*cnsr[0+0*nnx+k*nnxnny])/taunsrL+
      (cnsr[0+(0+1)*nnx+k*nnxnny]+cnsr[0+(0+1)*nnx+k*nnxnny]-2*cnsr[0+0*nnx+k*nnxnny])/taunsrT+
      (cnsr[0+0*nnx+(k+1)*nnxnny]+cnsr[0+0*nnx+(k-1)*nnxnny]-2*cnsr[0+0*nnx+k*nnxnny])/taunsrT;
    Icnsr[0+(nny-1)*nnx+k*nnxnny]=(cnsr[0+(nny-1)*nnx+k*nnxnny+1]+cnsr[0+(nny-1)*nnx+k*nnxnny+1]-2*cnsr[0+(nny-1)*nnx+k*nnxnny])/taunsrL+
      (cnsr[0+((nny-1)-1)*nnx+k*nnxnny]+cnsr[0+((nny-1)-1)*nnx+k*nnxnny]-2*cnsr[0+(nny-1)*nnx+k*nnxnny])/taunsrT+
      (cnsr[0+(nny-1)*nnx+(k+1)*nnxnny]+cnsr[0+(nny-1)*nnx+(k-1)*nnxnny]-2*cnsr[0+(nny-1)*nnx+k*nnxnny])/taunsrT;
    Icnsr[(nnx-1)+0*nnx+k*nnxnny]=(cnsr[(nnx-1)+0*nnx+k*nnxnny-1]+cnsr[(nnx-1)+0*nnx+k*nnxnny-1]-2*cnsr[(nnx-1)+0*nnx+k*nnxnny])/taunsrL+
      (cnsr[(nnx-1)+(0+1)*nnx+k*nnxnny]+cnsr[(nnx-1)+(0+1)*nnx+k*nnxnny]-2*cnsr[(nnx-1)+0*nnx+k*nnxnny])/taunsrT+
      (cnsr[(nnx-1)+0*nnx+(k+1)*nnxnny]+cnsr[(nnx-1)+0*nnx+(k-1)*nnxnny]-2*cnsr[(nnx-1)+0*nnx+k*nnxnny])/taunsrT;
    Icnsr[(nnx-1)+(nny-1)*nnx+k*nnxnny]=(cnsr[(nnx-1)+(nny-1)*nnx+k*nnxnny-1]+cnsr[(nnx-1)+(nny-1)*nnx+k*nnxnny-1]-2*cnsr[(nnx-1)+(nny-1)*nnx+k*nnxnny])/taunsrL+
      (cnsr[(nnx-1)+((nny-1)-1)*nnx+k*nnxnny]+cnsr[(nnx-1)+((nny-1)-1)*nnx+k*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+k*nnxnny])/taunsrT+
      (cnsr[(nnx-1)+(nny-1)*nnx+(k+1)*nnxnny]+cnsr[(nnx-1)+(nny-1)*nnx+(k-1)*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+k*nnxnny])/taunsrT;
  }
#endif
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int k=1;k<nnz-1;k++)
  {
    for (int j=1;j<nny-1;j++)
    {
#pragma ivdep
#pragma vector always
      for (int i=1;i<nnx-1;i++)
      {
        Icnsr[i+j*nnx+k*nnxnny]=(cnsr[i+j*nnx+k*nnxnny+1]+cnsr[i+j*nnx+k*nnxnny-1]-2*cnsr[i+j*nnx+k*nnxnny])/taunsrL+
          (cnsr[i+(j+1)*nnx+k*nnxnny]+cnsr[i+(j-1)*nnx+k*nnxnny]-2*cnsr[i+j*nnx+k*nnxnny])/taunsrT+
          (cnsr[i+j*nnx+(k+1)*nnxnny]+cnsr[i+j*nnx+(k-1)*nnxnny]-2*cnsr[i+j*nnx+k*nnxnny])/taunsrT;
      }
    }
  }
}


#ifdef ___NO_CS_BUFFER
void CSubcell::computecsmn(void)
{
#ifdef ___PERIODIC
  csmn[0+0*nnx+0*nnxnny]=(cs[0+0*nnx+0*nnxnny+1]+cs[(nnx-1)+0*nnx+0*nnxnny])/tausL+
    (cs[0+(0+1)*nnx+0*nnxnny]+cs[0+(nny-1)*nnx+0*nnxnny])/tausT+
    (cs[0+0*nnx+(0+1)*nnxnny]+cs[0+0*nnx+(nnz-1)*nnxnny])/tausT;
  csmn[0+(nny-1)*nnx+(nnz-1)*nnxnny]=(cs[0+(nny-1)*nnx+(nnz-1)*nnxnny+1]+cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tausL+
    (cs[0+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cs[0+(0)*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[0+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cs[0+(nny-1)*nnx+(0)*nnxnny])/tausT;
  csmn[0+(nny-1)*nnx+0*nnxnny]=(cs[0+(nny-1)*nnx+0*nnxnny+1]+cs[(nnx-1)+(nny-1)*nnx+0*nnxnny])/tausL+
    (cs[0+((nny-1)-1)*nnx+0*nnxnny]+cs[0+(0)*nnx+0*nnxnny])/tausT+
    (cs[0+(nny-1)*nnx+(0+1)*nnxnny]+cs[0+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT;
  csmn[0+0*nnx+(nnz-1)*nnxnny]=(cs[0+0*nnx+(nnz-1)*nnxnny+1]+cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/tausL+
    (cs[0+(0+1)*nnx+(nnz-1)*nnxnny]+cs[0+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[0+0*nnx+((nnz-1)-1)*nnxnny]+cs[0+0*nnx+(0)*nnxnny])/tausT;
  csmn[(nnx-1)+(nny-1)*nnx+0*nnxnny]=(cs[(nnx-1)+(nny-1)*nnx+0*nnxnny-1]+cs[(0)+(nny-1)*nnx+0*nnxnny])/tausL+
    (cs[(nnx-1)+((nny-1)-1)*nnx+0*nnxnny]+cs[(nnx-1)+(0)*nnx+0*nnxnny])/tausT+
    (cs[(nnx-1)+(nny-1)*nnx+(0+1)*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT;
  csmn[(nnx-1)+0*nnx+(nnz-1)*nnxnny]=(cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny-1]+cs[(0)+0*nnx+(nnz-1)*nnxnny])/tausL+
    (cs[(nnx-1)+(0+1)*nnx+(nnz-1)*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[(nnx-1)+0*nnx+((nnz-1)-1)*nnxnny]+cs[(nnx-1)+0*nnx+(0)*nnxnny])/tausT;
  csmn[(nnx-1)+0*nnx+0*nnxnny]=(cs[(nnx-1)+0*nnx+0*nnxnny-1]+cs[(0)+0*nnx+0*nnxnny])/tausL+
    (cs[(nnx-1)+(0+1)*nnx+0*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+0*nnxnny])/tausT+
    (cs[(nnx-1)+0*nnx+(0+1)*nnxnny]+cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/tausT;
  csmn[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]=(cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny-1]+cs[(0)+(nny-1)*nnx+(nnz-1)*nnxnny])/tausL+
    (cs[(nnx-1)+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cs[(nnx-1)+(0)*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[(nnx-1)+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+(0)*nnxnny])/tausT;
  //x fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int j=1;j<nny-1;j++)
  {
    csmn[0+j*nnx+0*nnxnny]=(cs[0+j*nnx+0*nnxnny+1]+cs[(nnx-1)+j*nnx+0*nnxnny])/tausL+
      (cs[0+(j+1)*nnx+0*nnxnny]+cs[0+(j-1)*nnx+0*nnxnny])/tausT+
      (cs[0+j*nnx+(0+1)*nnxnny]+cs[0+j*nnx+(nnz-1)*nnxnny])/tausT;
    csmn[(nnx-1)+j*nnx+0*nnxnny]=(cs[(nnx-1)+j*nnx+0*nnxnny-1]+cs[(0)+j*nnx+0*nnxnny])/tausL+
      (cs[(nnx-1)+(j+1)*nnx+0*nnxnny]+cs[(nnx-1)+(j-1)*nnx+0*nnxnny])/tausT+
      (cs[(nnx-1)+j*nnx+(0+1)*nnxnny]+cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/tausT;
    csmn[0+j*nnx+(nnz-1)*nnxnny]=(cs[0+j*nnx+(nnz-1)*nnxnny+1]+cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/tausL+
      (cs[0+(j+1)*nnx+(nnz-1)*nnxnny]+cs[0+(j-1)*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[0+j*nnx+((nnz-1)-1)*nnxnny]+cs[0+j*nnx+(0)*nnxnny])/tausT;
    csmn[(nnx-1)+j*nnx+(nnz-1)*nnxnny]=(cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny-1]+cs[(0)+j*nnx+(nnz-1)*nnxnny])/tausL+
      (cs[(nnx-1)+(j+1)*nnx+(nnz-1)*nnxnny]+cs[(nnx-1)+(j-1)*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[(nnx-1)+j*nnx+((nnz-1)-1)*nnxnny]+cs[(nnx-1)+j*nnx+(0)*nnxnny])/tausT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      csmn[0+j*nnx+k*nnxnny]=(cs[0+j*nnx+k*nnxnny+1]+cs[(nnx-1)+j*nnx+k*nnxnny])/tausL+
        (cs[0+(j+1)*nnx+k*nnxnny]+cs[0+(j-1)*nnx+k*nnxnny])/tausT+
        (cs[0+j*nnx+(k+1)*nnxnny]+cs[0+j*nnx+(k-1)*nnxnny])/tausT;
      csmn[(nnx-1)+j*nnx+k*nnxnny]=(cs[(nnx-1)+j*nnx+k*nnxnny-1]+cs[(0)+j*nnx+k*nnxnny])/tausL+
        (cs[(nnx-1)+(j+1)*nnx+k*nnxnny]+cs[(nnx-1)+(j-1)*nnx+k*nnxnny])/tausT+
        (cs[(nnx-1)+j*nnx+(k+1)*nnxnny]+cs[(nnx-1)+j*nnx+(k-1)*nnxnny])/tausT;
    }
  }
  //y fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=1;i<(nnx-1);i++)
  {
    csmn[i+0*nnx+0*nnxnny]=(cs[i+0*nnx+0*nnxnny+1]+cs[i+0*nnx+0*nnxnny-1])/tausL+
      (cs[i+(0+1)*nnx+0*nnxnny]+cs[i+(nny-1)*nnx+0*nnxnny])/tausT+
      (cs[i+0*nnx+(0+1)*nnxnny]+cs[i+0*nnx+(nnz-1)*nnxnny])/tausT;
    csmn[i+(nny-1)*nnx+0*nnxnny]=(cs[i+(nny-1)*nnx+0*nnxnny+1]+cs[i+(nny-1)*nnx+0*nnxnny-1])/tausL+
      (cs[i+((nny-1)-1)*nnx+0*nnxnny]+cs[i+(0)*nnx+0*nnxnny])/tausT+
      (cs[i+(nny-1)*nnx+(0+1)*nnxnny]+cs[i+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT;
    csmn[i+0*nnx+(nnz-1)*nnxnny]=(cs[i+0*nnx+(nnz-1)*nnxnny+1]+cs[i+0*nnx+(nnz-1)*nnxnny-1])/tausL+
      (cs[i+(0+1)*nnx+(nnz-1)*nnxnny]+cs[i+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[i+0*nnx+((nnz-1)-1)*nnxnny]+cs[i+0*nnx+(0)*nnxnny])/tausT;
    csmn[i+(nny-1)*nnx+(nnz-1)*nnxnny]=(cs[i+(nny-1)*nnx+(nnz-1)*nnxnny+1]+cs[i+(nny-1)*nnx+(nnz-1)*nnxnny-1])/tausL+
      (cs[i+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cs[i+(0)*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[i+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cs[i+(nny-1)*nnx+(0)*nnxnny])/tausT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      csmn[i+0*nnx+k*nnxnny]=(cs[i+0*nnx+k*nnxnny+1]+cs[i+0*nnx+k*nnxnny-1])/tausL+
        (cs[i+(0+1)*nnx+k*nnxnny]+cs[i+(nny-1)*nnx+k*nnxnny])/tausT+
        (cs[i+0*nnx+(k+1)*nnxnny]+cs[i+0*nnx+(k-1)*nnxnny])/tausT;
      csmn[i+(nny-1)*nnx+k*nnxnny]=(cs[i+(nny-1)*nnx+k*nnxnny+1]+cs[i+(nny-1)*nnx+k*nnxnny-1])/tausL+
        (cs[i+((nny-1)-1)*nnx+k*nnxnny]+cs[i+(0)*nnx+k*nnxnny])/tausT+
        (cs[i+(nny-1)*nnx+(k+1)*nnxnny]+cs[i+(nny-1)*nnx+(k-1)*nnxnny])/tausT;
    }
    //z fixed
#pragma ivdep
#pragma vector always
    for (int j=1;j<nny-1;j++)
    {
      csmn[i+j*nnx+0*nnxnny]=(cs[i+j*nnx+0*nnxnny+1]+cs[i+j*nnx+0*nnxnny-1])/tausL+
        (cs[i+(j+1)*nnx+0*nnxnny]+cs[i+(j-1)*nnx+0*nnxnny])/tausT+
        (cs[i+j*nnx+(0+1)*nnxnny]+cs[i+j*nnx+(nnz-1)*nnxnny])/tausT;
      csmn[i+j*nnx+(nnz-1)*nnxnny]=(cs[i+j*nnx+(nnz-1)*nnxnny+1]+cs[i+j*nnx+(nnz-1)*nnxnny-1])/tausL+
        (cs[i+(j+1)*nnx+(nnz-1)*nnxnny]+cs[i+(j-1)*nnx+(nnz-1)*nnxnny])/tausT+
        (cs[i+j*nnx+((nnz-1)-1)*nnxnny]+cs[i+j*nnx+(0)*nnxnny])/tausT;
    }
  }
#pragma ivdep
#pragma vector always
  for (int k=1;k<nnz-1;k++)
  {
    csmn[0+0*nnx+k*nnxnny]=(cs[0+0*nnx+k*nnxnny+1]+cs[(nnx-1)+0*nnx+k*nnxnny])/tausL+
      (cs[0+(0+1)*nnx+k*nnxnny]+cs[0+(nny-1)*nnx+k*nnxnny])/tausT+
      (cs[0+0*nnx+(k+1)*nnxnny]+cs[0+0*nnx+(k-1)*nnxnny])/tausT;
    csmn[0+(nny-1)*nnx+k*nnxnny]=(cs[0+(nny-1)*nnx+k*nnxnny+1]+cs[(nnx-1)+(nny-1)*nnx+k*nnxnny])/tausL+
      (cs[0+((nny-1)-1)*nnx+k*nnxnny]+cs[0+(0)*nnx+k*nnxnny])/tausT+
      (cs[0+(nny-1)*nnx+(k+1)*nnxnny]+cs[0+(nny-1)*nnx+(k-1)*nnxnny])/tausT;
    csmn[(nnx-1)+0*nnx+k*nnxnny]=(cs[(nnx-1)+0*nnx+k*nnxnny-1]+cs[(0)+0*nnx+k*nnxnny])/tausL+
      (cs[(nnx-1)+(0+1)*nnx+k*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+k*nnxnny])/tausT+
      (cs[(nnx-1)+0*nnx+(k+1)*nnxnny]+cs[(nnx-1)+0*nnx+(k-1)*nnxnny])/tausT;
    csmn[(nnx-1)+(nny-1)*nnx+k*nnxnny]=(cs[(nnx-1)+(nny-1)*nnx+k*nnxnny-1]+cs[0+(nny-1)*nnx+k*nnxnny])/tausL+
      (cs[(nnx-1)+((nny-1)-1)*nnx+k*nnxnny]+cs[(nnx-1)+(0)*nnx+k*nnxnny])/tausT+
      (cs[(nnx-1)+(nny-1)*nnx+(k+1)*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+(k-1)*nnxnny])/tausT;
  }
#else
  csmn[0+0*nnx+0*nnxnny]=(cs[0+0*nnx+0*nnxnny+1]+cs[0+0*nnx+0*nnxnny+1])/tausL+
    (cs[0+(0+1)*nnx+0*nnxnny]+cs[0+(0+1)*nnx+0*nnxnny])/tausT+
    (cs[0+0*nnx+(0+1)*nnxnny]+cs[0+0*nnx+(0+1)*nnxnny])/tausT;
  csmn[0+(nny-1)*nnx+(nnz-1)*nnxnny]=(cs[0+(nny-1)*nnx+(nnz-1)*nnxnny+1]+cs[0+(nny-1)*nnx+(nnz-1)*nnxnny+1])/tausL+
    (cs[0+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cs[0+((nny-1)-1)*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[0+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cs[0+(nny-1)*nnx+((nnz-1)-1)*nnxnny])/tausT;
  csmn[0+(nny-1)*nnx+0*nnxnny]=(cs[0+(nny-1)*nnx+0*nnxnny+1]+cs[0+(nny-1)*nnx+0*nnxnny+1])/tausL+
    (cs[0+((nny-1)-1)*nnx+0*nnxnny]+cs[0+((nny-1)-1)*nnx+0*nnxnny])/tausT+
    (cs[0+(nny-1)*nnx+(0+1)*nnxnny]+cs[0+(nny-1)*nnx+(0+1)*nnxnny])/tausT;
  csmn[0+0*nnx+(nnz-1)*nnxnny]=(cs[0+0*nnx+(nnz-1)*nnxnny+1]+cs[0+0*nnx+(nnz-1)*nnxnny+1])/tausL+
    (cs[0+(0+1)*nnx+(nnz-1)*nnxnny]+cs[0+(0+1)*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[0+0*nnx+((nnz-1)-1)*nnxnny]+cs[0+0*nnx+((nnz-1)-1)*nnxnny])/tausT;
  csmn[(nnx-1)+(nny-1)*nnx+0*nnxnny]=(cs[(nnx-1)+(nny-1)*nnx+0*nnxnny-1]+cs[(nnx-1)+(nny-1)*nnx+0*nnxnny-1])/tausL+
    (cs[(nnx-1)+((nny-1)-1)*nnx+0*nnxnny]+cs[(nnx-1)+((nny-1)-1)*nnx+0*nnxnny])/tausT+
    (cs[(nnx-1)+(nny-1)*nnx+(0+1)*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+(0+1)*nnxnny])/tausT;
  csmn[(nnx-1)+0*nnx+(nnz-1)*nnxnny]=(cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny-1]+cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny-1])/tausL+
    (cs[(nnx-1)+(0+1)*nnx+(nnz-1)*nnxnny]+cs[(nnx-1)+(0+1)*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[(nnx-1)+0*nnx+((nnz-1)-1)*nnxnny]+cs[(nnx-1)+0*nnx+((nnz-1)-1)*nnxnny])/tausT;
  csmn[(nnx-1)+0*nnx+0*nnxnny]=(cs[(nnx-1)+0*nnx+0*nnxnny-1]+cs[(nnx-1)+0*nnx+0*nnxnny-1])/tausL+
    (cs[(nnx-1)+(0+1)*nnx+0*nnxnny]+cs[(nnx-1)+(0+1)*nnx+0*nnxnny])/tausT+
    (cs[(nnx-1)+0*nnx+(0+1)*nnxnny]+cs[(nnx-1)+0*nnx+(0+1)*nnxnny])/tausT;
  csmn[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]=(cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny-1]+cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny-1])/tausL+
    (cs[(nnx-1)+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cs[(nnx-1)+((nny-1)-1)*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[(nnx-1)+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+((nnz-1)-1)*nnxnny])/tausT;
  //x fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int j=1;j<nny-1;j++)
  {
    csmn[0+j*nnx+0*nnxnny]=(cs[0+j*nnx+0*nnxnny+1]+cs[0+j*nnx+0*nnxnny+1])/tausL+
      (cs[0+(j+1)*nnx+0*nnxnny]+cs[0+(j-1)*nnx+0*nnxnny])/tausT+
      (cs[0+j*nnx+(0+1)*nnxnny]+cs[0+j*nnx+(0+1)*nnxnny])/tausT;
    csmn[(nnx-1)+j*nnx+0*nnxnny]=(cs[(nnx-1)+j*nnx+0*nnxnny-1]+cs[(nnx-1)+j*nnx+0*nnxnny-1])/tausL+
      (cs[(nnx-1)+(j+1)*nnx+0*nnxnny]+cs[(nnx-1)+(j-1)*nnx+0*nnxnny])/tausT+
      (cs[(nnx-1)+j*nnx+(0+1)*nnxnny]+cs[(nnx-1)+j*nnx+(0+1)*nnxnny])/tausT;
    csmn[0+j*nnx+(nnz-1)*nnxnny]=(cs[0+j*nnx+(nnz-1)*nnxnny+1]+cs[0+j*nnx+(nnz-1)*nnxnny+1])/tausL+
      (cs[0+(j+1)*nnx+(nnz-1)*nnxnny]+cs[0+(j-1)*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[0+j*nnx+((nnz-1)-1)*nnxnny]+cs[0+j*nnx+((nnz-1)-1)*nnxnny])/tausT;
    csmn[(nnx-1)+j*nnx+(nnz-1)*nnxnny]=(cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny-1]+cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny-1])/tausL+
      (cs[(nnx-1)+(j+1)*nnx+(nnz-1)*nnxnny]+cs[(nnx-1)+(j-1)*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[(nnx-1)+j*nnx+((nnz-1)-1)*nnxnny]+cs[(nnx-1)+j*nnx+((nnz-1)-1)*nnxnny])/tausT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      csmn[0+j*nnx+k*nnxnny]=(cs[0+j*nnx+k*nnxnny+1]+cs[0+j*nnx+k*nnxnny+1])/tausL+
        (cs[0+(j+1)*nnx+k*nnxnny]+cs[0+(j-1)*nnx+k*nnxnny])/tausT+
        (cs[0+j*nnx+(k+1)*nnxnny]+cs[0+j*nnx+(k-1)*nnxnny])/tausT;
      csmn[(nnx-1)+j*nnx+k*nnxnny]=(cs[(nnx-1)+j*nnx+k*nnxnny-1]+cs[(nnx-1)+j*nnx+k*nnxnny-1])/tausL+
        (cs[(nnx-1)+(j+1)*nnx+k*nnxnny]+cs[(nnx-1)+(j-1)*nnx+k*nnxnny])/tausT+
        (cs[(nnx-1)+j*nnx+(k+1)*nnxnny]+cs[(nnx-1)+j*nnx+(k-1)*nnxnny])/tausT;
    }
  }
  //y fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=1;i<(nnx-1);i++)
  {
    csmn[i+0*nnx+0*nnxnny]=(cs[i+0*nnx+0*nnxnny+1]+cs[i+0*nnx+0*nnxnny-1])/tausL+
      (cs[i+(0+1)*nnx+0*nnxnny]+cs[i+(0+1)*nnx+0*nnxnny])/tausT+
      (cs[i+0*nnx+(0+1)*nnxnny]+cs[i+0*nnx+(0+1)*nnxnny])/tausT;
    csmn[i+(nny-1)*nnx+0*nnxnny]=(cs[i+(nny-1)*nnx+0*nnxnny+1]+cs[i+(nny-1)*nnx+0*nnxnny-1])/tausL+
      (cs[i+((nny-1)-1)*nnx+0*nnxnny]+cs[i+((nny-1)-1)*nnx+0*nnxnny])/tausT+
      (cs[i+(nny-1)*nnx+(0+1)*nnxnny]+cs[i+(nny-1)*nnx+(0+1)*nnxnny])/tausT;
    csmn[i+0*nnx+(nnz-1)*nnxnny]=(cs[i+0*nnx+(nnz-1)*nnxnny+1]+cs[i+0*nnx+(nnz-1)*nnxnny-1])/tausL+
      (cs[i+(0+1)*nnx+(nnz-1)*nnxnny]+cs[i+(0+1)*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[i+0*nnx+((nnz-1)-1)*nnxnny]+cs[i+0*nnx+((nnz-1)-1)*nnxnny])/tausT;
    csmn[i+(nny-1)*nnx+(nnz-1)*nnxnny]=(cs[i+(nny-1)*nnx+(nnz-1)*nnxnny+1]+cs[i+(nny-1)*nnx+(nnz-1)*nnxnny-1])/tausL+
      (cs[i+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cs[i+((nny-1)-1)*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[i+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cs[i+(nny-1)*nnx+((nnz-1)-1)*nnxnny])/tausT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      csmn[i+0*nnx+k*nnxnny]=(cs[i+0*nnx+k*nnxnny+1]+cs[i+0*nnx+k*nnxnny-1])/tausL+
        (cs[i+(0+1)*nnx+k*nnxnny]+cs[i+(0+1)*nnx+k*nnxnny])/tausT+
        (cs[i+0*nnx+(k+1)*nnxnny]+cs[i+0*nnx+(k-1)*nnxnny])/tausT;
      csmn[i+(nny-1)*nnx+k*nnxnny]=(cs[i+(nny-1)*nnx+k*nnxnny+1]+cs[i+(nny-1)*nnx+k*nnxnny-1])/tausL+
        (cs[i+((nny-1)-1)*nnx+k*nnxnny]+cs[i+((nny-1)-1)*nnx+k*nnxnny])/tausT+
        (cs[i+(nny-1)*nnx+(k+1)*nnxnny]+cs[i+(nny-1)*nnx+(k-1)*nnxnny])/tausT;
    }
    //z fixed
#pragma ivdep
#pragma vector always
    for (int j=1;j<nny-1;j++)
    {
      csmn[i+j*nnx+0*nnxnny]=(cs[i+j*nnx+0*nnxnny+1]+cs[i+j*nnx+0*nnxnny-1])/tausL+
        (cs[i+(j+1)*nnx+0*nnxnny]+cs[i+(j-1)*nnx+0*nnxnny])/tausT+
        (cs[i+j*nnx+(0+1)*nnxnny]+cs[i+j*nnx+(0+1)*nnxnny])/tausT;
      csmn[i+j*nnx+(nnz-1)*nnxnny]=(cs[i+j*nnx+(nnz-1)*nnxnny+1]+cs[i+j*nnx+(nnz-1)*nnxnny-1])/tausL+
        (cs[i+(j+1)*nnx+(nnz-1)*nnxnny]+cs[i+(j-1)*nnx+(nnz-1)*nnxnny])/tausT+
        (cs[i+j*nnx+((nnz-1)-1)*nnxnny]+cs[i+j*nnx+((nnz-1)-1)*nnxnny])/tausT;
    }
  }
#pragma ivdep
#pragma vector always
  for (int k=1;k<nnz-1;k++)
  {
    csmn[0+0*nnx+k*nnxnny]=(cs[0+0*nnx+k*nnxnny+1]+cs[0+0*nnx+k*nnxnny+1])/tausL+
      (cs[0+(0+1)*nnx+k*nnxnny]+cs[0+(0+1)*nnx+k*nnxnny])/tausT+
      (cs[0+0*nnx+(k+1)*nnxnny]+cs[0+0*nnx+(k-1)*nnxnny])/tausT;
    csmn[0+(nny-1)*nnx+k*nnxnny]=(cs[0+(nny-1)*nnx+k*nnxnny+1]+cs[0+(nny-1)*nnx+k*nnxnny+1])/tausL+
      (cs[0+((nny-1)-1)*nnx+k*nnxnny]+cs[0+((nny-1)-1)*nnx+k*nnxnny])/tausT+
      (cs[0+(nny-1)*nnx+(k+1)*nnxnny]+cs[0+(nny-1)*nnx+(k-1)*nnxnny])/tausT;
    csmn[(nnx-1)+0*nnx+k*nnxnny]=(cs[(nnx-1)+0*nnx+k*nnxnny-1]+cs[(nnx-1)+0*nnx+k*nnxnny-1])/tausL+
      (cs[(nnx-1)+(0+1)*nnx+k*nnxnny]+cs[(nnx-1)+(0+1)*nnx+k*nnxnny])/tausT+
      (cs[(nnx-1)+0*nnx+(k+1)*nnxnny]+cs[(nnx-1)+0*nnx+(k-1)*nnxnny])/tausT;
    csmn[(nnx-1)+(nny-1)*nnx+k*nnxnny]=(cs[(nnx-1)+(nny-1)*nnx+k*nnxnny-1]+cs[(nnx-1)+(nny-1)*nnx+k*nnxnny-1])/tausL+
      (cs[(nnx-1)+((nny-1)-1)*nnx+k*nnxnny]+cs[(nnx-1)+((nny-1)-1)*nnx+k*nnxnny])/tausT+
      (cs[(nnx-1)+(nny-1)*nnx+(k+1)*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+(k-1)*nnxnny])/tausT;
  }
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int k=1;k<nnz-1;k++)
  {
    for (int j=1;j<nny-1;j++)
    {
#pragma ivdep
#pragma vector always
      for (int i=1;i<nnx-1;i++)
      {
        csmn[i+j*nnx+k*nnxnny]=(cs[i+j*nnx+k*nnxnny+1]+cs[i+j*nnx+k*nnxnny-1])/tausL+
          (cs[i+(j+1)*nnx+k*nnxnny]+cs[i+(j-1)*nnx+k*nnxnny])/tausT+
          (cs[i+j*nnx+(k+1)*nnxnny]+cs[i+j*nnx+(k-1)*nnxnny])/tausT;
      }
    }
  }
}
#else
void CSubcell::computeIcs(void)
{

#ifdef ___PERIODIC
  Ics[0+0*nnx+0*nnxnny]=(cs[0+0*nnx+0*nnxnny+1]+cs[(nnx-1)+0*nnx+0*nnxnny]-2*cs[0+0*nnx+0*nnxnny])/tausL+
    (cs[0+(0+1)*nnx+0*nnxnny]+cs[0+(nny-1)*nnx+0*nnxnny]-2*cs[0+0*nnx+0*nnxnny])/tausT+
    (cs[0+0*nnx+(0+1)*nnxnny]+cs[0+0*nnx+(nnz-1)*nnxnny]-2*cs[0+0*nnx+0*nnxnny])/tausT;
  Ics[0+(nny-1)*nnx+(nnz-1)*nnxnny]=(cs[0+(nny-1)*nnx+(nnz-1)*nnxnny+1]+cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cs[0+(nny-1)*nnx+(nnz-1)*nnxnny])/tausL+
    (cs[0+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cs[0+(0)*nnx+(nnz-1)*nnxnny]-2*cs[0+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[0+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cs[0+(nny-1)*nnx+(0)*nnxnny]-2*cs[0+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT;
  Ics[0+(nny-1)*nnx+0*nnxnny]=(cs[0+(nny-1)*nnx+0*nnxnny+1]+cs[(nnx-1)+(nny-1)*nnx+0*nnxnny]-2*cs[0+(nny-1)*nnx+0*nnxnny])/tausL+
    (cs[0+((nny-1)-1)*nnx+0*nnxnny]+cs[0+(0)*nnx+0*nnxnny]-2*cs[0+(nny-1)*nnx+0*nnxnny])/tausT+
    (cs[0+(nny-1)*nnx+(0+1)*nnxnny]+cs[0+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cs[0+(nny-1)*nnx+0*nnxnny])/tausT;
  Ics[0+0*nnx+(nnz-1)*nnxnny]=(cs[0+0*nnx+(nnz-1)*nnxnny+1]+cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny]-2*cs[0+0*nnx+(nnz-1)*nnxnny])/tausL+
    (cs[0+(0+1)*nnx+(nnz-1)*nnxnny]+cs[0+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cs[0+0*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[0+0*nnx+((nnz-1)-1)*nnxnny]+cs[0+0*nnx+(0)*nnxnny]-2*cs[0+0*nnx+(nnz-1)*nnxnny])/tausT;
  Ics[(nnx-1)+(nny-1)*nnx+0*nnxnny]=(cs[(nnx-1)+(nny-1)*nnx+0*nnxnny-1]+cs[(0)+(nny-1)*nnx+0*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+0*nnxnny])/tausL+
    (cs[(nnx-1)+((nny-1)-1)*nnx+0*nnxnny]+cs[(nnx-1)+(0)*nnx+0*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+0*nnxnny])/tausT+
    (cs[(nnx-1)+(nny-1)*nnx+(0+1)*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+0*nnxnny])/tausT;
  Ics[(nnx-1)+0*nnx+(nnz-1)*nnxnny]=(cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny-1]+cs[(0)+0*nnx+(nnz-1)*nnxnny]-2*cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/tausL+
    (cs[(nnx-1)+(0+1)*nnx+(nnz-1)*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[(nnx-1)+0*nnx+((nnz-1)-1)*nnxnny]+cs[(nnx-1)+0*nnx+(0)*nnxnny]-2*cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/tausT;
  Ics[(nnx-1)+0*nnx+0*nnxnny]=(cs[(nnx-1)+0*nnx+0*nnxnny-1]+cs[(0)+0*nnx+0*nnxnny]-2*cs[(nnx-1)+0*nnx+0*nnxnny])/tausL+
    (cs[(nnx-1)+(0+1)*nnx+0*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+0*nnxnny]-2*cs[(nnx-1)+0*nnx+0*nnxnny])/tausT+
    (cs[(nnx-1)+0*nnx+(0+1)*nnxnny]+cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny]-2*cs[(nnx-1)+0*nnx+0*nnxnny])/tausT;
  Ics[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]=(cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny-1]+cs[(0)+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tausL+
    (cs[(nnx-1)+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cs[(nnx-1)+(0)*nnx+(nnz-1)*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[(nnx-1)+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+(0)*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT;
  //x fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int j=1;j<nny-1;j++)
  {
    Ics[0+j*nnx+0*nnxnny]=(cs[0+j*nnx+0*nnxnny+1]+cs[(nnx-1)+j*nnx+0*nnxnny]-2*cs[0+j*nnx+0*nnxnny])/tausL+
      (cs[0+(j+1)*nnx+0*nnxnny]+cs[0+(j-1)*nnx+0*nnxnny]-2*cs[0+j*nnx+0*nnxnny])/tausT+
      (cs[0+j*nnx+(0+1)*nnxnny]+cs[0+j*nnx+(nnz-1)*nnxnny]-2*cs[0+j*nnx+0*nnxnny])/tausT;
    Ics[(nnx-1)+j*nnx+0*nnxnny]=(cs[(nnx-1)+j*nnx+0*nnxnny-1]+cs[(0)+j*nnx+0*nnxnny]-2*cs[(nnx-1)+j*nnx+0*nnxnny])/tausL+
      (cs[(nnx-1)+(j+1)*nnx+0*nnxnny]+cs[(nnx-1)+(j-1)*nnx+0*nnxnny]-2*cs[(nnx-1)+j*nnx+0*nnxnny])/tausT+
      (cs[(nnx-1)+j*nnx+(0+1)*nnxnny]+cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny]-2*cs[(nnx-1)+j*nnx+0*nnxnny])/tausT;
    Ics[0+j*nnx+(nnz-1)*nnxnny]=(cs[0+j*nnx+(nnz-1)*nnxnny+1]+cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny]-2*cs[0+j*nnx+(nnz-1)*nnxnny])/tausL+
      (cs[0+(j+1)*nnx+(nnz-1)*nnxnny]+cs[0+(j-1)*nnx+(nnz-1)*nnxnny]-2*cs[0+j*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[0+j*nnx+((nnz-1)-1)*nnxnny]+cs[0+j*nnx+(0)*nnxnny]-2*cs[0+j*nnx+(nnz-1)*nnxnny])/tausT;
    Ics[(nnx-1)+j*nnx+(nnz-1)*nnxnny]=(cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny-1]+cs[(0)+j*nnx+(nnz-1)*nnxnny]-2*cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/tausL+
      (cs[(nnx-1)+(j+1)*nnx+(nnz-1)*nnxnny]+cs[(nnx-1)+(j-1)*nnx+(nnz-1)*nnxnny]-2*cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[(nnx-1)+j*nnx+((nnz-1)-1)*nnxnny]+cs[(nnx-1)+j*nnx+(0)*nnxnny]-2*cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/tausT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      Ics[0+j*nnx+k*nnxnny]=(cs[0+j*nnx+k*nnxnny+1]+cs[(nnx-1)+j*nnx+k*nnxnny]-2*cs[0+j*nnx+k*nnxnny])/tausL+
        (cs[0+(j+1)*nnx+k*nnxnny]+cs[0+(j-1)*nnx+k*nnxnny]-2*cs[0+j*nnx+k*nnxnny])/tausT+
        (cs[0+j*nnx+(k+1)*nnxnny]+cs[0+j*nnx+(k-1)*nnxnny]-2*cs[0+j*nnx+k*nnxnny])/tausT;
      Ics[(nnx-1)+j*nnx+k*nnxnny]=(cs[(nnx-1)+j*nnx+k*nnxnny-1]+cs[(0)+j*nnx+k*nnxnny]-2*cs[(nnx-1)+j*nnx+k*nnxnny])/tausL+
        (cs[(nnx-1)+(j+1)*nnx+k*nnxnny]+cs[(nnx-1)+(j-1)*nnx+k*nnxnny]-2*cs[(nnx-1)+j*nnx+k*nnxnny])/tausT+
        (cs[(nnx-1)+j*nnx+(k+1)*nnxnny]+cs[(nnx-1)+j*nnx+(k-1)*nnxnny]-2*cs[(nnx-1)+j*nnx+k*nnxnny])/tausT;
    }
  }
  //y fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=1;i<(nnx-1);i++)
  {
    Ics[i+0*nnx+0*nnxnny]=(cs[i+0*nnx+0*nnxnny+1]+cs[i+0*nnx+0*nnxnny-1]-2*cs[i+0*nnx+0*nnxnny])/tausL+
      (cs[i+(0+1)*nnx+0*nnxnny]+cs[i+(nny-1)*nnx+0*nnxnny]-2*cs[i+0*nnx+0*nnxnny])/tausT+
      (cs[i+0*nnx+(0+1)*nnxnny]+cs[i+0*nnx+(nnz-1)*nnxnny]-2*cs[i+0*nnx+0*nnxnny])/tausT;
    Ics[i+(nny-1)*nnx+0*nnxnny]=(cs[i+(nny-1)*nnx+0*nnxnny+1]+cs[i+(nny-1)*nnx+0*nnxnny-1]-2*cs[i+(nny-1)*nnx+0*nnxnny])/tausL+
      (cs[i+((nny-1)-1)*nnx+0*nnxnny]+cs[i+(0)*nnx+0*nnxnny]-2*cs[i+(nny-1)*nnx+0*nnxnny])/tausT+
      (cs[i+(nny-1)*nnx+(0+1)*nnxnny]+cs[i+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cs[i+(nny-1)*nnx+0*nnxnny])/tausT;
    Ics[i+0*nnx+(nnz-1)*nnxnny]=(cs[i+0*nnx+(nnz-1)*nnxnny+1]+cs[i+0*nnx+(nnz-1)*nnxnny-1]-2*cs[i+0*nnx+(nnz-1)*nnxnny])/tausL+
      (cs[i+(0+1)*nnx+(nnz-1)*nnxnny]+cs[i+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cs[i+0*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[i+0*nnx+((nnz-1)-1)*nnxnny]+cs[i+0*nnx+(0)*nnxnny]-2*cs[i+0*nnx+(nnz-1)*nnxnny])/tausT;
    Ics[i+(nny-1)*nnx+(nnz-1)*nnxnny]=(cs[i+(nny-1)*nnx+(nnz-1)*nnxnny+1]+cs[i+(nny-1)*nnx+(nnz-1)*nnxnny-1]-2*cs[i+(nny-1)*nnx+(nnz-1)*nnxnny])/tausL+
      (cs[i+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cs[i+(0)*nnx+(nnz-1)*nnxnny]-2*cs[i+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[i+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cs[i+(nny-1)*nnx+(0)*nnxnny]-2*cs[i+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      Ics[i+0*nnx+k*nnxnny]=(cs[i+0*nnx+k*nnxnny+1]+cs[i+0*nnx+k*nnxnny-1]-2*cs[i+0*nnx+k*nnxnny])/tausL+
        (cs[i+(0+1)*nnx+k*nnxnny]+cs[i+(nny-1)*nnx+k*nnxnny]-2*cs[i+0*nnx+k*nnxnny])/tausT+
        (cs[i+0*nnx+(k+1)*nnxnny]+cs[i+0*nnx+(k-1)*nnxnny]-2*cs[i+0*nnx+k*nnxnny])/tausT;
      Ics[i+(nny-1)*nnx+k*nnxnny]=(cs[i+(nny-1)*nnx+k*nnxnny+1]+cs[i+(nny-1)*nnx+k*nnxnny-1]-2*cs[i+(nny-1)*nnx+k*nnxnny])/tausL+
        (cs[i+((nny-1)-1)*nnx+k*nnxnny]+cs[i+(0)*nnx+k*nnxnny]-2*cs[i+(nny-1)*nnx+k*nnxnny])/tausT+
        (cs[i+(nny-1)*nnx+(k+1)*nnxnny]+cs[i+(nny-1)*nnx+(k-1)*nnxnny]-2*cs[i+(nny-1)*nnx+k*nnxnny])/tausT;
    }
    //z fixed
#pragma ivdep
#pragma vector always
    for (int j=1;j<nny-1;j++)
    {
      Ics[i+j*nnx+0*nnxnny]=(cs[i+j*nnx+0*nnxnny+1]+cs[i+j*nnx+0*nnxnny-1]-2*cs[i+j*nnx+0*nnxnny])/tausL+
        (cs[i+(j+1)*nnx+0*nnxnny]+cs[i+(j-1)*nnx+0*nnxnny]-2*cs[i+j*nnx+0*nnxnny])/tausT+
        (cs[i+j*nnx+(0+1)*nnxnny]+cs[i+j*nnx+(nnz-1)*nnxnny]-2*cs[i+j*nnx+0*nnxnny])/tausT;
      Ics[i+j*nnx+(nnz-1)*nnxnny]=(cs[i+j*nnx+(nnz-1)*nnxnny+1]+cs[i+j*nnx+(nnz-1)*nnxnny-1]-2*cs[i+j*nnx+(nnz-1)*nnxnny])/tausL+
        (cs[i+(j+1)*nnx+(nnz-1)*nnxnny]+cs[i+(j-1)*nnx+(nnz-1)*nnxnny]-2*cs[i+j*nnx+(nnz-1)*nnxnny])/tausT+
        (cs[i+j*nnx+((nnz-1)-1)*nnxnny]+cs[i+j*nnx+(0)*nnxnny]-2*cs[i+j*nnx+(nnz-1)*nnxnny])/tausT;
    }
  }
#pragma ivdep
#pragma vector always
  for (int k=1;k<nnz-1;k++)
  {
    Ics[0+0*nnx+k*nnxnny]=(cs[0+0*nnx+k*nnxnny+1]+cs[(nnx-1)+0*nnx+k*nnxnny]-2*cs[0+0*nnx+k*nnxnny])/tausL+
      (cs[0+(0+1)*nnx+k*nnxnny]+cs[0+(nny-1)*nnx+k*nnxnny]-2*cs[0+0*nnx+k*nnxnny])/tausT+
      (cs[0+0*nnx+(k+1)*nnxnny]+cs[0+0*nnx+(k-1)*nnxnny]-2*cs[0+0*nnx+k*nnxnny])/tausT;
    Ics[0+(nny-1)*nnx+k*nnxnny]=(cs[0+(nny-1)*nnx+k*nnxnny+1]+cs[(nnx-1)+(nny-1)*nnx+k*nnxnny]-2*cs[0+(nny-1)*nnx+k*nnxnny])/tausL+
      (cs[0+((nny-1)-1)*nnx+k*nnxnny]+cs[0+(0)*nnx+k*nnxnny]-2*cs[0+(nny-1)*nnx+k*nnxnny])/tausT+
      (cs[0+(nny-1)*nnx+(k+1)*nnxnny]+cs[0+(nny-1)*nnx+(k-1)*nnxnny]-2*cs[0+(nny-1)*nnx+k*nnxnny])/tausT;
    Ics[(nnx-1)+0*nnx+k*nnxnny]=(cs[(nnx-1)+0*nnx+k*nnxnny-1]+cs[(0)+0*nnx+k*nnxnny]-2*cs[(nnx-1)+0*nnx+k*nnxnny])/tausL+
      (cs[(nnx-1)+(0+1)*nnx+k*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+k*nnxnny]-2*cs[(nnx-1)+0*nnx+k*nnxnny])/tausT+
      (cs[(nnx-1)+0*nnx+(k+1)*nnxnny]+cs[(nnx-1)+0*nnx+(k-1)*nnxnny]-2*cs[(nnx-1)+0*nnx+k*nnxnny])/tausT;
    Ics[(nnx-1)+(nny-1)*nnx+k*nnxnny]=(cs[(nnx-1)+(nny-1)*nnx+k*nnxnny-1]+cs[(0)+(nny-1)*nnx+k*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+k*nnxnny])/tausL+
      (cs[(nnx-1)+((nny-1)-1)*nnx+k*nnxnny]+cs[(nnx-1)+(0)*nnx+k*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+k*nnxnny])/tausT+
      (cs[(nnx-1)+(nny-1)*nnx+(k+1)*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+(k-1)*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+k*nnxnny])/tausT;
  }
#else
  Ics[0+0*nnx+0*nnxnny]=(cs[0+0*nnx+0*nnxnny+1]+cs[0+0*nnx+0*nnxnny+1]-2*cs[0+0*nnx+0*nnxnny])/tausL+
    (cs[0+(0+1)*nnx+0*nnxnny]+cs[0+(0+1)*nnx+0*nnxnny]-2*cs[0+0*nnx+0*nnxnny])/tausT+
    (cs[0+0*nnx+(0+1)*nnxnny]+cs[0+0*nnx+(0+1)*nnxnny]-2*cs[0+0*nnx+0*nnxnny])/tausT;
  Ics[0+(nny-1)*nnx+(nnz-1)*nnxnny]=(cs[0+(nny-1)*nnx+(nnz-1)*nnxnny+1]+cs[0+(nny-1)*nnx+(nnz-1)*nnxnny+1]-2*cs[0+(nny-1)*nnx+(nnz-1)*nnxnny])/tausL+
    (cs[0+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cs[0+((nny-1)-1)*nnx+(nnz-1)*nnxnny]-2*cs[0+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[0+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cs[0+(nny-1)*nnx+((nnz-1)-1)*nnxnny]-2*cs[0+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT;
  Ics[0+(nny-1)*nnx+0*nnxnny]=(cs[0+(nny-1)*nnx+0*nnxnny+1]+cs[0+(nny-1)*nnx+0*nnxnny+1]-2*cs[0+(nny-1)*nnx+0*nnxnny])/tausL+
    (cs[0+((nny-1)-1)*nnx+0*nnxnny]+cs[0+((nny-1)-1)*nnx+0*nnxnny]-2*cs[0+(nny-1)*nnx+0*nnxnny])/tausT+
    (cs[0+(nny-1)*nnx+(0+1)*nnxnny]+cs[0+(nny-1)*nnx+(0+1)*nnxnny]-2*cs[0+(nny-1)*nnx+0*nnxnny])/tausT;
  Ics[0+0*nnx+(nnz-1)*nnxnny]=(cs[0+0*nnx+(nnz-1)*nnxnny+1]+cs[0+0*nnx+(nnz-1)*nnxnny+1]-2*cs[0+0*nnx+(nnz-1)*nnxnny])/tausL+
    (cs[0+(0+1)*nnx+(nnz-1)*nnxnny]+cs[0+(0+1)*nnx+(nnz-1)*nnxnny]-2*cs[0+0*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[0+0*nnx+((nnz-1)-1)*nnxnny]+cs[0+0*nnx+((nnz-1)-1)*nnxnny]-2*cs[0+0*nnx+(nnz-1)*nnxnny])/tausT;
  Ics[(nnx-1)+(nny-1)*nnx+0*nnxnny]=(cs[(nnx-1)+(nny-1)*nnx+0*nnxnny-1]+cs[(nnx-1)+(nny-1)*nnx+0*nnxnny-1]-2*cs[(nnx-1)+(nny-1)*nnx+0*nnxnny])/tausL+
    (cs[(nnx-1)+((nny-1)-1)*nnx+0*nnxnny]+cs[(nnx-1)+((nny-1)-1)*nnx+0*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+0*nnxnny])/tausT+
    (cs[(nnx-1)+(nny-1)*nnx+(0+1)*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+(0+1)*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+0*nnxnny])/tausT;
  Ics[(nnx-1)+0*nnx+(nnz-1)*nnxnny]=(cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny-1]+cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny-1]-2*cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/tausL+
    (cs[(nnx-1)+(0+1)*nnx+(nnz-1)*nnxnny]+cs[(nnx-1)+(0+1)*nnx+(nnz-1)*nnxnny]-2*cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[(nnx-1)+0*nnx+((nnz-1)-1)*nnxnny]+cs[(nnx-1)+0*nnx+((nnz-1)-1)*nnxnny]-2*cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/tausT;
  Ics[(nnx-1)+0*nnx+0*nnxnny]=(cs[(nnx-1)+0*nnx+0*nnxnny-1]+cs[(nnx-1)+0*nnx+0*nnxnny-1]-2*cs[(nnx-1)+0*nnx+0*nnxnny])/tausL+
    (cs[(nnx-1)+(0+1)*nnx+0*nnxnny]+cs[(nnx-1)+(0+1)*nnx+0*nnxnny]-2*cs[(nnx-1)+0*nnx+0*nnxnny])/tausT+
    (cs[(nnx-1)+0*nnx+(0+1)*nnxnny]+cs[(nnx-1)+0*nnx+(0+1)*nnxnny]-2*cs[(nnx-1)+0*nnx+0*nnxnny])/tausT;
  Ics[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]=(cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny-1]+cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny-1]-2*cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tausL+
    (cs[(nnx-1)+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cs[(nnx-1)+((nny-1)-1)*nnx+(nnz-1)*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[(nnx-1)+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+((nnz-1)-1)*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT;
  //x fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int j=1;j<nny-1;j++)
  {
    Ics[0+j*nnx+0*nnxnny]=(cs[0+j*nnx+0*nnxnny+1]+cs[0+j*nnx+0*nnxnny+1]-2*cs[0+j*nnx+0*nnxnny])/tausL+
      (cs[0+(j+1)*nnx+0*nnxnny]+cs[0+(j-1)*nnx+0*nnxnny]-2*cs[0+j*nnx+0*nnxnny])/tausT+
      (cs[0+j*nnx+(0+1)*nnxnny]+cs[0+j*nnx+(0+1)*nnxnny]-2*cs[0+j*nnx+0*nnxnny])/tausT;
    Ics[(nnx-1)+j*nnx+0*nnxnny]=(cs[(nnx-1)+j*nnx+0*nnxnny-1]+cs[(nnx-1)+j*nnx+0*nnxnny-1]-2*cs[(nnx-1)+j*nnx+0*nnxnny])/tausL+
      (cs[(nnx-1)+(j+1)*nnx+0*nnxnny]+cs[(nnx-1)+(j-1)*nnx+0*nnxnny]-2*cs[(nnx-1)+j*nnx+0*nnxnny])/tausT+
      (cs[(nnx-1)+j*nnx+(0+1)*nnxnny]+cs[(nnx-1)+j*nnx+(0+1)*nnxnny]-2*cs[(nnx-1)+j*nnx+0*nnxnny])/tausT;
    Ics[0+j*nnx+(nnz-1)*nnxnny]=(cs[0+j*nnx+(nnz-1)*nnxnny+1]+cs[0+j*nnx+(nnz-1)*nnxnny+1]-2*cs[0+j*nnx+(nnz-1)*nnxnny])/tausL+
      (cs[0+(j+1)*nnx+(nnz-1)*nnxnny]+cs[0+(j-1)*nnx+(nnz-1)*nnxnny]-2*cs[0+j*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[0+j*nnx+((nnz-1)-1)*nnxnny]+cs[0+j*nnx+((nnz-1)-1)*nnxnny]-2*cs[0+j*nnx+(nnz-1)*nnxnny])/tausT;
    Ics[(nnx-1)+j*nnx+(nnz-1)*nnxnny]=(cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny-1]+cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny-1]-2*cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/tausL+
      (cs[(nnx-1)+(j+1)*nnx+(nnz-1)*nnxnny]+cs[(nnx-1)+(j-1)*nnx+(nnz-1)*nnxnny]-2*cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[(nnx-1)+j*nnx+((nnz-1)-1)*nnxnny]+cs[(nnx-1)+j*nnx+((nnz-1)-1)*nnxnny]-2*cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/tausT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      Ics[0+j*nnx+k*nnxnny]=(cs[0+j*nnx+k*nnxnny+1]+cs[0+j*nnx+k*nnxnny+1]-2*cs[0+j*nnx+k*nnxnny])/tausL+
        (cs[0+(j+1)*nnx+k*nnxnny]+cs[0+(j-1)*nnx+k*nnxnny]-2*cs[0+j*nnx+k*nnxnny])/tausT+
        (cs[0+j*nnx+(k+1)*nnxnny]+cs[0+j*nnx+(k-1)*nnxnny]-2*cs[0+j*nnx+k*nnxnny])/tausT;
      Ics[(nnx-1)+j*nnx+k*nnxnny]=(cs[(nnx-1)+j*nnx+k*nnxnny-1]+cs[(nnx-1)+j*nnx+k*nnxnny-1]-2*cs[(nnx-1)+j*nnx+k*nnxnny])/tausL+
        (cs[(nnx-1)+(j+1)*nnx+k*nnxnny]+cs[(nnx-1)+(j-1)*nnx+k*nnxnny]-2*cs[(nnx-1)+j*nnx+k*nnxnny])/tausT+
        (cs[(nnx-1)+j*nnx+(k+1)*nnxnny]+cs[(nnx-1)+j*nnx+(k-1)*nnxnny]-2*cs[(nnx-1)+j*nnx+k*nnxnny])/tausT;
    }
  }
  //y fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=1;i<(nnx-1);i++)
  {
    Ics[i+0*nnx+0*nnxnny]=(cs[i+0*nnx+0*nnxnny+1]+cs[i+0*nnx+0*nnxnny-1]-2*cs[i+0*nnx+0*nnxnny])/tausL+
      (cs[i+(0+1)*nnx+0*nnxnny]+cs[i+(0+1)*nnx+0*nnxnny]-2*cs[i+0*nnx+0*nnxnny])/tausT+
      (cs[i+0*nnx+(0+1)*nnxnny]+cs[i+0*nnx+(0+1)*nnxnny]-2*cs[i+0*nnx+0*nnxnny])/tausT;
    Ics[i+(nny-1)*nnx+0*nnxnny]=(cs[i+(nny-1)*nnx+0*nnxnny+1]+cs[i+(nny-1)*nnx+0*nnxnny-1]-2*cs[i+(nny-1)*nnx+0*nnxnny])/tausL+
      (cs[i+((nny-1)-1)*nnx+0*nnxnny]+cs[i+((nny-1)-1)*nnx+0*nnxnny]-2*cs[i+(nny-1)*nnx+0*nnxnny])/tausT+
      (cs[i+(nny-1)*nnx+(0+1)*nnxnny]+cs[i+(nny-1)*nnx+(0+1)*nnxnny]-2*cs[i+(nny-1)*nnx+0*nnxnny])/tausT;
    Ics[i+0*nnx+(nnz-1)*nnxnny]=(cs[i+0*nnx+(nnz-1)*nnxnny+1]+cs[i+0*nnx+(nnz-1)*nnxnny-1]-2*cs[i+0*nnx+(nnz-1)*nnxnny])/tausL+
      (cs[i+(0+1)*nnx+(nnz-1)*nnxnny]+cs[i+(0+1)*nnx+(nnz-1)*nnxnny]-2*cs[i+0*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[i+0*nnx+((nnz-1)-1)*nnxnny]+cs[i+0*nnx+((nnz-1)-1)*nnxnny]-2*cs[i+0*nnx+(nnz-1)*nnxnny])/tausT;
    Ics[i+(nny-1)*nnx+(nnz-1)*nnxnny]=(cs[i+(nny-1)*nnx+(nnz-1)*nnxnny+1]+cs[i+(nny-1)*nnx+(nnz-1)*nnxnny-1]-2*cs[i+(nny-1)*nnx+(nnz-1)*nnxnny])/tausL+
      (cs[i+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cs[i+((nny-1)-1)*nnx+(nnz-1)*nnxnny]-2*cs[i+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[i+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cs[i+(nny-1)*nnx+((nnz-1)-1)*nnxnny]-2*cs[i+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      Ics[i+0*nnx+k*nnxnny]=(cs[i+0*nnx+k*nnxnny+1]+cs[i+0*nnx+k*nnxnny-1]-2*cs[i+0*nnx+k*nnxnny])/tausL+
        (cs[i+(0+1)*nnx+k*nnxnny]+cs[i+(0+1)*nnx+k*nnxnny]-2*cs[i+0*nnx+k*nnxnny])/tausT+
        (cs[i+0*nnx+(k+1)*nnxnny]+cs[i+0*nnx+(k-1)*nnxnny]-2*cs[i+0*nnx+k*nnxnny])/tausT;
      Ics[i+(nny-1)*nnx+k*nnxnny]=(cs[i+(nny-1)*nnx+k*nnxnny+1]+cs[i+(nny-1)*nnx+k*nnxnny-1]-2*cs[i+(nny-1)*nnx+k*nnxnny])/tausL+
        (cs[i+((nny-1)-1)*nnx+k*nnxnny]+cs[i+((nny-1)-1)*nnx+k*nnxnny]-2*cs[i+(nny-1)*nnx+k*nnxnny])/tausT+
        (cs[i+(nny-1)*nnx+(k+1)*nnxnny]+cs[i+(nny-1)*nnx+(k-1)*nnxnny]-2*cs[i+(nny-1)*nnx+k*nnxnny])/tausT;
    }
    //z fixed
#pragma ivdep
#pragma vector always
    for (int j=1;j<nny-1;j++)
    {
      Ics[i+j*nnx+0*nnxnny]=(cs[i+j*nnx+0*nnxnny+1]+cs[i+j*nnx+0*nnxnny-1]-2*cs[i+j*nnx+0*nnxnny])/tausL+
        (cs[i+(j+1)*nnx+0*nnxnny]+cs[i+(j-1)*nnx+0*nnxnny]-2*cs[i+j*nnx+0*nnxnny])/tausT+
        (cs[i+j*nnx+(0+1)*nnxnny]+cs[i+j*nnx+(0+1)*nnxnny]-2*cs[i+j*nnx+0*nnxnny])/tausT;
      Ics[i+j*nnx+(nnz-1)*nnxnny]=(cs[i+j*nnx+(nnz-1)*nnxnny+1]+cs[i+j*nnx+(nnz-1)*nnxnny-1]-2*cs[i+j*nnx+(nnz-1)*nnxnny])/tausL+
        (cs[i+(j+1)*nnx+(nnz-1)*nnxnny]+cs[i+(j-1)*nnx+(nnz-1)*nnxnny]-2*cs[i+j*nnx+(nnz-1)*nnxnny])/tausT+
        (cs[i+j*nnx+((nnz-1)-1)*nnxnny]+cs[i+j*nnx+((nnz-1)-1)*nnxnny]-2*cs[i+j*nnx+(nnz-1)*nnxnny])/tausT;
    }
  }
#pragma ivdep
#pragma vector always
  for (int k=1;k<nnz-1;k++)
  {
    Ics[0+0*nnx+k*nnxnny]=(cs[0+0*nnx+k*nnxnny+1]+cs[0+0*nnx+k*nnxnny+1]-2*cs[0+0*nnx+k*nnxnny])/tausL+
      (cs[0+(0+1)*nnx+k*nnxnny]+cs[0+(0+1)*nnx+k*nnxnny]-2*cs[0+0*nnx+k*nnxnny])/tausT+
      (cs[0+0*nnx+(k+1)*nnxnny]+cs[0+0*nnx+(k-1)*nnxnny]-2*cs[0+0*nnx+k*nnxnny])/tausT;
    Ics[0+(nny-1)*nnx+k*nnxnny]=(cs[0+(nny-1)*nnx+k*nnxnny+1]+cs[0+(nny-1)*nnx+k*nnxnny+1]-2*cs[0+(nny-1)*nnx+k*nnxnny])/tausL+
      (cs[0+((nny-1)-1)*nnx+k*nnxnny]+cs[0+((nny-1)-1)*nnx+k*nnxnny]-2*cs[0+(nny-1)*nnx+k*nnxnny])/tausT+
      (cs[0+(nny-1)*nnx+(k+1)*nnxnny]+cs[0+(nny-1)*nnx+(k-1)*nnxnny]-2*cs[0+(nny-1)*nnx+k*nnxnny])/tausT;
    Ics[(nnx-1)+0*nnx+k*nnxnny]=(cs[(nnx-1)+0*nnx+k*nnxnny-1]+cs[(nnx-1)+0*nnx+k*nnxnny-1]-2*cs[(nnx-1)+0*nnx+k*nnxnny])/tausL+
      (cs[(nnx-1)+(0+1)*nnx+k*nnxnny]+cs[(nnx-1)+(0+1)*nnx+k*nnxnny]-2*cs[(nnx-1)+0*nnx+k*nnxnny])/tausT+
      (cs[(nnx-1)+0*nnx+(k+1)*nnxnny]+cs[(nnx-1)+0*nnx+(k-1)*nnxnny]-2*cs[(nnx-1)+0*nnx+k*nnxnny])/tausT;
    Ics[(nnx-1)+(nny-1)*nnx+k*nnxnny]=(cs[(nnx-1)+(nny-1)*nnx+k*nnxnny-1]+cs[(nnx-1)+(nny-1)*nnx+k*nnxnny-1]-2*cs[(nnx-1)+(nny-1)*nnx+k*nnxnny])/tausL+
      (cs[(nnx-1)+((nny-1)-1)*nnx+k*nnxnny]+cs[(nnx-1)+((nny-1)-1)*nnx+k*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+k*nnxnny])/tausT+
      (cs[(nnx-1)+(nny-1)*nnx+(k+1)*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+(k-1)*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+k*nnxnny])/tausT;
  }

#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int k=1;k<nnz-1;k++)
  {
    for (int j=1;j<nny-1;j++)
    {
#pragma ivdep
#pragma vector always
      for (int i=1;i<nnx-1;i++)
      {
        Ics[i+j*nnx+k*nnxnny]=(cs[i+j*nnx+k*nnxnny+1]+cs[i+j*nnx+k*nnxnny-1]-2*cs[i+j*nnx+k*nnxnny])/tausL+
          (cs[i+(j+1)*nnx+k*nnxnny]+cs[i+(j-1)*nnx+k*nnxnny]-2*cs[i+j*nnx+k*nnxnny])/tausT+
          (cs[i+j*nnx+(k+1)*nnxnny]+cs[i+j*nnx+(k-1)*nnxnny]-2*cs[i+j*nnx+k*nnxnny])/tausT;
      }
    }
  }
}

#endif
