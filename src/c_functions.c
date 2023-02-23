#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>

int int_max(int a, int b)
{
  if(a>b)
    return(a);
  else
    return(b);
}

SEXP likelihood(SEXP Y,SEXP X,SEXP Delta,SEXP beta,SEXP sigma0,SEXP lambda,SEXP tau,SEXP theta,SEXP tilde_theta)
{
  // Variable definitions
  SEXP out;
  SEXP LH;
  SEXP beta_derivative;
  SEXP theta_derivative;
  SEXP theta_tilde_derivative;
  SEXP lambda_derivative;
  double* sigmasq;
  int n;
  int p;
  int M,m,tilde_m;
  int lag_order;
  int sum_order;
  int i,j,a,b,k,k1,k2,i1,i2;
  int* binom;
  int* factorials;
  double Hlambda;
  double laguerre,laguerre2;
  double likeli;
  double real_contribution;
  double* xbeta;
  double* H;
  double* evaluation_points;
  double* lag;
  double* I;
  double* density;
  double* support_vector;
  double* dtheta;
  double* dtheta_tilde;
  double* dsigmasq;

  // Initialize data
  n=LENGTH(Y);
  p=LENGTH(X)/n-1;
  if(p>0)
    M=LENGTH(lambda)/p-1;
  else
    M=0;
  m=LENGTH(theta)-1;
  tilde_m=LENGTH(tilde_theta)-1;
  lag_order=int_max(m,tilde_m);
  sum_order=2*lag_order;

  sigmasq=(double*)malloc(sizeof(double)*n);
  H=(double*)malloc(sizeof(double)*(M+1)*n);
  evaluation_points=(double*)malloc(sizeof(double)*n);
  lag=(double*)malloc(sizeof(double)*n*(lag_order+1));
  I=(double*)malloc(sizeof(double)*n*(sum_order+1));
  binom=(int*)malloc(sizeof(int)*(lag_order+1)*(lag_order+1));
  factorials=(int*)malloc(sizeof(int)*(lag_order+1));
  xbeta=(double*)malloc(sizeof(double)*n);
  density=(double*)malloc(sizeof(double)*n);
  support_vector=(double*)malloc(sizeof(double)*n*(p+1));
  dtheta=(double*)malloc(sizeof(double)*(m+1));
  dtheta_tilde=(double*)malloc(sizeof(double)*(tilde_m+1));
  dsigmasq=(double*)malloc(sizeof(double)*(M+1)*p*n);



  // Compute sigma^2(X) for each individual covariate other than the intercept
  for(i=0;i<=n-1;i++)
    sigmasq[i]=1;

  if(p>0)
    for(i=0;i<=n-1;i++)
      for(a=1;a<=p;a++)
        for(k=0;k<=M;k++)
          dsigmasq[i+(a-1)*n+k*n*p]=1;

  if(p>0)
  {
    for(a=1;a<=p;a++)
    {
      // First Hermite Polynomial equals one
      for(i=0;i<=n-1;i++)
        H[i]=1;

      // Compute higher order Hermite polynomials for covariate a by using a recurrence relation
      if(M+1>=2)
      {
        for(i=0;i<=n-1;i++)
          H[i+n]=REAL(X)[i+a*n];

        if(M+1>=3)
        {
          for(k=2;k<=M;k++)
          {
            for(i=0;i<=n-1;i++)
              H[i+k*n]=REAL(X)[i+a*n]/sqrt(k)*H[i+(k-1)*n]-sqrt((double)(k-1)/(double)k)*H[i+(k-2)*n];
          }
        }
      }


      // Multiply the result to sigmasq according to its definition and compute dsigmasq
      for(i=0;i<=n-1;i++)
      {
        Hlambda=0;
        for(k=0;k<=M;k++)
          Hlambda=Hlambda+H[i+k*n]*REAL(lambda)[k+(a-1)*(M+1)];

        sigmasq[i]=sigmasq[i]*Hlambda*Hlambda;

        // Compute dsigmasq
        for(k=0;k<=M;k++)
          for(b=1;b<=p;b++)
            if(a!=b)
              dsigmasq[i+(b-1)*n+k*n*p]=dsigmasq[i+(b-1)*n+k*n*p]*Hlambda*Hlambda;
            else
              dsigmasq[i+(b-1)*n+k*n*p]=dsigmasq[i+(b-1)*n+k*n*p]*2*Hlambda*H[i+k*n];
      }
    }
  }

  // Add sigma0
  for(i=0;i<=n-1;i++)
    sigmasq[i]=REAL(sigma0)[0]+sigmasq[i];


  // Compute xbeta
  for(i=0;i<=n-1;i++)
  {
    xbeta[i]=0;
    for(j=0;j<=p;j++)
      xbeta[i]=xbeta[i]+REAL(beta)[j]*REAL(X)[i+j*n];
  }


  // Compute evaluation points
  for(i=0;i<=n-1;i++)
  {
    if(REAL(Y)[i]>xbeta[i])
      evaluation_points[i]=REAL(tau)[0]*(REAL(Y)[i]-xbeta[i])/sigmasq[i];
    else
      evaluation_points[i]=-(1-REAL(tau)[0])*(REAL(Y)[i]-xbeta[i])/sigmasq[i];
  }

  // Compute Laguerre Polynomials using a recurrence relation
  for(i=0;i<=n-1;i++)
  {
    lag[i]=1;
    if(lag_order>=1)
    {
      lag[i+n]=1-evaluation_points[i];
      if(lag_order>=2)
      {
        for(k=2;k<=lag_order;k++)
          lag[i+k*n]=((2*k-1-evaluation_points[i])*lag[i+(k-1)*n]-(k-1)*lag[i+(k-2)*n])/k;
      }
    }
  }

  // Compute exponential integral using a recurrence relation
  for(i=0;i<=n-1;i++)
  {
    I[i]=exp(-evaluation_points[i]);
    if(sum_order>=1)
    {
      for(k=1;k<=sum_order;k++)
        I[i+k*n]=pow(-evaluation_points[i],(double)k)*exp(-evaluation_points[i])-(double)k*I[i+(k-1)*n];
    }
  }

  // Compute binomial coefficients
  for(k1=0;k1<=lag_order;k1++)
  {
    for(k2=0;k2<=lag_order;k2++)
    {
      if(k2==0)
        binom[k1+k2*(lag_order+1)]=1;
      else if(k1==0)
        binom[k1+k2*(lag_order+1)]=0;
      else
        binom[k1+k2*(lag_order+1)]=binom[(k1-1)+k2*(lag_order+1)]+binom[(k1-1)+(k2-1)*(lag_order+1)];
    }
  }

  // Compute factorials
  factorials[0]=1;
  if(lag_order>=1)
    for(k=1;k<=lag_order;k++)
      factorials[k]=k*factorials[k-1];

  // Compute density and its derivative
  for(i=0;i<=n-1;i++) {
    laguerre=0;
    laguerre2=0;
    if(REAL(Y)[i]<xbeta[i])
    {
      // We are on the left branch: Use tilde:
      for(k=0;k<=tilde_m;k++)
      {
        laguerre=laguerre+REAL(tilde_theta)[k]*lag[i+k*n];
        if(k>0)
          laguerre2=laguerre2+REAL(tilde_theta)[k]*(double)k*(lag[i+k*n]-lag[i+(k-1)*n])/evaluation_points[i];
      }
    } else {
      // We are on the right branch: Use no tilde:
      for(k=0;k<=m;k++)
      {
        laguerre=laguerre+REAL(theta)[k]*lag[i+k*n];
        if(k>0)
          laguerre2=laguerre2+REAL(theta)[k]*(double)k*(lag[i+k*n]-lag[i+(k-1)*n])/evaluation_points[i];
      }
    }
    density[i]=REAL(tau)[0]*(1-REAL(tau)[0])*exp(-evaluation_points[i])*laguerre*laguerre/sigmasq[i];
    if(REAL(Y)[i]<xbeta[i])
      support_vector[i]=(1-REAL(tau)[0])*(density[i]-2*REAL(tau)[0]*(1-REAL(tau)[0])*exp(-evaluation_points[i])*laguerre*laguerre2/sigmasq[i]);
    else
      support_vector[i]=  -REAL(tau)[0] *(density[i]-2*REAL(tau)[0]*(1-REAL(tau)[0])*exp(-evaluation_points[i])*laguerre*laguerre2/sigmasq[i]);
  }

  // Compute Likelihood
  likeli=0;
  for(i=0;i<=n-1;i++)
  {
    real_contribution=0;
    if(INTEGER(Delta)[i]==1)
    {
      // Uncensored observation: Use density
      real_contribution=density[i];
    } else {
      // Censored Observation: Compute distribution
      laguerre=0.0;
      if(REAL(Y)[i]<xbeta[i])
      {
        // We are in the left branch: Use tilde
        for(k1=0;k1<=tilde_m;k1++)
          for(k2=0;k2<=tilde_m;k2++)
            for(i1=0;i1<=k1;i1++)
              for(i2=0;i2<=k2;i2++)
                laguerre=laguerre+1.0*binom[k1+i1*(lag_order+1)]*binom[k2+i2*(lag_order+1)]/factorials[i1]/factorials[i2]*I[i+(i1+i2)*n]*REAL(tilde_theta)[k1]*REAL(tilde_theta)[k2];

        real_contribution=1-REAL(tau)[0]*laguerre;
      } else {
        // We are in the left branch: Use not tilde
        for(k1=0;k1<=m;k1++)
          for(k2=0;k2<=m;k2++)
            for(i1=0;i1<=k1;i1++)
              for(i2=0;i2<=k2;i2++)
                laguerre=laguerre+1.0*binom[k1+i1*(lag_order+1)]*binom[k2+i2*(lag_order+1)]/factorials[i1]/factorials[i2]*I[i+(i1+i2)*n]*REAL(theta)[k1]*REAL(theta)[k2];

        real_contribution=(1-REAL(tau)[0])*laguerre;
      }
      support_vector[i]=real_contribution;
    }
//    if(real_contribution<=0) {
 //     fprintf(stderr,"ERROR!!!\n");
//      fprintf(stderr,"real_contribution=%f\n",real_contribution);
//      fprintf(stderr,"i=%d, Y=%f, xbeta=%f, x=%f, Delta=%d, %f\n",i,REAL(Y)[i],xbeta[i],REAL(X)[i+n],INTEGER(Delta)[i],1.0*(i+1)/i);
//    }
    likeli=likeli+log(real_contribution);
  }

  // Compute derivative with respect to theta
  for(j=0;j<=m;j++)
  {
    dtheta[j]=0;
    for(i=0;i<=n-1;i++)
      if(REAL(Y)[i]>xbeta[i])
      {
        if(INTEGER(Delta)[i]==1)
        {
          for(k=0;k<=m;k++)
            dtheta[j]=dtheta[j]+2*REAL(tau)[0]*(1-REAL(tau)[0])*exp(-evaluation_points[i])*REAL(theta)[k]*lag[i+k*n]*lag[i+j*n]/sigmasq[i]/density[i];
        }
        else
        {
          for(k=0;k<=m;k++)
            for(i1=0;i1<=k;i1++)
              for(i2=0;i2<=j;i2++)
                dtheta[j]=dtheta[j]+2*(1-REAL(tau)[0])*REAL(theta)[k]*binom[k+i1*(lag_order+1)]*binom[j+i2*(lag_order+1)]/factorials[i1]/factorials[i2]*I[i+(i1+i2)*n]/support_vector[i];
        }
      }
  }

  // Compute derivative with respect to theta_tilde
  for(j=0;j<=tilde_m;j++)
  {
    dtheta_tilde[j]=0;
    for(i=0;i<=n-1;i++)
      if(REAL(Y)[i]<=xbeta[i])
      {
        if(INTEGER(Delta)[i]==1)
        {
          for(k=0;k<=tilde_m;k++)
            dtheta_tilde[j]=dtheta_tilde[j]+2*REAL(tau)[0]*(1-REAL(tau)[0])*exp(-evaluation_points[i])*REAL(tilde_theta)[k]*lag[i+k*n]*lag[i+j*n]/sigmasq[i]/density[i];
        }
        else
        {
          for(k=0;k<=tilde_m;k++)
            for(i1=0;i1<=k;i1++)
              for(i2=0;i2<=j;i2++)
                dtheta_tilde[j]=dtheta_tilde[j]-2*REAL(tau)[0]*REAL(tilde_theta)[k]*binom[k+i1*(lag_order+1)]*binom[j+i2*(lag_order+1)]/factorials[i1]/factorials[i2]*I[i+(i1+i2)*n]/support_vector[i];
        }
      }
  }

  // Prepare output
  // Objective
  LH=PROTECT(allocVector(REALSXP,1));
  REAL(LH)[0]=likeli/n;

  // Derivative with respect to beta
  beta_derivative=PROTECT(allocVector(REALSXP,p+1));
  for(k=0;k<=p;k++)
  {
    REAL(beta_derivative)[k]=0;
    for(i=0;i<=n-1;i++)
    {
      if(INTEGER(Delta)[i]==0)
        REAL(beta_derivative)[k]=REAL(beta_derivative)[k]+REAL(X)[i+k*n]*density[i]/support_vector[i]/n;
      else
        REAL(beta_derivative)[k]=REAL(beta_derivative)[k]-REAL(X)[i+k*n]*support_vector[i]/density[i]/sigmasq[i]/n;
    }
  }

  // Derivatives with respect to theta and tilde_theta
  theta_derivative=PROTECT(allocVector(REALSXP,m+1));
  for(k=0;k<=m;k++)
    REAL(theta_derivative)[k]=dtheta[k]/n;

  theta_tilde_derivative=PROTECT(allocVector(REALSXP,tilde_m+1));
  for(k=0;k<=tilde_m;k++)
    REAL(theta_tilde_derivative)[k]=dtheta_tilde[k]/n;

  // Derivative with respect to lambda
  if(p>0)
  {
    lambda_derivative=PROTECT(allocVector(REALSXP,(M+1)*p));
    for(b=1;b<=p;b++)
    {
      for(k=0;k<=M;k++)
      {
        REAL(lambda_derivative)[k+(b-1)*(M+1)]=0;
        for(i=0;i<=n-1;i++)
          if(INTEGER(Delta)[i]==1)
            REAL(lambda_derivative)[k+(b-1)*(M+1)]=REAL(lambda_derivative)[k+(b-1)*(M+1)]-dsigmasq[i+(b-1)*n+k*n*p]*(1/sigmasq[i]+(REAL(Y)[i]-xbeta[i])*support_vector[i]/sigmasq[i]/sigmasq[i]/density[i])/n;
          else
            REAL(lambda_derivative)[k+(b-1)*(M+1)]=REAL(lambda_derivative)[k+(b-1)*(M+1)]+dsigmasq[i+(b-1)*n+k*n*p]*(REAL(Y)[i]-xbeta[i])*density[i]/support_vector[i]/sigmasq[i]/n;
      }
    }
  }

  if(p>0)
    out=PROTECT(allocVector(VECSXP,5));
  else
    out=PROTECT(allocVector(VECSXP,4));

  SET_VECTOR_ELT(out,0,LH);
  SET_VECTOR_ELT(out,1,beta_derivative);
  SET_VECTOR_ELT(out,2,theta_derivative);
  SET_VECTOR_ELT(out,3,theta_tilde_derivative);

  if(p>0)
    SET_VECTOR_ELT(out,4,lambda_derivative);

  // Free used data
  free(xbeta);
  free(H);
  free(evaluation_points);
  free(sigmasq);
  free(lag);
  free(I);
  free(binom);
  free(factorials);
  free(density);
  free(support_vector);
  free(dtheta);
  free(dtheta_tilde);
  free(dsigmasq);

  if(p>0)
    UNPROTECT(6);
  else
    UNPROTECT(5);

  return(out);

}
