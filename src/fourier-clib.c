/**********************************************************************
 * Integral Transform Method for Sufficient Dimension Reduction
 * Rcmd SHLIB fourier-clib.c -lRblas (Windows)
 * Reference: 
 *   Peng Zeng and Yu Zhu (2010). An integral transform method for 
 *      estimating the central mean and central Subspaces.
 *      Journal of Multivariate Analysis 101, 271--290.
 *
 * Peng Zeng @ Auburn University
 * 03-30-2010
 **********************************************************************/
  
#include <R.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>

/**********************************************************************
 * Calculate the candidate matrix for the Central Mean Subspace.
 * Weight function W(x, u): Gaussian function with covariance h * I_p.
 *
 * Input:  
 *    x  --  n-by-p matrix  (predictors)
 *    y  --  n-by-1 vector  (response)
 *    gk --  n-by-p matrix  (derivative of log f(x))
 *    h  --  scalar         (sigma^2)
 *    out--  1 or 0         (1 = leave-i-out, 0 = not)
 *
 * Output:
 *    Mmat -- p-by-p symmetric matrix  (candidate matrix)
 *            notice: only the lower triangle of Mmat are calculated. 
 **********************************************************************/

void FM_mean(const double* x, const double* gk, const double* y,
              const double* h, const int* n, const int* p, 
              const int* out, double* Mmat)
{
  const char uplo = 'L';
  const double sc1 = (*h), sc2 = -0.5*(*h);
  const double a = 1.0, b = -1.0;
  const int pp = (*p) * (*p), one = 1;
  double alpha, wij, wsum; 
  double xij[*p], gkxij[*p]; 
  int i, j;

  for(i = 0; i < pp; i++)  Mmat[i] = 0.0;

  wsum = 0.0; 
  if((*out) == 0)      /* exclude i = j from summation */
    for(i = 0; i < (*n); i++)
      {
        alpha = y[i] * y[i];   
        F77_CALL(dsyr)(&uplo, p, &alpha, gk + i, n, Mmat, p);
        wsum += (alpha * sc1);  /* y[i] * y[i] / (2 * h) */
      }

  for(i = 0; i < (*n); i++)
    for(j = i + 1; j < (*n); j++)
      {
        /* calculate xij = x[i, ] - x[j, ] */
        F77_CALL(dcopy)(p, x + i, n, xij, &one); /* DCOPY - copy x to y */
        F77_CALL(daxpy)(p, &b, x + j, n, xij, &one); /* DAXPY - replace y by da*x + y */

        wij = F77_CALL(ddot)(p, xij, &one, xij, &one); /* DDOT - inner product of x and y */
        alpha = y[i] * y[j] * exp(sc2 * wij);
        wsum += (alpha / (*h));

        /* gkxij = gk[i, ] - xij / (2 * h) */
        F77_CALL(dscal)(p, &sc1, xij, &one); /* DSCAL - scale a one-dimensional array */
        F77_CALL(dcopy)(p, gk + i, n, gkxij, &one); /* DCOPY - copy x to y */
        F77_CALL(daxpy)(p, &b, xij, &one, gkxij, &one); /* DAXPY - replace y by da*x + y */
        /* xij = gk[j, ] + xij / (2 * h) */
        F77_CALL(daxpy)(p, &a, gk + j, n, xij, &one); /* DAXPY - replace y by da*x + y */

        F77_CALL(dsyr2)(&uplo, p, &alpha, gkxij, &one, xij, &one, Mmat, p); /* DSYR2 - perform the symmetric rank 2 operation */
      }

  for(i = 0; i < pp; i+=((*p)+1)) Mmat[i] += wsum;
  alpha = 1.0 / (double)( (*n) * ((*n) - (*out)) );
  F77_CALL(dscal)(&pp, &alpha, Mmat, &one); /* DSCAL - scale a one-dimensional array */
}


/**********************************************************************
 * Calculate the candidate matrix for the Central Subspace.
 * Weight function W(x, u): Gaussian function with covariance h * I_p.
 * Weight function w(y, v): Gaussian function with variance h.
 *
 * Input:  
 *    x  --  n-by-p matrix  (predictors)
 *    y  --  n-by-1 vector  (response)
 *    gk --  n-by-p matrix  (derivative of log f(x))
 *    hx --  scalar         (sigma^2 for x)
 *    hy --  scalr          (sigma^2 for y)
 *    out--  1 or 0         (1 = leave-i-out, 0 = not)
 *
 * Output:
 *    Mmat -- p-by-p symmetric matrix  (candidate matrix)
 *            notice: only the lower triangle of Mmat are calculated. 
 **********************************************************************/


void FM_pdf(const double* x, const double* gk, const double* y, 
             const double* hx, const double* hy,
             const int* n, const int* p, const int* out, double* Mmat)
{
  const char uplo = 'L', trans = 'T';
  const double sc1x = (*hx), sc2x = -0.5*(*hx);
  const double sc2y = -0.5*(*hy);
  const double a = 1.0, b = -1.0;
  const int pp = (*p) * (*p), one = 1;
  double alpha, wij, wsum; 
  double xij[*p], gkxij[*p];
  int i, j;

  for(i = 0; i < pp; i++)  Mmat[i] = 0.0;

  wsum = 0.0;
  if((*out) == 0)        /* exclude i = j from summation */
  {
    F77_CALL(dsyrk)(&uplo, &trans, p, n, &a, gk, n, &a, Mmat, p); /* DSYRK - perform one of the symmetric rank k operations */
    wsum = ( sc1x * ((double)(*n)) );
  }

  for(i = 0; i < (*n); i++)
    for(j = i+1; j < (*n); j++)
      {
        /* calculate xij = x[i, ] - x[j, ] */
        F77_CALL(dcopy)(p, x+i, n, xij, &one);
        F77_CALL(daxpy)(p, &b, x+j, n, xij, &one); 

        wij = F77_CALL(ddot)(p, xij, &one, xij, &one);
        alpha = y[i] - y[j];
        alpha = exp(alpha * alpha * sc2y + wij * sc2x);
        wsum += (alpha / (*hx));

        /* gkxij = gk[i, ] - xij / (2 * hx) */
        F77_CALL(dscal)(p, &sc1x, xij, &one);
        F77_CALL(dcopy)(p, gk + i, n, gkxij, &one);
        F77_CALL(daxpy)(p, &b, xij, &one, gkxij, &one);
        /* xij = gk[j, ] + xij / (2 * hx) */
        F77_CALL(daxpy)(p, &a, gk + j, n, xij, &one);

        F77_CALL(dsyr2)(&uplo, p, &alpha, gkxij, &one, xij, &one, Mmat, p);
      }

  for(i = 0; i < pp; i+=((*p)+1)) Mmat[i] += wsum;
  alpha = 1.0 / (double)( (*n) * ((*n) - (*out)) );
  F77_CALL(dscal)(&pp, &alpha, Mmat, &one);
}


/**********************************************************************
 * Calculate the candidate matrix for the Central Mean Subspace.
 * Weight function W(x, u): Gaussian function with covariance h * I_p.
 * Assume the density of x is multivariate normal N(0, I_p). 
 * Therefore, gk = -x
 *
 * Input:  
 *    x  --  n-by-p matrix  (predictors)
 *    y  --  n-by-1 vector  (response)
 *    h  --  scalar         (sigma^2)
 *    out--  1 or 0         (1 = leave-i-out, 0 = not)
 *
 * Output:
 *    Mmat -- p-by-p symmetric matrix  (candidate matrix)
 *            notice: only the lower triangle of Mmat are calculated. 
 **********************************************************************/


void FM_mean_norm(const double* x, const double* y, const double* h,
                   const int* n, const int* p, const int* out, 
                   double* Mmat)
{
  const char uplo = 'L';
  const double sc1 = -(*h)/2;
  const double b = (1.0 + (*h))*2 * (*h);
  const double a = b + 1.0;
  const double oned = -1.0;
  const int pp = (*p) * (*p), one = 1;
  double aij, alpha, Asum;
  double Arow[*n], xij[*p];
  int i, j;

  for(i = 0; i < pp; i++) Mmat[i] = 0.0;
  if((*out) == 0) for(i = 0; i < (*n); i++) Arow[i] = (y[i] * y[i]);
  else for(i = 0; i < (*n); i++) Arow[i] = 0.0;

  for(i = 0; i < (*n); i++)
    for(j = i+1; j < (*n); j++)
      {
        /* calculate xij = x[i, ] - x[j, ] */
        F77_CALL(dcopy)(p, x + i, n, xij, &one);
        F77_CALL(daxpy)(p, &oned, x + j, n, xij, &one);

        aij = F77_CALL(ddot)(p, xij, &one, xij, &one);
        aij = exp(sc1 * aij) * (y[i] * y[j]);
        Arow[i] += aij; Arow[j] += aij;

        alpha = aij * a;
        F77_CALL(dsyr2)(&uplo, p, &alpha, x+i, n, x+j, n, Mmat, p);
      }

  if((*out) == 0) for(i = 0; i < (*n); i++)
    {
      alpha = a * y[i] * y[i] - Arow[i] * b; 
      F77_CALL(dsyr)(&uplo, p, &alpha, x + i, n, Mmat, p);
    }
  else for(i = 0; i < (*n); i++)
    {
      alpha = -b * Arow[i];
      F77_CALL(dsyr)(&uplo, p, &alpha, x + i, n, Mmat, p);
    }

  Asum = 0.0;
  for(i = 0; i < (*n); i++) Asum += Arow[i];
  alpha = Asum * (*h);
  for(i = 0; i < pp; i+=((*p)+1)) Mmat[i] += alpha;
  alpha = 1.0 / (double)( (*n) * ((*n) - (*out)) );
  F77_CALL(dscal)(&pp, &alpha, Mmat, &one);
}


/**********************************************************************
 * Calculate the candidate matrix for the Central Subspace.
 * Weight function W(x, u): Gaussian function with covariance h * I_p.
 * Weight function w(y, v): Gaussian function with variance h.
 * Assume the density of x is multivariate normal N(0, I_p). 
 * Therefore, gk = -x
 *
 * Input:  
 *    x  --  n-by-p matrix  (predictors)
 *    y  --  n-by-1 vector  (response)
 *    hx --  scalar         (sigma^2 for x)
 *    hy --  scalar         (sigma^2 for y)
 *    out--  1 or 0         (1 = leave-i-out, 0 = not)
 *
 * Output:
 *    Mmat -- p-by-p symmetric matrix  (candidate matrix)
 *            notice: only the lower triangle of Mmat are calculated. 
 **********************************************************************/

void FM_pdf_norm(const double* x, const double* y, const double* hx,
                  const double* hy, const int* n, const int* p, 
                  const int* out, double* Mmat)
{
  const char uplo = 'L';
  const double sc1 = -0.5*(*hx);
  const double scy = -0.5*(*hy);
  const double b = (1.0 + (*hx)) *2* (*hx);
  const double a = b + 1.0;
  const double oned = -1.0;
  const int pp = (*p) * (*p), one = 1;
  double aij, alpha, Asum;
  double Arow[*n], xij[*p];
  int i, j;

  for(i = 0; i < pp; i++) Mmat[i] = 0.0;
  if((*out) == 0) for(i = 0; i < (*n); i++) Arow[i] = 1.0;
  else for(i = 0; i < (*n); i++) Arow[i] = 0.0;

  for(i = 0; i < (*n); i++)
    for(j = i+1; j < (*n); j++)
      {
        /* calculate xij = x[i, ] - x[j, ] */
        F77_CALL(dcopy)(p, x + i, n, xij, &one);
        F77_CALL(daxpy)(p, &oned, x + j, n, xij, &one);

        aij = F77_CALL(ddot)(p, xij, &one, xij, &one);
        alpha = y[i] - y[j];
        aij = exp(sc1 * aij + scy * alpha * alpha);
        Arow[i] += aij; Arow[j] += aij;

        alpha = aij * a;
        F77_CALL(dsyr2)(&uplo, p, &alpha, x+i, n, x+j, n, Mmat, p);
      }

  if((*out) == 0) for(i = 0; i < (*n); i++)
    {
      alpha = a - Arow[i] * b; 
      F77_CALL(dsyr)(&uplo, p, &alpha, x+i, n, Mmat, p);
    }
  else for(i = 0; i < (*n); i++)
    {
      alpha = -b * Arow[i];
      F77_CALL(dsyr)(&uplo, p, &alpha, x+i, n, Mmat, p);
    }

  Asum = 0.0;
  for(i = 0; i < (*n); i++) Asum += Arow[i];
  alpha = Asum*(*hx) ;
  for(i = 0; i < pp; i+=((*p)+1)) Mmat[i] += alpha;
  alpha = 1.0 / (double)( (*n) * ((*n) - (*out)) );
  F77_CALL(dscal)(&pp, &alpha, Mmat, &one);
}

/**********************************************************************
 * Estimate density and derivative of log density using kernel method.
 * Use Gaussian kernel function.
 *
 * Input:
 *    x    --  n-by-p matrix  (data)
 *    h    --  scalar (bandwidth)
 *    out  --  1 or 0 (1: exclue i = j from summation, 0: not)
 *
 * Output:
 *    den   --  n-by-1 vector  (estimated density at x)
 *    dlogf --  n-by-p matrix  (estimated d-log-den at x)
 **********************************************************************/

void Fdlogden1(const double* x, const int* n, const int* p,
              const double* h, const int* out, 
              double* den, double* dlogf)
{
  const int np = (*n) * (*p); 
  const double sc = R_pow_di(M_1_SQRT_2PI / (*h), *p) / (double)((*n) - (*out));
  double a, b, wi0;
  double xixi[*n];
  int i, j;

  for(i = 0; i < (*n); i++)
    /* DDOT - inner product of x and y */
    xixi[i] = F77_CALL(ddot)(p, x + i, n, x + i, n);

  if((*out) == 0)/*in C all 0, false and none are same*/
    {
      for(i = 0; i < (*n); i++) den[i] = 1.0;
      for(i = 0; i < np; i++) dlogf[i] = x[i];
    }
  else
    {
      for(i = 0; i < (*n); i++) den[i] = 0.0;
      for(i = 0; i < np; i++) dlogf[i] = 0.0;
    }

  a = -0.5 / (*h) / (*h);
  for(i = 0; i < (*n); i++)
    for(j = i + 1; j < (*n); j++)
      {
        wi0 = F77_CALL(ddot)(p, x + i, n, x + j, n);
        wi0 = exp(a * (xixi[i] + xixi[j] - 2.0 * wi0));
        den[i] += wi0;
        den[j] += wi0;
        /* DAXPY - replace y by da*x + y */
        F77_CALL(daxpy)(p, &wi0, x + i, n, dlogf + j, n);
        F77_CALL(daxpy)(p, &wi0, x + j, n, dlogf + i, n);
      }

  b = 1.0 / (*h) / (*h);
  for(i = 0; i < (*n); i++)
    {  
      a = b / den[i];
    /* DSCAL - scale a one-dimensional array */
      F77_CALL(dscal)(p, &a, dlogf + i, n);
      den[i] *= sc;
    }

  b = -1.0 / (*h) / (*h); 
  for(i = 0; i < np; i++) dlogf[i] += (x[i] * b);
}


/**********************************************************************
 * Estimate density and derivative of log density using kernel method.
 * Use Gaussian kernel function. 
 * X follows elliptically contoured distribution with mean 0.
 *
 * Input:
 *    x    --  n-by-p matrix  (data)
 *    h    --  scalar (bandwidth)
 *    out  --  1 or 0 (1: exclue i = j from summation, 0: not)
 *
 * Output:
 *    den    --  n-by-1 vector  (estimated density at r)
 *    rdlogf --  n-by-1 vector  (estimated density at r)
 *    dlogf  --  n-by-p matrix  (estimated d-log-den at x)
 **********************************************************************/

void Fdlogden3(const double* x, const int* n, const int* p,
              const double* h, const int* out,
              double* den, double* rdlogf, double* dlogf)
{
  const int np = (*n) * (*p);
  const double sc = M_1_SQRT_2PI / (*h) / (double)((*n) - (*out));
  double a, b, wi0;
  double r[*n];
  int i, j;

  for(i = 0; i < (*n); i++)
    r[i] = sqrt(F77_CALL(ddot)(p, x + i, n, x + i, n));

  if((*out) == 0) for(i = 0; i < (*n); i++) 
    { 
      den[i] = 1.0; 
      rdlogf[i] = r[i]; 
    }
  else for(i = 0; i < (*n); i++) 
    { 
      den[i] = 0.0; 
      rdlogf[i] = 0.0; 
    }

  a = -0.5 / (*h) / (*h);
  for(i = 0; i < (*n); i++)
    for(j = i + 1; j < (*n); j++)
      {
        wi0 = r[i] - r[j];
        wi0 = exp(a * wi0 * wi0);
        den[i] += wi0;
        den[j] += wi0;

        rdlogf[i] += (wi0 * r[j]);
        rdlogf[j] += (wi0 * r[i]);
      }

  for(i = 0; i < np; i++) dlogf[i] = x[i];

  a = 1.0 / (*h) / (*h);
  b = (double)(*p) - 1.0;
  for(i = 0; i < (*n); i++)
    {
      rdlogf[i] = a * (rdlogf[i] / den[i] - r[i]);       
      den[i] *= sc;
      wi0 = (rdlogf[i] - b / r[i]) / r[i];
      F77_CALL(dscal)(p, &wi0, dlogf + i, n);
    }
}



/**********************************************************************
 *   END OF THE FILE
 **********************************************************************/
