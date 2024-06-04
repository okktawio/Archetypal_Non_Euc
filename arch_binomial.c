#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>


#define __iter__ 200

#define infinity 1.7976931348623157e308

double logit(double p) {
  if (p <= 0) p = 1e-12;
  if (p >= 1) p = 1 - 1e-12;
  return log(p / (1 - p));
}

double expit(double x) {
  return 1.0 / (1.0 + exp(-x));
}

void transfer_m(double **M, double *A, int x, int y)
{
  for(int i = 0; i < x; ++i)
    M[i] = &A[i * y];
}

double binomial_like(double k, double n, double p) {
    // Calculate the likelihood function
  p = fmax(fmin(p, 1-1e-9),1e-9);
  return -2 * lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1) + k * log(p) + (n - k) * log(1 - p);
}

double distancia_i(double P, double X, double N)
{
  return binomial_like(X, N, P / N);
}

double distancia_binomial(double *X, double *N, double *P, int size) {
  double da = 0, di = 0;
  for (int i = 0; i < size; i++)
    {
      di = -1;
      if(X[i] >= 0)
	{
	  di = binomial_like(X[i], N[i], P[i] / N[i]);
	  da += di;
	}
    }
  return da;
}

double distanciadi(double * A1, int * A2, int size) {
  double da = 0;
  for (int i = 0; i < size; i++)
    if(A1[i] > 0 || A2[i] > 0)
      da += pow(A1[i] - A2[i], 2);
  return da;
}

double distanciadd(double * A1, double * A2, int size) {
  double da = 0;
  for (int i = 0; i < size; i++)
    {
      if(A1[i] >= 0 && A2[i] >= 0)
	da += pow(A1[i] - A2[i], 2);
    }
  return da;
}


void estima_Si(double **XC, double Se[], double E[], int m, int n)
{
  for(uint k = 0; k < n; ++k)
    {
      E[k]  = 0;
      for(uint j = 0; j < m; ++j)
	E[k] += Se[j] * XC[j][k];
    }
}

double error_t(double **XC, double **S, double **A, double *N, int m, int n, int x) {
  double et = 0.0;
  double E[n];
  for (int k = 0; k < x; k++)
    {
      estima_Si(XC, S[k], E, m, n);
      double ei = distancia_binomial(A[k], N, E, n);
      et += ei;
    }
  return et;
}

double distancia_media(double **XC, double *A, int m, int n) {
  double SS = 0.0;
  for(int j = 0; j < m; ++j)
      SS += distanciadd(XC[j], A, n);
  return SS / m;
}

int ya_esta(int ix[], int ii, int x)
{
  for(int i = 0; i < x; ++i)
    if(ix[i] == ii) return 1;
  return 0;
}

int bucle_lejano(double **XC, double **A, int *ix, int m, int n, int x, int *validos) {
  double dmax = 0.0, di;
  int ii = 0;
  for(int i = 0; i < x; ++i)
    {
      if(validos[i] == 0)
	{
	  di = distancia_media(XC, A[i], m, n);
	  if(di > dmax && !(ya_esta(ix, i, m)))
	    {
	      ii = i;
	      dmax = di;
	    }
	}
    }
  return ii;
}

void update_XC(double **XC, double **A, int *ix, int m, int n, int ii) {
  // Copy XC[1:,:] to XC[:-1,:]
  if(m > 1)
    {
      for(int i = 1; i < m; ++i)
	memcpy(XC[i - 1], XC[i], sizeof(double) * n);
      memcpy(ix, ix + 1, sizeof(uint) * (m - 1));
    }
  // Copy data[i,:] to XC[-1,:]
  memcpy(XC[m - 1], A[ii], sizeof(double) * n);
  ix[m-1] = ii;
}


void init_XC(double **XC, double **A, int *ix, int m, int n, int x, int * validos)
{
  int ii;
  for(int i = 0; i < m; ++i)
    {
      do
	{
	  ii = (int) (drand48() * x);
	}
      while(validos[ii] || ya_esta(ix, ii, m));
      memcpy(XC[i], A[ii], sizeof(double) * n);
      ix[i] = ii;
    }
}

void furthest_sum(double **XC, double **A, int *ix, int m, int n, int x, int * validos)
{
  for(int i = 0; i < m; ++i) ix[i] = -99999;
  init_XC(XC, A, ix, m, n, x, validos);
  //if(drand48() < 0.01)
    {
      for(int i = 0; i < m; ++i)
	{
	  int ii = bucle_lejano(XC, A, ix, m, n, x, validos);
	  update_XC(XC, A, ix, m, n, ii);
	}
    }
    //printf("\nix: ");
    //for(int i = 0; i < m; ++i) printf("%d ",ix[i]);
    //printf("\n");
}

void init_C(double **C, int *ix, int m, int x)
{
  for(int i = 0; i < m; ++i)
    for(int j = 0; j < x; ++j)
      C[i][j] = (double) 0.1 / (x - 1);

  //printf("IX ");
  for(int i = 0; i < m; ++i)
    {
      //printf("%d ", ix[i]);
      C[i][ix[i]] = (double) 0.9;
    }
  //printf("\n");
}
  

void retorna_pesos(double **XC, int m, int n, double *Ai, double *pesos) {
  for (int i = 0; i < m; i++)
    pesos[i] = 1.0 / m;
}

double inicia_S(double **XC, double **A, double *N, double **S, int m, int n, int x) {
  int i;
  for (i = 0; i < x; ++i)
    retorna_pesos(XC, m, n, A[i], S[i]);
  
  return error_t(XC, S, A, N, m, n, x);
}

double error_g(double Se[], double **XC, double *Ai, double * N, int m, int n)
{
  double E[n];
  estima_Si(XC, Se, E, m, n);
  return distancia_binomial(Ai, N, E, n);
}

double desciende_s(double delta, double *gr, double **XC, double *Se, double *Si, double *Ai, double *N,
		 int m, int n)
{
  memcpy(Si, Se, m * sizeof(double));

  double sum = 0;
  for(int i = 0; i < m; ++i)
    {
      Si[i] = expit(logit(Si[i]) - delta * gr[i]);
      sum += Si[i];
    }
  if(sum > 0)
    {
      for(int i = 0; i < m; ++i)
	Si[i] /= sum;
      return error_g(Si, XC, Ai, N, m, n);
    }
  else return 1e38;
}

double golden_linesearchs(double delta, double *gr, double **XC, double *Se, double *Si, double *Ai, double *N,
			  int m, int n)
{
  double gor = (sqrt(5) + 1) / 2;
  double a = 0, b = delta, tol=1e-8;
  double c, d, fb;
  double e0 = error_g(Se, XC, Ai, N, m, n);
  fb = desciende_s(b, gr, XC, Se, Si, Ai, N, m, n);

   while(fb < e0 && b < 100.0)
    {
      b *= gor;
      fb = desciende_s(b, gr, XC, Se, Si, Ai, N, m, n);
     }

  while (fabs(b - a) > tol)
    {
      c = b - (b - a) / gor;
      d = a + (b - a) / gor;
      double fc = desciende_s(c, gr, XC, Se, Si, Ai, N, m, n);
      double fd = desciende_s(d, gr, XC, Se, Si, Ai, N, m, n);
       if (fc < fd) b = d;
      else         a = c;
    }

  delta = (b + a) / 2;
  double ei = desciende_s(delta, gr, XC, Se, Si, Ai, N, m, n);
   
  if (ei < e0)
    {
      memcpy(Se, Si, m * sizeof(double));
      e0 = ei;
    }
  return e0;
}

double gradient_s(double Se[], double gr[], double **XC, double *Ai, double *N, int m, int n)
{
  double df = 1e-4;
  double norma = 0;
  double Si[m];
  double e0 = error_g(Se, XC, Ai, N, m, n);
   
  for(int i = 0; i < m; ++i)
    {
      memcpy(Si, Se, m * sizeof(double));
      Si[i] = expit(logit(Si[i]) + df);
      for(int j = 0; j < m; ++j)
	Si[j] /= 1 + df;
      double ei = error_g(Si, XC, Ai, N, m, n);
      gr[i] = (ei - e0) / df;

      memcpy(Si, Se, m * sizeof(double));
      Si[i] = expit(logit(Si[i]) - df);
      for(int j = 0; j < m; ++j)
	Si[j] /= 1 - df;
      ei = error_g(Si, XC, Ai, N, m, n);
      gr[i] += (e0 - ei) / df;
      gr[i] /= 2.0;
      norma += pow(gr[i], 2);
    }
  return norma;
}

double gd_s(double **XC, double *Ai, double *N, double *Se, int m, int n)
{
  double e0, e4;
  double Si[m];
  double gr[m];
  double norma = 1e9;

  e4 = error_g(Se, XC, Ai, N, m, n);
  e0 = e4 + e4;
  int itera = 0;
  double delta = 1.0;

  while((norma > 1e-6) && (itera < 1000) && (fabs(e0 - e4) > 1e-8))
    {
      e0 = e4;
      norma = gradient_s(Se, gr, XC, Ai, N, m, n);
      //delta = 1e-4; //fmax(fmin(0.5/maxgrad, 1), 1e-6);
      e4 = golden_linesearchs(delta, gr, XC, Se, Si, Ai, N, m, n);
      ++itera;
    }
  return e4;
}

double error_ij(double **E, double Xo, double Xd, double **Se, double **A, double *N, int m, int n, int x, int i, int j)
{
  double dx, ds, SS = 0, so = 0;
  for (int k = 0; k < x; k++)
    {
      if(A[k][i] >= 0 && E[k][i] >= 0)
	{
	  double Eo = E[k][i];
	  //printf("Error i %d j %d k %3d Eo %f Aij %f Ni %f ", i,j,k,Eo,A[k][i], N[i]);
	  so =  distancia_i(Eo, A[k][i], N[i]);
	  dx = Xd - Xo;
	  ds = Se[k][j] * dx;
	  Eo -= ds;
	  SS += distancia_i(Eo, A[k][i], N[i]) - so;
	}
    }
  return SS;
}


double update_S(double **XC, double **A, double *N, double **S, int m, int n, int x)
{
  int i;
  double e0, et = 0, E[n];
  if(m <= 1)
    {
      for (i = 0; i < x; ++i)
	{
	  S[i][0] = 1.0;
	  estima_Si(XC, S[i], E, m, n);
	  et += distancia_binomial(A[i], N, E, n);
	}
      return et;
    }
  else
    {
      for (i = 0; i < x; ++i)
	{
	  e0 = gd_s(XC, A[i], N, S[i], m, n);
	  et += e0;
	}
      return et;
    }
}

void estima_E(double **E, double **XC, double **S, int m, int n, int x) {
  for (int k = 0; k < x; k++)
      estima_Si(XC, S[k], E[k], m, n);
} 

void printXC(double **XC, int m, int n)
{
  printf("\nXC\n");
  for(int i = 0; i < m; ++i)
    {
      for(int j = 0; j < n; ++j)
	{
	  printf("%12.6f ", XC[i][j]);
	}
      printf("\n");
    }
  printf("\n");
}


double gradient_XCS(double **GXC, double **XC, double **A, double *N, double **S, int m, int n, int x)
{
  double delta = 1e-4;
  //double eo    = error_t(XC, S, A, N, m, n, x);
  double ei, ej;
  double norma = 0;
  
  double *Ea = (double *) malloc(n * x * sizeof(double));
  double **E = (double **) malloc(x * sizeof(double *));
  

  transfer_m(E, Ea, x, n);
  estima_E(E, XC, S, m, n, x);

  for(uint i = 0; i < n; ++i)
    for(uint j = 0; j < m; ++j)
      {
	double Xo = XC[j][i];
	XC[j][i] += delta;
	ei = error_ij(E, Xo, XC[j][i], S, A, N, m, n, x, i, j);
	GXC[j][i] = ei / delta;
	
	XC[j][i] = Xo - delta;
	ej = error_ij(E, Xo, XC[j][i], S, A, N, m, n, x, i, j);
	GXC[j][i] -= ej / delta;
	
	GXC[j][i] /= 2.0;
	XC[j][i] = Xo;
	norma += pow(GXC[j][i], 2);
	//printf("i %d j %d eo %f ei %f ej %f\n", i, j, eo, ei, ej);
      }
  free(Ea);
  free(E);
  return norma;
}

double desciendeXC(double delta, double **GX, double **Xd, double **XC, double **S, double **A,
		 double *N, int m, int n, int x)
{
  memcpy(Xd[0], XC[0], n * m * sizeof(double));
  for(uint i = 0; i < n; ++i)
    for(uint j = 0; j < m; ++j)
      Xd[j][i] = fmax(fmin(Xd[j][i] - delta * GX[j][i], N[i]), 0);
  return error_t(Xd, S, A, N, m, n, x);
}

double golden_linesearchxc(double delta, double e0, double **GX, double **Xd, double **XC, double **S, double **A,
			   double *N, int m, int n, int x)
{
  double gr = (sqrt(5) + 1) / 2;
  double a = 0, b = delta, tol=1e-8;
  double c, d, fb;
  /*
    """Golden-section search
    to find the minimum of f on [a,b]
  */
  fb = desciendeXC(b, GX, Xd, XC, S, A, N, m, n, x);
  while(fb < e0 && b < 1)
    {
      b *= gr;
      fb = desciendeXC(b, GX, Xd, XC, S, A, N, m, n, x);
    }

  while (fabs(b - a) > tol)
    {
      c = b - (b - a) / gr;
      d = a + (b - a) / gr;
      double fc = desciendeXC(c, GX, Xd, XC, S, A, N, m, n, x);
      double fd = desciendeXC(d, GX, Xd, XC, S, A, N, m, n, x);
      if (fc < fd) b = d;
      else         a = c;
    }
  delta = (b + a) / 2;
  double ei = desciendeXC(delta, GX, Xd, XC, S, A, N, m, n, x);
  if (ei < e0)
    {
      memcpy(XC[0], Xd[0], m * n * sizeof(double));
      e0 = ei;
    }
  return e0;
}
  
void gd_XCS(double **XC, double **A, double *N, double **S, int m, int n, int x)
{
  double *GXC = (double *) malloc(n * m * sizeof(double));
  double **GX = (double **) malloc(m * sizeof(double *));
  double *XCd = (double *) malloc(n * m * sizeof(double));
  double **Xd = (double **) malloc(m * sizeof(double *));
  
  transfer_m(GX, GXC, m, n);
  transfer_m(Xd, XCd, m, n);
  
  double e0 = error_t(XC, S, A, N, m, n, x), ei;
  int itera = 0;
  double delta = 1e-6;
  int acepta = 1;
  double norma = 1e9;
  while(norma > 1e-4 && itera < 100 && acepta == 1)
    {
      acepta = 0;
      norma = gradient_XCS(GX, XC, A, N, S, m, n, x);
      //printf("GR XC\n");
      //printXC(GX, m, n);
      ei = golden_linesearchxc(delta, e0, GX, Xd, XC, S, A, N, m, n, x);
      if(ei < e0)
	{
	  e0 = ei;
	  acepta = 1;
	}
      ++itera;
    }
  free(GXC);
  free(GX);
  free(XCd);
  free(Xd);
}


double calcula_XCS(double **XC, double **A, double *N, double **S, int m, int n, int x)
{
  gd_XCS(XC, A, N, S, m, n, x);
  return error_t(XC, S, A, N, m, n, x);
}

void estima_XCi(double **A, double Ce[], double E[], int n, int x)
{
  for(int k = 0; k < n; ++k)
    {
      E[k]  = 0;
      for(int j = 0; j < x; ++j)
	  E[k] += Ce[j] * A[j][k];
    }
}

double error_c(double *Ce, double **XC, double **A, int m, int n, int x, int ix)
{
  double E[n];
  estima_XCi(A, Ce, E, n, x);
  return distanciadd(E, XC[ix], n);
}

void actualiza_XCi(double **A, double Ce, double Ci, double E[], int n, int xj)
{
  for(uint k = 0; k < n; ++k)
    //primero resto la estimacion anterior
    E[k] -= Ce * A[xj][k];
  for(uint k = 0; k < n; ++k)
    //luego le sumo la estimacion actual
    E[k] -= Ce * A[xj][k];
}


double error_cg(double Ce, double Ci, double E[], double **XC, double **A, int m, int n, int xj, int ix)
{
  actualiza_XCi(A, Ce, Ci, E, n, xj);
  return distanciadd(E, XC[ix], n);
}


double desciende_c(double delta, double *gr, double **XC, double *Ce, double *Ci, double **A,
		   int m, int n, int x, int ix)
{
  memcpy(Ci, Ce, x * sizeof(double));
  double sum = 0;
  for(int i = 0; i < x; ++i)
    {
      Ci[i] = expit(logit(Ci[i]) - delta * gr[i]); //fmax(fmin(Ci[i] - delta * gr[i], 1), 0);
      sum += Ci[i];
    }
  for(int i = 0; i < x; ++i)
    Ci[i] /= sum;
  
  return error_c(Ci, XC, A, m, n, x, ix);
}

double golden_linesearchx(double delta, double *gr, double **XC, double *Ce, double *Ci, double **A,
			      int m, int n, int x, int ix)
{
  double gor = (sqrt(5) + 1) / 2;
  double a = 0, b = delta, tol=1e-12;
  double c, d, fb;
  double e0 = error_c(Ce, XC, A, m, n, x, ix);
  fb = desciende_c(b, gr, XC, Ce, Ci, A, m, n, x, ix);
  while(fb < e0 && b < 0.1)
    {
      b *= gor;
      fb = desciende_c(b, gr, XC, Ce, Ci, A, m, n, x, ix);
    }

  while (fabs(b - a) > tol)
    {
      c = b - (b - a) / gor;
      d = a + (b - a) / gor;
      double fc = desciende_c(c, gr, XC, Ce, Ci, A, m, n, x, ix);
      double fd = desciende_c(d, gr, XC, Ce, Ci, A, m, n, x, ix);
      if (fc < fd) b = d;
      else         a = c;
    }
  delta = (b + a) / 2;
  double ei = desciende_c(delta, gr, XC, Ce, Ci, A, m, n, x, ix);
  if (ei < e0)
    {
      memcpy(Ce, Ci, x * sizeof(double));
      e0 = ei;
    }
  return e0;
}

double gradient_c(double *Ce, double gr[], double **XC, double **Ai, int m, int n, int x, int ix)
{
  double df = 1e-4, dg = 0;
  double norma = 0;
  double Ci[x];
  double E[n];
  estima_XCi(Ai, Ce, E, n, x);
  double e0 = distanciadd(E, XC[ix], n);
  
  for(int i = 0; i < x; ++i)
    {
      memcpy(Ci, Ce, x * sizeof(double));
      dg = expit(logit(Ci[i]) + df) - Ci[i];
      Ci[i] += dg;
      //Ci[i] = fmin(Ci[i], 1);
      for(int j = 0; j < x; ++j)
      	Ci[j] /= 1 + dg;
      double ei = error_c(Ci, XC, Ai, m, n, x, ix);
      gr[i] = (ei - e0) / df; //error_c(Ce[i], Ci[i], E, XC, Ai, m, n, i, ix) - e0) / df;
      memcpy(Ci, Ce, x * sizeof(double));
      dg = expit(logit(Ci[i]) - df) - Ci[i];
      Ci[i] += dg;
      //Ci[i] = fmax(Ci[i], 0);
      for(int j = 0; j < x; ++j)
	Ci[j] /= 1 + dg;
      double ej = error_c(Ci, XC, Ai, m, n, x, ix);
      gr[i] += (e0 - ej) / df;  //(e0 - error_cg(Ce[i], Ci[i], E, XC, Ai, m, n, i, ix)) / df;
      gr[i] /= 2.0;
      norma += pow(gr[i], 2);
    }
  return norma;
}


double gd_c(double **XC, double **Ai, double *Ce, int m, int n, int x, int ix)
{
  double e0, e4;
  double gr[x];
  double norma = 1e9;
  double Ci[x];

  e4 = error_c(Ce, XC, Ai, m, n, x, ix);
  
  e0 = e4 + 1;
  int itera = 0;
  double delta = 0.1;

  //printf("\nix %d ei %f \n", ix, e4);
  while(norma > 1e-4 && itera < 1000 && fabs(e0 - e4) > 1e-8)
    {
      e0 = e4;
      norma = gradient_c(Ce, gr, XC, Ai, m, n, x, ix);
      e4 = golden_linesearchx(delta, gr, XC, Ce, Ci, Ai, m, n, x, ix);
      //printf("norma %f ei %f \n", norma, e4);
      ++itera;
    }
  //printf("ef %f \n", e4);
  return e4;
}


double update_C(double **XC, double **A, double **C, int m, int n, int x)
{
  int i;
  double E[n];
  double e0 = 0;
  for (i = 0; i < m; ++i)
    {
      estima_XCi(A, C[i], E, n, x);
      e0 += gd_c(XC, A, C[i], m, n, x, i);
      estima_XCi(A, C[i], XC[i], n, x);
    }
  return e0;
}

double inf_index(double like, int n, int k, double* results) {
  like = -like / 2;
  double bic = -like + k * log(n);
  double aicc = -2 * like + 2 * k + 2 * k * (k + 1) / (n - k - 1);
  
  results[0] = like;
  results[1] = bic;
  results[2] = aicc;
  return like;
}


void check_validos(double **A, int *validos, int n, int x)
{
  for(int j = 0; j < x; ++j)
    {
      validos[j] = 0;
      for(int i = 0; i < n; ++i)
	{
	  validos[j] += ((int) A[j][i] < 0);
	  if(validos[j] > 0) break;
	}
    }
} 


double Archetypal(double *Ae, double *N, double *XCe, double *Ce, double *Se, int m, int n, int x, double *IE)
{
  double llx = -1e38, lli = 0, lls, llz,  aic = 1e9, bic = 1e9, ec = 0, ex = 1; //likelihood
  int validos[x];

  double *XC[m], *C[m], *S[x], *A[x];
  transfer_m(A,   Ae, x, n);
  transfer_m(XC, XCe, m, n);
  transfer_m(C,   Ce, m, x);
  transfer_m(S,   Se, x, m);
  srand48(time(NULL));

  check_validos(A, validos, n, x);
  {
    int ix[m];
    double *XCr = (double *) malloc(n * m * sizeof(double));
    double *Sr  = (double *) malloc(x * m * sizeof(double));
    
    double **XCi = (double **) malloc(m * sizeof(double *));
    double **Si  = (double **) malloc(x * sizeof(double *));

    transfer_m(XCi, XCr, m, n);
    transfer_m(Si,   Sr, x, m);
    int i = 10;
    //initiate the archetypes
    while(i > 0)
      {
	furthest_sum(XCi, A, ix, m, n, x, validos);
	lli = inicia_S(XCi, A, N, Si, m, n, x); //ok
	lli = update_S(XCi, A, N, Si, m, n, x); //ok
	lli = inf_index(lli, x * n, m * x, IE);
	if(lli > llx)
	  {
	    memcpy(Se,  Sr,  x * m * sizeof(double));
	    memcpy(XCe, XCr, n * m * sizeof(double));
	    init_C(C, ix, m, x);
	    llx = lli;
	    aic = IE[2];
	    bic = IE[1];
	    i += 1;
	  }
	else --i;
	printf("\rNOC = %2d P %3d N %4d K %8d Error FS %4d  Llike %12.2f Llike max %12.2f AICc %12.2f BIC %12.2f",
	       m,  n, n * x, m * x, i, lli, llx, aic, bic);
	fflush(stdout);
      }
    free(XCr); free(XCi);
    free(Sr);  free(Si);
  }
  //printf("\nFin de Loop\n");
  //printXC(XC, m, n);
  int i = 0;
  lls = update_S(XC, A, N, S, m, n, x);
  lls = inf_index(lls, x * n, m * x, IE);
  lls = IE[0];
  llz = lls;
  printf("\rNOC = %2d P %3d N %4d K %8d Error Iteracion %4d Llike max %12.4f Llike maz %12.4f dflikex %12.6f dflikez %12.6f  AICc %12.2f BIC %12.2f", m,  n, n * x, m * x, i, llx, llz, lli - llx, lls - llz, IE[2], IE[1]);

  lli = calcula_XCS(XC, A, N, S, m, n, x); // ok
  lli = inf_index(lli, x * n, m * x, IE);
  lli = IE[0];
  llx = lli;
  printf("\rNOC = %2d P %3d N %4d K %8d Error Iteracion %4d Llike max %12.4f Llike maz %12.4f dflikex %12.6f dflikez %12.6f  AICc %12.2f BIC %12.2f", m,  n, n * x, m * x, i, llx, llz, lli - llx, lls - llz, IE[2], IE[1]);
  //printf("\nCalculado XCS \n");
  //printXC(XC, m, n);
  
  do
    {
      llx = lli;
      llz = lls;
      ex = ec;
      ec = update_C(XC, A, C, m, n, x);
      printf("\rNOC = %2d P %3d N %4d K %8d Error Iteracion %4d Llike max %12.4f Llike maz %12.4f dflikex %12.6f dflikez %12.6f  EC %12.6f AICc %12.2f BIC %12.2f", m,  n, n * x, m * x, i, lli, lls, lli - llx, lls - llz, ec, IE[2], IE[1]);
      //printf("\nCalculado C\n");
      //printXC(XC, m, n);
      lls = update_S(XC, A, N, S, m, n, x);
      lls = inf_index(lls, x * n, m * x, IE);
      printf("\rNOC = %2d P %3d N %4d K %8d Error Iteracion %4d Llike max %12.4f Llike maz %12.4f dflikex %12.6f dflikez %12.6f  EC %12.6f AICc %12.2f BIC %12.2f", m,  n, n * x, m * x, i, lli, lls, lli - llx, lls - llz, ec, IE[2], IE[1]);

      lli = calcula_XCS(XC, A, N, S, m, n, x);
      lli = inf_index(lli, x * n, m * x, IE);
      printf("\rNOC = %2d P %3d N %4d K %8d Error Iteracion %4d Llike max %12.4f Llike maz %12.4f dflikex %12.6f dflikez %12.6f  EC %12.6f dfEC %12.6f AICc %12.2f BIC %12.2f", m,  n, n * x, m * x, i, lli, lls, lli - llx, lls - llz, ec, ex - ec, IE[2], IE[1]);
      //printf("\nCalculado XCS \n");
      //printXC(XC, m, n);

      fflush(stdout);
      ++i;
    }
  while((fabs(lli - llx) > 1e-4 || fabs(lls - llz) > 1e-4 || fabs(ec - ex) > 1e-6) && i < 100000);
  printf("\n");fflush(stdout);
  return lli;
}


