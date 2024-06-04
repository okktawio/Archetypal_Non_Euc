#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>


#define __iter__ 400

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

double distancia_i(double A1, double A2, int distribution)
{
  double Di;
  if (distribution == 0)
    {
      Di = fabs(A1 - A2);
      Di = fmin(Di, 360 - Di);
    }
  else
    Di = A1 - A2;
  return pow(Di, 2);
}


double distancia_angular(double* A1, double* A2, int size) {
  double Di, da = 0;
  for (int i = 0; i < size; i++){
    Di = fabs(A1[i] - A2[i]);
    Di = fmin(Di, 360 - Di);
    da += pow(Di, 2);
  }
  return da;
}

double distancia_lineal(double* A1, double* A2, int size) {
  double da = 0;
  for (int i = 0; i < size; i++)
    da += pow(A1[i] - A2[i], 2);
  return da;
}

double distancia(double* A1, double* A2, int size, int distribution)
{
  switch(distribution){
  case 0:
    return distancia_angular(A1, A2, size); //angulos
  case 1:
    return distancia_lineal(A1, A2, size); //lineal
  default:
    return distancia_angular(A1, A2, size); //angulos
  }
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

double error_t(double **XC, double **S, double **A, int m, int n, int x, int distribution) {
  double SS = 0.0;
  double E[n];

  for (int k = 0; k < x; k++)
    {
      //printf("Estima Si %6d\n",k); fflush(stdout);
      estima_Si(XC, S[k], E, m, n);
      //printf("distancia\n"); fflush(stdout);
      SS += distancia(E, A[k], n, distribution);
    }
  return SS;
}

void estima_E(double **E, double **XC, double **S, int m, int n, int x) {
  for (int k = 0; k < x; k++)
    {
      //printf("Estima E %6d\n",k); fflush(stdout);
      estima_Si(XC, S[k], E[k], m, n);
    }
} 

double error_ij(double **E, double Xo, double Xd, double **Se, double **A, int m, int n, int x,
		int distribution, int i, int j)
{
  double dx, ds, SS = 0, so = 0;
  for (int k = 0; k < x; k++)
    {
      so = distancia_i(E[k][i], A[k][i], distribution);
      dx = Xo - Xd;
      ds = Se[k][j] * dx;
      E[k][i] -= ds;
      SS += distancia_i(E[k][i], A[k][i], distribution) - so;
      E[k][i] += ds;
    }
  return SS;
}

double distancia_media(double **XC, double *A, int m, int n, int distribution) {
  double SS = 0.0;
  for(uint j = 0; j < m; ++j)
    SS += distancia(XC[j], A, n, distribution);
  return SS / m;
}


int ya_esta(int ix[], int ii, int x)
{
  for(int i = 0; i < x; ++i)
    if(ix[i] == ii) return 1;
  return 0;
}


int bucle_lejano(double **XC, double **A, int *ix, int m, int n, int x, int distribution) {
  double dmax = 0.0, di;
  int ii = 0;
  for(int i = 0; i < x; ++i)
    {
      di = distancia_media(XC, A[i], m, n, distribution);
      if(di > dmax && !(ya_esta(ix, i, m)))
	{
	  ii = i;
	  dmax = di;
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

void append_XC(double **XC, double **A, int *ix, int m, int n, int ii) {
  memcpy(XC[m - 1], A[ii], sizeof(double) * n);
  ix[m-1] = ii;
}



void init_XC(double **XC, double **A, int *ix, int m, int n, int x)
{
  int ii;
  for(int i = 0; i < m; ++i)
    {
      do
	{
	  ii = (int) (drand48() * x);
	}
      while(ya_esta(ix, ii, m));
      memcpy(XC[i], A[ii], sizeof(double) * n);
      ix[i] = ii;
    }
}
  
void new_XC(double **XC, double **A, int *ix, int m, int n, int x)
{
  for(int i = 0; i < m; ++i)
    {
      int ii = ix[i];
      memcpy(XC[i], A[ii], sizeof(double) * n);
    }
}


void retorna_pesos(double **XC, int m, int n, double *Ai, double *pesos, int distribution) {
  int i, j;
  double da, sp = 0, dmax = 0;
  if(m == 1)
    pesos[0] = 1.0;
  else
    {
      for (i = 0; i < m - 1; i++)
	for (j = i + 1; j < m; j++)
	  {
	    da = distancia(XC[i], XC[j], n, distribution);
	    if (da > dmax)
	      dmax = da;
	  }
      for (i = 0; i < m; i++)
	{
	  da = 1 - fmin(distancia(XC[i], Ai, n, distribution) / dmax, 1.0);
	  //1.0 / fmax(fmin(distancia(XC[i], Ai, n, distribution), 1e24), 1e-24);
	  pesos[i] = da;
	  sp += da;
	}
      if(sp > 0)
	for (i = 0; i < m; i++)
	  pesos[i] /= sp;
      else
	for (i = 0; i < m; i++)
	  pesos[i] = 1.0 / m;
    }
}

double inicia_S(double **XC, double **A, double **S, int m, int n, int x, int distribution) {
  int i;
  for (i = 0; i < x; ++i)
      retorna_pesos(XC, m, n, A[i], S[i], distribution);
  return error_t(XC, S, A, m, n, x, distribution);
}

double error_g(double Se[], double **XC, double *Ai, int m, int n, int distribution)
{
  double E[n];
  estima_Si(XC, Se, E, m, n);
  return distancia(E, Ai, n, distribution);
}

double desciende_s(double delta, double *gr, double **XC, double *Se, double *Si, double *Ai,
		 int m, int n, int distribution)
{
  memcpy(Si, Se, m * sizeof(double));

  double sum = 0;
  for(int i = 0; i < m; ++i)
    {
      Si[i] = fmax(Si[i] - delta * gr[i], 0);
      sum += Si[i];
    }
  for(int i = 0; i < m; ++i)
    {
      Si[i] /= sum;
      //printf("Si %12.10f Se %12.10f\n", Si[i], Se[i]);
    }
  return error_g(Si, XC, Ai, m, n, distribution);
}

double golden_linesearchs(double delta, double *gr, double **XC, double *Se, double *Si, double *Ai,
			  int m, int n, int distribution)
{
  double gor = (sqrt(5) + 1) / 2;
  double a = 0, b = delta, tol=1e-8;
  double c, d, fb;
  double e0 = error_g(Se, XC, Ai, m, n, distribution);
  fb = desciende_s(b, gr, XC, Se, Si, Ai, m, n, distribution);
  
  //printf("delta %12.10f e0 %12.10f b %12.10f fb %12.10f\n", delta, e0, b, fb); fflush(stdout);

  while(fb < e0 && b < 100)
    {
      b *= gor;
      fb = desciende_s(b, gr, XC, Se, Si, Ai, m, n, distribution);
      //printf("e0 %12.10f b %12.10f fb %12.10f\n", e0, b, fb); fflush(stdout);
    }

  while (fabs(b - a) > tol)
    {
      c = b - (b - a) / gor;
      d = a + (b - a) / gor;
      double fc = desciende_s(c, gr, XC, Se, Si, Ai, m, n, distribution);
      double fd = desciende_s(d, gr, XC, Se, Si, Ai, m, n, distribution);
      if (fc < fd) b = d;
      else         a = c;
      //printf("a %12.10f b %12.10f c %12.10f d %12.10f eb %12.10f ec %12.10f ed %12.10f e0 %12.10f\n",
      //	     a, b, c, d, fb, fc, fd, e0); fflush(stdout);
    }
  delta = (b + a) / 2;
  double ei = desciende_s(delta, gr, XC, Se, Si, Ai, m, n, distribution);
  //printf("e0 %12.10f delta %12.10f ei %12.10f\n", e0, b, ei); fflush(stdout);
  if (ei < e0)
    {
      memcpy(Se, Si, m * sizeof(double));
      e0 = ei;
    }
  return e0;
}

double gradient_s(double Se[], double gr[], double **XC, double *Ai, int m, int n, int distribution)
{
  double df = 1e-6;
  double norma = 0;
  double Si[m];
  
  for(int i = 0; i < m; ++i)
    {
      double e0 = error_g(Se, XC, Ai, m, n, distribution);
      memcpy(Si, Se, m * sizeof(double));
      Si[i] += df;
      for(int j = 0; j < m; ++j)
	Si[j] /= 1 + df;
      gr[i] = (error_g(Si, XC, Ai, m, n, distribution) - e0) / df;

      memcpy(Si, Se, m * sizeof(double));
      Si[i] -= df;
      for(int j = 0; j < m; ++j)
	Si[j] /= 1 - df;
      gr[i] += (e0 - error_g(Si, XC, Ai, m, n, distribution)) / df;
      gr[i] /= 2.0;
      norma += pow(gr[i], 2);
    }
  return norma;
}

double gd_s(double **XC, double *Ai, double *Se, int m, int n, int distribution)
{
  double e0, e4;
  double Si[m];
  double gr[m];
  double norma = 1e9;
  
  e4 = error_g(Se, XC, Ai, m, n, distribution);
  e0 = e4 + 1;
  //printf("\n\nDescenso de gradiente de S e0 = %8.4f\n", e4);fflush(stdout);
  int itera = 0;
  double delta = 1e-6;

  while((norma > 1e-6) && (itera < 1000) && (fabs(e0 - e4) > 1e-8))
    {
      e0 = e4;
      norma = gradient_s(Se, gr, XC, Ai, m, n, distribution);
      
      double maxgrad = gr[0]; 
      for(int i = 0; i < m; ++i)
	{
	  if(fabs(gr[i]) > maxgrad)
	    maxgrad = fabs(gr[i]);
	}

      delta = fmin(0.5/maxgrad, 100);
      e4 = golden_linesearchs(delta, gr, XC, Se, Si, Ai, m, n, distribution);
      //printf("Norma %12.8f Itera %3d maxgrad %8.4f e0 %8.4f e4 = %8.4f df %8.4f\n ", norma, itera, maxgrad, e0, e4,
      //       log10(fabs(e0 - e4))); fflush(stdout);
      ++itera;
    }
  //printf("Error final = %8.4f\n", e0);fflush(stdout);
  return e0;
}

double update_S(double **XC, double **A, double **S, int m, int n, int x, int distribution)
{
  int i;
  double e0, et = 0, E[n];
  if(m <= 1)
    {
      for (i = 0; i < x; ++i)
	{
	  S[i][0] = 1.0;
	  estima_Si(XC, S[i], E, m, n);
	  et += distancia(E, A[i], n, distribution);
	}
      return et;
    }
  else
    {
      for (i = 0; i < x; ++i)
	{
	  //printf("Ajusta %4d %8.4f\n",i, S[i][0]);
	  e0 = gd_s(XC, A[i], S[i], m, n, distribution);
	  et += e0;
	}
      return et;
    }
}

void estima_XCi(double **A, double Ce[], double E[], int n, int x)
{
  for(uint k = 0; k < n; ++k)
    {
      E[k]  = 0;
      for(uint j = 0; j < x; ++j)
	  E[k] += Ce[j] * A[j][k];
    }
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

double error_c(double *Ce, double **XC, double **A, int m, int n, int x, int ix, int distribution)
{
  double E[n];
  estima_XCi(A, Ce, E, n, x);
  return distancia(E, XC[ix], n, distribution);
}

double error_cg(double Ce, double Ci, double E[], double **XC, double **A, int m, int n, int xj, int ix, int distribution)
{
  actualiza_XCi(A, Ce, Ci, E, n, xj);
  return distancia(E, XC[ix], n, distribution);
}

double desciende_c(double delta, double *gr, double **XC, double *Ce, double *Ci, double **A,
		   int m, int n, int x, int ix, int distribution)
{
  memcpy(Ci, Ce, x * sizeof(double));
  double sum = 0;
  for(int i = 0; i < x; ++i)
    {
      Ci[i] = fmax(fmin(Ci[i] - delta * gr[i], 1), 0);
      sum += Ci[i];
    }
  for(int i = 0; i < x; ++i)
    Ci[i] /= sum;
  
  return error_c(Ci, XC, A, m, n, x, ix, distribution);
}

double golden_linesearchx(double delta, double *gr, double **XC, double *Ce, double *Ci, double **A,
			      int m, int n, int x, int ix, int distribution)
{
  double gor = (sqrt(5) + 1) / 2;
  double a = 0, b = delta, tol=1e-8;
  double c, d, fb;
  double e0 = error_c(Ce, XC, A, m, n, x, ix, distribution);
  fb = desciende_c(b, gr, XC, Ce, Ci, A, m, n, x, ix, distribution);
  while(fb < e0 && b < 100.0)
    {
      b *= gor;
      fb = desciende_c(b, gr, XC, Ce, Ci, A, m, n, x, ix, distribution);
    }

  while (fabs(b - a) > tol)
    {
      c = b - (b - a) / gor;
      d = a + (b - a) / gor;
      double fc = desciende_c(c, gr, XC, Ce, Ci, A, m, n, x, ix, distribution);
      double fd = desciende_c(d, gr, XC, Ce, Ci, A, m, n, x, ix, distribution);
      if (fc < fd) b = d;
      else         a = c;
    }
  delta = (b + a) / 2;
  double ei = desciende_c(delta, gr, XC, Ce, Ci, A, m, n, x, ix, distribution);
  if (ei < e0)
    {
      memcpy(Ce, Ci, x * sizeof(double));
      e0 = ei;
    }
  return e0;
}

double gradient_c(double *Ce, double gr[], double **XC, double **Ai, int m, int n, int x, int ix, int distribution)
{
  double df = 1e-6;
  double norma = 0;
  double Ci[x];
  double E[n];
  estima_XCi(Ai, Ce, E, n, x);
  
  for(int i = 0; i < x; ++i)
    {
      double e0 = distancia(E, XC[ix], n, distribution);
      memcpy(Ci, Ce, x * sizeof(double));
      //Ci[i] += delta;
      Ci[i] += df;
      for(int j = 0; j < x; ++j)
      	Ci[j] /= 1 + df;
      gr[i] = (error_cg(Ce[i], Ci[i], E, XC, Ai, m, n, i, ix, distribution) - e0) / df;
      memcpy(Ci, Ce, x * sizeof(double));
      Ci[i] -= df;
      Ci[i] = fmax(Ci[i], 0);
      for(int j = 0; j < x; ++j)
	Ci[j] /= 1 - df;
      gr[i] += (e0 - error_cg(Ce[i], Ci[i], E, XC, Ai, m, n, i, ix, distribution)) / df;
      gr[i] /= 2.0;
      norma += pow(gr[i], 2);
    }
  return norma;
}


double gd_c(double **XC, double **Ai, double *Ce, int m, int n, int x, int ix, int distribution)
{
  double e0, e4;
  double gr[x];
  double norma = 1e9;
  double Ci[x];
  
  e4 = error_c(Ce, XC, Ai, m, n, x, ix, distribution);
  e0 = e4 + 1;
  int itera = 0;
  double delta = 1e-6;
  
  while(norma > 1e-4 && itera < 1000 && fabs(e0 - e4) > 1e-8)
    {
      e0 = e4;
      norma = gradient_c(Ce, gr, XC, Ai, m, n, x, ix, distribution);

      double maxgrad = gr[0]; 
      for(int i = 0; i < m; ++i)
	{
	  if(fabs(gr[i]) > maxgrad)
	    maxgrad = fabs(gr[i]);
	}
      
      delta  = fmin(0.1/maxgrad, 10);
      e4 = golden_linesearchx(delta, gr, XC, Ce, Ci, Ai, m, n, x, ix, distribution);
      //printf("norma %10.8f e0 %10.8f e4 %10.8f df %10.8f\n", norma, e0, e4, e0 - e4); fflush(stdout);
      ++itera;
    }
  //printf("e0 = %8.4f\n", e0);fflush(stdout);
  return e0;
}


double update_C(double **XC, double **A, double **C, int m, int n, int x, int distribution)
{
  int i;
  double E[n];
  double e0 = 0;
  
  for (i = 0; i < m; ++i)
    {
      estima_XCi(A, C[i], E, n, x);
      e0 += gd_c(XC, A, C[i], m, n, x, i, distribution);
      estima_XCi(A, C[i], XC[i], n, x);
    }
  return e0;
}

double corrige_angulo(double X, double ang_max)
{
  while((X < 0) || (X > ang_max))
    {
      if (X < 0)   X = -X;
      if (X > ang_max) X = X - ang_max;
    }
  return X;
}

void matmul(double **A, double **B, double **C, int x, int y, int z) {
  for (int i = 0; i < x; i++)
    for (int j = 0; j < z; j++)
      {
	C[i][j] = 0;
	for (int k = 0; k < y; k++)
	  C[i][j] += A[i][k] * B[k][j];
      }
}

void transpose(double **A, double **AT, int m, int n) {
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      AT[j][i] = A[i][j];
}

int inverse (double **a, double **b, int n)
{
  int i, j, k, p;
  double f, g, tol;
  if (n < 1) return -1;  /* Function Body */
  f = 0.;  /* Frobenius norm of a */
  for (i = 0; i < n; ++i)
    for (j = 0; j < n; ++j)
      {
	g = a[i][j];
	f += g * g;
      }
  f = sqrt(f);
  tol = f * 2.2204460492503131e-016;
  
  /* Set b to identity matrix. */
  for (i = 0; i < n; ++i)
    for (j = 0; j < n; ++j)
      {
	b[i][j] = (i == j) ? 1. : 0.;
      }
  for (k = 0; k < n; ++k)
    {  /* Main loop */
      f = fabs(a[k][k]);  /* Find pivot. */
      p = k;
      for (i = k+1; i < n; ++i)
	{
	  g = fabs(a[i][k]);
	  if (g > f)
	    {
	      f = g;
	      p = i;
	    }
	}
      if (f < tol) return 1;  /* Matrix is singular. */
      if (p != k)
	{  /* Swap rows. */
	  for (j = k; j < n; ++j)
	    {
	      f = a[k][j];
	      a[k][j] = a[p][j];
	      a[p][j] = f;
	    }
	  for (j = 0; j < n; ++j) {
	    f = b[k][j];
	    b[k][j] = b[p][j];
	    b[p][j] = f;
	  }
	}
      f = 1. / a[k][k];  /* Scale row so pivot is 1. */
      for (j = k; j < n; ++j) a[k][j] *= f;
      for (j = 0; j < n; ++j) b[k][j] *= f;
      for (i = 0; i < n; ++i)
	{  /* Subtract to get zeros. */
	  if (i == k) continue;
	  f = a[i][k];
	  for (j = k; j < n; ++j) a[i][j] -= a[k][j] * f;
	  for (j = 0; j < n; ++j) b[i][j] -= b[k][j] * f;
	}
    }
  return 0;
} /* end of gjinv */

void solve_linear_system(double **S, double **A, double **XC, int m, int n, int x)
{
  double *STm     = (double *)  malloc(m * x * sizeof(double));
  double **ST     = (double **) malloc(m * sizeof(double *));

  double *STAm    = (double *)  malloc(m * m * sizeof(double));
  double **STA    = (double **) malloc(m * sizeof(double *));
  
  double *STIm    = (double *)  malloc(m * m * sizeof(double));
  double **STI    = (double **) malloc(m * sizeof(double *));

  double *STXm    = (double *)  malloc(m * n * sizeof(double));
  double **STX    = (double **) malloc(n * sizeof(double *));

  transfer_m(ST,  STm,  m, x);
  transfer_m(STA, STAm, m, m);
  transfer_m(STI, STIm, m, m);
  transfer_m(STX, STXm, m, n);
  
  transpose(S, ST, x, m);
  matmul (ST,  S, STA, m, x, m);
  inverse(STA, STI, m);  
  matmul (ST,   A, STX, m, x, n);
  matmul (STI, STX, XC, m, m, n);
    
  free(ST);
  free(STm);
  
  free(STA);
  free(STAm);

  free(STI);
  free(STIm);

  free(STX);
  free(STXm);
}

double gradient_XCS(double **GXC, double **XC, double **A, double **S, int m, int n, int x, int distribution)
{
  double delta = 1e-6;
  double eo    = error_t(XC, S, A, m, n, x, distribution);
  double ei;
  double norma = 0;
  
  double *Ea = (double *) malloc(n * x * sizeof(double));
  double **E = (double **) malloc(x * sizeof(double *));
  
  transfer_m(E, Ea, x, n);
  estima_E(E, XC, S, m, n, x);
  //printf("\n");
  
  for(uint i = 0; i < n; ++i)
    for(uint j = 0; j < m; ++j)
      {
	double Xo = XC[j][i];
	XC[j][i] += delta;
	if (distribution == 0)
	  XC[j][i] = corrige_angulo(XC[j][i], 360);
	//ei = error_t(XC, S, A, m, n, x, distribution);
	ei = eo + error_ij(E, Xo, XC[j][i], S, A, m, n, x, distribution, i, j);
	//printf("Error Original %12.8f Error nuevo %12.8f, Diff %12.8f\n", ei, ej, ei - ej);fflush(stdout);
	GXC[j][i] = (ei - eo) / delta;
	XC[j][i] = Xo - delta;
	if (distribution == 0)
	  XC[j][i] = corrige_angulo(XC[j][i], 360);
	//ei = error_t(XC, S, A, m, n, x, distribution);
	ei = eo + error_ij(E, Xo, XC[j][i], S, A, m, n, x, distribution, i, j);
	//printf("Error Original %12.8f Error nuevo %12.8f, Diff %12.8f\n", ei, ej, ei - ej);fflush(stdout);
	GXC[j][i] += (eo - ei) / delta;
	GXC[j][i] /= 2.0;
	XC[j][i] = Xo;
	norma += pow(GXC[j][i], 2);
      }
  free(Ea);
  free(E);
  return norma;
}

double desciende(double delta, double **GX, double **Xd, double **XC, double **S, double **A,
		 int m, int n, int x, int distribution)
{
  memcpy(Xd[0], XC[0], n * m * sizeof(double));
  for(uint i = 0; i < n; ++i)
    for(uint j = 0; j < m; ++j)
      {
	Xd[j][i] -= delta * GX[j][i];
	if (distribution == 0)
	  Xd[j][i] = corrige_angulo(Xd[j][i], 360);
      }
  return error_t(Xd, S, A, m, n, x, distribution);
}

double golden_linesearchxc(double delta, double e0, double **GX, double **Xd, double **XC, double **S, double **A,
			  int m, int n, int x, int distribution)
{
  double gr = (sqrt(5) + 1) / 2;
  double a = 0, b = delta, tol=1e-6;
  double c, d, fb;
  /*
    """Golden-section search
    to find the minimum of f on [a,b]
  */
  fb = desciende(b, GX, Xd, XC, S, A, m, n, x, distribution);
  while(fb < e0)
    {
      b *= gr;
      fb = desciende(b, GX, Xd, XC, S, A, m, n, x, distribution);
    }

  while (fabs(b - a) > tol)
    {
      c = b - (b - a) / gr;
      d = a + (b - a) / gr;
      double fc = desciende(c, GX, Xd, XC, S, A, m, n, x, distribution);
      double fd = desciende(d, GX, Xd, XC, S, A, m, n, x, distribution);
      if (fc < fd) b = d;
      else         a = c;
      //printf("a %8.4f b %8.4f c %8.4f d %8.4f eb %8.4f ec %8.4f ed %8.4f e0 %8.4f\n",
      //log10(a), log10(b), log10(c), log10(d), fb, fc, fd, e0); fflush(stdout);
    }
  delta = (b + a) / 2;
  double ei = desciende(delta, GX, Xd, XC, S, A, m, n, x, distribution);
  if (ei < e0)
    {
      memcpy(XC[0], Xd[0], m * n * sizeof(double));
      e0 = ei;
    }
  return e0;
}
  
void gd_XCS(double **XC, double **A, double **S, int m, int n, int x, int distribution)
{
  double *GXC = (double *) malloc(n * m * sizeof(double));
  double **GX = (double **) malloc(m * sizeof(double *));
  double *XCd = (double *) malloc(n * m * sizeof(double));
  double **Xd = (double **) malloc(m * sizeof(double *));
  
  transfer_m(GX, GXC, m, n);
  transfer_m(Xd, XCd, m, n);
  
  double e0 = error_t(XC, S, A, m, n, x, distribution), ei;
  int itera = 0;
  double delta = 1e-6;
  int acepta = 1;
  double norma = 1e9;
  while(norma > 1e-4 && itera < 100 && acepta == 1)
    {
      acepta = 0;
      norma = gradient_XCS(GX, XC, A, S, m, n, x, distribution);
      double maxgrad = GX[0][0]; 
      for(uint i = 0; i < n; ++i)
	for(uint j = 0; j < m; ++j)
	  {
	    if(fabs(GX[j][i]) > maxgrad)
	      maxgrad = fabs(GX[j][i]);
	  }
      delta  = fmin(1e-2/maxgrad, 0.1);
      ei = golden_linesearchxc(delta, e0, GX, Xd, XC, S, A, m, n, x, distribution);
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


double calcula_XCS(double **XC, double **A, double **S, int m, int n, int x, int distribution)
{
  if(distribution == 0)
      gd_XCS(XC, A, S, m, n, x, distribution);
  else
    solve_linear_system(S, A, XC, m, n, x);
  return error_t(XC, S, A, m, n, x, distribution);
}

double normal_log_likelihood(double SSE, int n) {
  //int dof = n - k;
    double sigma_squared = SSE / n;
    double log_likelihood = -(0.5 * n) * (log((2 * M_PI) * sigma_squared) + 1);
    //return -(n/2) * log((2*) sigma_squared);
    return log_likelihood;
}

void like_rss(double SSE, int n, int k, double* results) {
    double normal_like = normal_log_likelihood(SSE, n);
    double bic  = -2 * normal_like + k * log(n);
    double aicc = -2 * normal_like + 2 * k + 2 * k * (k + 1) / (n - k - 1);
    
    results[0] = normal_like;
    results[1] = bic;
    results[2] = aicc;
}


  
double AllCombinations(double **XC, double **A, double **S, int *ix, int *ixf, int m2, int m, int n, int x, int distribution)
{
  int indices[m], ij[m];
  double SSx = 1e38;
  
  // Initialize the first combination
  for (int i = 0; i < m; i++) indices[i] = i;
  
  while (1) {
    for (int i = 0; i < m; i++)
	ij[i] = ix[indices[i]];
        
    new_XC(XC, A, ij, m, n, x);
    double SSi = inicia_S(XC, A, S, m, n, x, distribution);
    if (SSi < SSx)
      {
	memcpy(ixf, ij, m * sizeof(int));
	SSx = SSi;
      }
	
    // Find the rightmost index that can be incremented
    int i;
    for (i = m - 1; i >= 0; i--)
      {
	if (indices[i] < m2 - m + i) break;
      }

    // If no such index is found, all combinations have been generated
    if (i < 0)
      {
	new_XC(XC, A, ixf, m, n, x);
	return SSx;
      }
    // Increment the current index and reset the indices to the right
    indices[i]++;
    for (int j = i + 1; j < m; j++) indices[j] = indices[j - 1] + 1;
  }
  return SSx;
}

int mini(int x, int y)
{
  return (x < y) ? x : y;
}

void furthest_sum(double **XC, double **A, int *ix, int m, int n, int x, int distribution, int r, int ii)
{
  if(r == 1 || m == 1)
    {
      if (m == 1)
	{
	  memcpy(XC[0], A[ii], sizeof(double) * n);
	  ix[0] = ii;
	}
      else
	init_XC(XC, A, ix, m, n, x);
    }
  else
    {
      //printf("Inicia con FS\n"); fflush(stdout);
      for(int i = 0; i < m; ++i) ix[i] = 0;
      init_XC(XC, A, ix, m, n, x);
      
      for(int i = 0; i < m; ++i)
	{
	  ii = bucle_lejano(XC, A, ix, m, n, x, distribution);
	  update_XC(XC, A, ix, m, n, ii);
	}
    }
}

void init_C(double **C, int *ix, int m, int x)
{
  for(int i = 0; i < m; ++i)
    for(int j = 0; j < x; ++j)
      C[i][j] = (double) 0.01 / x;

  for(int i = 0; i < m; ++i)
    C[i][ix[i]] += (double) 0.99;
}


double Archetypal(double *Ae, double *XCe, double *Ce, double *Se,
		  int m, int n, int x, int distribution, double *IE,
		  int inicializado, int tipo_inicia)
{
  //int tipo_inicia = 1;
  double SSe = 1e99, SSi = 1e99, SSs = 1e99;
  //double IE[3]; //information indexes
  double llx = -1e38, lli = -1e38, lls = -1e38, llz = -1e38, aic = 1e38, bic = 1e38;
  double *A[x], *XC[m], *C[m], *S[x];
  transfer_m(A,   Ae, x, n);
  transfer_m(XC, XCe, m, n);
  transfer_m(C,   Ce, m, x);
  transfer_m(S,   Se, x, m);
  srand48(time(NULL));
  
  if(!(inicializado))
  {
    double *XCr = (double *) malloc(n * m * sizeof(double));
    double *Sr  = (double *) malloc(x * m * sizeof(double));
    double **XCi = (double **) malloc(m * sizeof(double *));
    double **Si  = (double **) malloc(x * sizeof(double *));
    int ix[m];
    transfer_m(XCi, XCr, m, n);
    transfer_m(Si,   Sr, x, m);

    //initiate the archetypes
    int i = (x - 1) * (tipo_inicia == 0) + (1000 * m) * (tipo_inicia == 1);
    //printf("i %d tipo_inicia %d\n",i,tipo_inicia); fflush(stdout);
    while(i >= 0)
      {
	furthest_sum(XCi, A, ix, m, n, x, distribution, 1 - tipo_inicia, i);
	SSi = inicia_S(XCi, A, Si, m, n, x, distribution);
	//if (tipo_inicia == 1) SSi = update_S(XCi, A, Si, m, n, x, distribution);
	like_rss(SSi, x * n, m * x, IE);
	lli = IE[0];
	
	if(SSi < SSe)
	  {
	    memcpy(Se,  Sr,  x * m * sizeof(double));
	    memcpy(XCe, XCr, n * m * sizeof(double));
	    init_C(C, ix, m, x);
	    SSe = SSi;
	    llx = lli;
	    aic = IE[2];
	    bic = IE[1];
	    if(m > 1)
	      i += 1000 * (tipo_inicia == 0) + 100 * (tipo_inicia == 1);
	  }
	printf("\rNOC = %2d P %3d N %4d K %8d Error FS %6d  Llike ite %13.2f Llike max %13.2f AICc %13.2f BIC %13.2f",
	       m,  n, n * x, m * x, i, lli, llx, aic, bic);
	fflush(stdout);
	--i;
      }
    free(XCr); free(XCi);
    free(Sr);  free(Si);
    SSs = update_S(XC, A, S, m, n, x, distribution);
    like_rss(SSs, n * x, m * x, IE);
    lls = IE[0];
    llz = lls;
    printf("\rNOC = %2d P %3d N %4d K %8d Error FS %6d  Llike upt %13.2f Llike max %13.2f AICc %13.2f BIC %13.2f",
	       m,  n, n * x, m * x, -1, lls, llx, aic, bic);
    fflush(stdout);
  }
  int i = 0;
  SSi = calcula_XCS(XC, A, S, m, n, x, distribution);
  like_rss(SSi, n * x, m * x, IE);
  lli = IE[0];
  llx = lli;
  llz = lls;
  printf("\rNOC = %2d P %3d N %4d K %8d Error Iteracion %4d Llike max %13.2f Llike maz %13.2f dflikex %13.4f dflikez %13.4f  AICc %13.2f BIC %13.2f", m,  n, n * x, m * x, -1, lls, llz, lli - llx, lls - llz, IE[2], IE[1]);
  do
    {
      SSe = SSi;
      llx = lli;
      llz = lls;
      update_C(XC, A, C, m, n, x, distribution);
      SSs = update_S(XC, A, S, m, n, x, distribution);
      SSi = calcula_XCS(XC, A, S, m, n, x, distribution);
      like_rss(SSi, n * x, m * x, IE);
      lli = IE[0];
      like_rss(SSs, n * x, m * x, IE);
      lls = IE[0];
      printf("\rNOC = %2d P %3d N %4d K %8d Error Iteracion %4d Llike max %13.2f Llike maz %13.2f dflikex %13.4f dflikez %13.4f  AICc %13.2f BIC %13.2f", m,  n, n * x, m * x, i, lls, llz, lli - llx, lls - llz, IE[2], IE[1]);
      fflush(stdout);
      ++i;
    }
  while((fabs(lli - llx) > 1e-4 || fabs(lls - llz) > 1e-4) && i < 1000);
  printf("\n");fflush(stdout);
  return SSi;
}




double Archetypal_E(double *Ai, double *XCi, double *Ci, double *Si,
		    int m, int n, int x, int distribution, double *IE)
{

  double SSi = 99;
  double *A[x], *XC[m], *C[m], *S[x];
  
  transfer_m(A,   Ai, x, n);
  transfer_m(XC, XCi, m, n);
  transfer_m(C,   Ci, m, x);
  transfer_m(S,   Si, x, m);
      
  update_C(XC, A, C, m, n, x, distribution);
  SSi = update_S(XC, A, S, m, n, x, distribution);
  //SSi = calcula_XCS(XC, A, S, m, n, x, distribution);
  like_rss(SSi, x * n, m * x, IE);
  printf("\rNOC = %2d P %3d N %4d K %8d Error %13.2f AICc %13.2f BIC %13.2f\n",
	 m,  n, n * x, m * x, IE[0], IE[2], IE[1]);

  return SSi;
}
