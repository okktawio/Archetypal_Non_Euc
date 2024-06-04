from numpy import *
set_printoptions(precision = 4, suppress = True, linewidth = 300, threshold = 100000)


import ctypes
import sys
import time
import fcntl

#if "FS" in sys.argv:
#    # Load the shared library
#    #print("FS")
#lib = ctypes.CDLL('./arch_auxfs.so')
#else:
#    #print("RA")
lib = ctypes.CDLL('./arch_aux.so')


# Define the function argument and return types
lib.Archetypal.argtypes = [ctypes.POINTER(ctypes.c_double),
                           ctypes.POINTER(ctypes.c_double),
                           ctypes.POINTER(ctypes.c_double),
                           ctypes.POINTER(ctypes.c_double),
                           ctypes.c_int,
                           ctypes.c_int,
                           ctypes.c_int,
                           ctypes.c_int,
                           ctypes.POINTER(ctypes.c_double),
                           ctypes.c_int,
                           ctypes.c_int]
lib.Archetypal.restype = ctypes.c_double

lib.Archetypal_E.argtypes = [ctypes.POINTER(ctypes.c_double),
                           ctypes.POINTER(ctypes.c_double),
                           ctypes.POINTER(ctypes.c_double),
                           ctypes.POINTER(ctypes.c_double),
                           ctypes.c_int,
                           ctypes.c_int,
                           ctypes.c_int,
                           ctypes.c_int,
                           ctypes.POINTER(ctypes.c_double)]
lib.Archetypal_E.restype = ctypes.c_double

def normal_log_likelihood(SSE, n):
    sigma_squared  = SSE / n
    log_likelihood = -(0.5 * n) * (log((2 * pi) * sigma_squared) + 1)
    return log_likelihood

def like_rss(SSE, M, S, C):
    n    = M.size
    k    = S.size
    like = normal_log_likelihood(SSE, n)
    bic  = -like + k * log(n)
    aicc = -2 * like + 2 * k + 2 * k * (k + 1) / (n - k - 1)
    return like, bic, aicc

def Arch_submuestra(datos, noc = 3, distribution = 0, aic0 = 1e99, n = 1000):
    vrv = ["fases", "potencias"]
    archivo = "arquetipos_" + vrv[distribution] + "_%03d_%02d.npz"%(6, noc)
    if distribution == 0:
        datos[where(datos < 0)] += 360
    
    A = float64(datos[::100,:].copy())
    XC = zeros((noc, A.shape[1]))
    Cr = zeros((noc, A.shape[0]))
    Sr = zeros((A.shape[0], noc))
    IE = zeros((3,))
    
    XC_ptr = XC.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    A_ptr  =  A.ctypes.data_as (ctypes.POINTER(ctypes.c_double))
    C_ptr  = Cr.ctypes.data_as (ctypes.POINTER(ctypes.c_double))
    S_ptr  = Sr.ctypes.data_as (ctypes.POINTER(ctypes.c_double))
    IE_ptr = IE.ctypes.data_as (ctypes.POINTER(ctypes.c_double))

    SSe = lib.Archetypal(A_ptr, XC_ptr, C_ptr, S_ptr,
                         XC.shape[0], XC.shape[1], A.shape[0], distribution, IE_ptr, 0)

    A = float64(datos.copy())
    ix = arange(A.shape[0]) 
    C =  zeros((noc, A.shape[0]))
    C[:,ix[::100]] += Cr
    S =  zeros((A.shape[0], noc)) + 1/noc
    S[ix[::100],:] = Sr
    IE = zeros((3,))
    
    XC_ptr = XC.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    A_ptr  = A.ctypes.data_as (ctypes.POINTER(ctypes.c_double))
    C_ptr  = C.ctypes.data_as (ctypes.POINTER(ctypes.c_double))
    S_ptr  = S.ctypes.data_as (ctypes.POINTER(ctypes.c_double))
    IE_ptr = IE.ctypes.data_as (ctypes.POINTER(ctypes.c_double))

    SSe = lib.Archetypal(A_ptr, XC_ptr, C_ptr, S_ptr,
                         XC.shape[0], XC.shape[1], A.shape[0], distribution, IE_ptr, 1)

    return IE

def Arch(datos, noc = 3, distribution = 0, aic0 = 1e99, FS = False):
    vrv = ["fases", "potencias"]
    datos = float64(datos)
    archivo = "arquetipos_" + vrv[distribution] + "_%03d_%02d.npz"%(6, noc)
    
    if distribution == 0:
        datos[where(datos < 0)] += 360

    XC = zeros((noc, datos.shape[1]))
    C = zeros((noc, datos.shape[0]))
    S = zeros((datos.shape[0], noc))
    IE = zeros((3,))
    
    XC_ptr = XC.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    A_ptr  = datos.ctypes.data_as (ctypes.POINTER(ctypes.c_double))
    C_ptr  = C.ctypes.data_as (ctypes.POINTER(ctypes.c_double))
    S_ptr  = S.ctypes.data_as (ctypes.POINTER(ctypes.c_double))
    IE_ptr = IE.ctypes.data_as (ctypes.POINTER(ctypes.c_double))

    #print("FS", FS)
    #if not FS:
    SSe = lib.Archetypal(A_ptr, XC_ptr, C_ptr, S_ptr,
                         XC.shape[0], XC.shape[1], datos.shape[0], distribution, IE_ptr, 0, 0)

    #else:
    #    SSe = lib.Archetypal(A_ptr, XC_ptr, C_ptr, S_ptr,
    #                         XC.shape[0], XC.shape[1], datos.shape[0], distribution, IE_ptr, 0, 1)

    
    lk, bic, aic = IE

    aic1 = 1e99
    #print("Terminado", IE, aic0)
    if aic < aic0:
        try:
            fl = open(archivo, "rb")
            fcntl.flock(fl, fcntl.LOCK_EX)
            aic1 = load(fl)["ajuste"][2]
            fcntl.flock(fl, fcntl.LOCK_UN)
            fl.close()
        except IOError: aic1 = 1e99
        #print("Obtenido AIC1", IE, aic0, aic1)
        
        open("arquetipos_" + vrv[distribution] + "_%02d.log"%(noc),"a").write("aic %12.4f aic1 %12.4f\n"%(aic, aic1))
                                                                                     
        if aic < aic1:
            #print("Guardando", IE, aic0, aic1)
            fl = open(archivo, "wb")
            fcntl.flock(fl, fcntl.LOCK_EX)
            savez_compressed(fl, XC = XC, S = S, C = C, ajuste = array([SSe, lk, aic, bic]))
            fcntl.flock(fl, fcntl.LOCK_UN)
            fl.close()
    aic = min(aic0, min(aic, aic1))
    #print("Retornando nuevo AIC", aic)
    return lk, aic, bic

def main(n = -1, distribution = 0, nocmin = 0, nocmax = 5, FS = False):
    if   distribution == 0:
        vfases = load("lagunablanca_vtransformadas_fases.npz")
    elif distribution == 1:
        vfases = load("lagunablanca_vtransformadas_potencias.npz")
    vtr = vfases["vtr"]
    aicm = zeros((10))
    bicm = zeros((10))
    for i in range(nocmin, nocmax):
        aic0 = 1e99
        for j in range(1):
            lk, aic, bic = Arch(vtr.copy(), noc = i + 1, distribution = distribution, FS = FS)
            aic0 = min(aic, aic0)
            aicm[i] = aic0
            bicm[i] = bic

def main_debug(n = -1, distribution = 0, nocmin = 0):
    if   distribution == 0:
        vfases = load("lagunablanca_vtransformadas_fases.npz")
    elif distribution == 1:
        vfases = load("lagunablanca_vtransformadas_potencias.npz")
    vtr = vfases["vtr"]
    print(vtr.shape)
    lk, aic, bic = Arch(vtr[::1000,:2].copy(), noc = nocmin, distribution = distribution)

def bucle_movil(n = 24, distribution = 0):
    if distribution == 0:
        vfases = load("lagunablanca_vtransformadas_fases.npz")
    elif distribution == 1:
        vfases = load("lagunablanca_vtransformadas_potencias.npz")
    vtr = vfases["vtr"]
    aicm = zeros((10, vtr.shape[1] - n))
    bicm = zeros((10, vtr.shape[1] - n))
    for j in range(aicm.shape[1]):
        for i in range(10):
            lk, aic, bic = Arch(vtr[:,j:j+n].copy(), noc = i + 1, nnn = n, distribution = distribution)
            aicm[i,j] = aic
            bicm[i,j] = bic
            print(aicm[:,j].min(), end = " ")
            print(where(aicm[:,j] == aicm[:,j].min()))

def AA(distribution = 0, noc = 5, FS = False):
    if   distribution == 0:
        vfases = load("lagunablanca_vtransformadas_fases.npz")
    elif distribution == 1:
        vfases = load("lagunablanca_vtransformadas_potencias.npz")
    vtr = vfases["vtr"]
    Arch(vtr.copy(), noc = noc, distribution = distribution, FS = FS)
            
if __name__ == '__main__':
    distribution = int(sys.argv[-2])
    noc = int(sys.argv[-1])
    FS = "FS" in sys.argv
    #print(distribution, noc, FS)
    AA(distribution = distribution, noc = noc, FS = FS)



