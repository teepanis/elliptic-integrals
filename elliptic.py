#!/usr/bin/python3
import numpy as np

# Teepanis Chachiyo, Simple and accurate complete elliptic 
# integrals for the full range of modulus. 2025
#

# AGM K and E to machine precision from Abramowitz_Stegun Eq.17.6.3-4
# If you want to check with scipy be aware that you need to use 
# scipy.special.ellipk(k**2) and scipy.special.ellipe(k**2)
def AGM(k):

    # special cases
    if k==0:
        K = np.pi/2; E = np.pi/2; dKdk = 0; dEdk = 0
    if k==1:
        K = np.inf; E = 1; dKdk = np.inf; dEdk = -np.inf

    if 0 < k and k < 1:
        a0 = 1; b0 = np.sqrt(1-k**2)
        tn = 1; S  = k**2
        while True:
            a  = (a0+b0)/2; b  = np.sqrt(a0*b0)
            tn = 2*tn; S = S + tn*(a0-b0)**2/4

            if a == a0: break
            a0 = a; b0 = b

        K = np.pi/(2*a); E = K*(1-S/2)

        dKdk = (E-(1-k**2)*K)/(k*(1-k**2))
        dEdk = (E-K)/k

    return K, E, dKdk, dEdk

#
# simple and accurate series
#

def K(k):
    #n = (np.log(4) - np.log(np.pi))/(np.pi/2 - np.log(4))
    #b = np.exp(n*np.pi/2) - 4**n
    # -- pre-computed for speed ---
    n = 1.3092785997521463
    b = 1.678061276031407

    # avoid numerical divide-by-zero
    if k==1: return np.inf
    
    return 1/n*np.log( (4/np.sqrt(1-k**2))**n + b)

def E(k):
    #n = np.log(3*np.pi/2 - 4)/(np.log(4) - np.pi + 3/2)
    #b = np.exp(n*(np.pi-2)) - (4/np.sqrt(np.e))**n

    # --- pre-computed for speed ---
    n = 1.3283723627880784
    b = 1.3103755722411727

    # avoid numerical divide-by-zero
    if k==1: return 1

    return 1 + 1/2*(1-k**2)/n*np.log( (4/np.sqrt(np.e)/np.sqrt(1-k**2))**n + b  )

def invK(K):
    #n = (np.log(4) - np.log(np.pi))/(np.pi/2 - np.log(4))
    #b = np.exp(n*np.pi/2) - 4**n
    # -- pre-computed for speed ---
    n = 1.3092785997521463
    b = 1.678061276031407

    # avoid numerical divide-by-zero
    if K==np.pi/2: return 0

    return np.sqrt(1-16/(np.exp(n*K)-b)**(2/n))

def invK_Newton(K):

    # special case
    if K==np.pi/2: return 0,0
    
    # initial guess
    k = invK(K)
    
    # safeguard
    MAX=100

    for nIter in range(1,MAX):
        Kn,En,dKdk_n,dEdk_n = AGM(k)
        dk = -(Kn-K)/dKdk_n
        k = k + dk
        if np.abs(dk)<1e-14: break

    return k, nIter

############
### Main ###
############

# test the AGM method
k = 0.5
K_exact, E_exact, dKdk_exact, dEdk_exact = AGM(k)
print("k = %15.10f" % (k))
print("Exact K, E            : %15.10f %15.10f" % (K_exact, E_exact))
print("Simple & accurate K, E: %15.10f %15.10f" % (K(k), E(k)) )
print("Exact dKdk, dEdk      : %15.10f %15.10f" % (dKdk_exact, dEdk_exact))

# test Newton root finding
print("");
k, nIter = invK_Newton(K_exact)
print("Exact inverse of K    : %15.10f" % (k))
print("Simple & accurate invK: %15.10f" % (invK(K_exact)))