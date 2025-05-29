#include <math.h>
#include <stdio.h>

// Teepanis Chachiyo, Simple and accurate complete elliptic 
// integrals for the full range of modulus. 2025
//

//
// Compute the exact elliptic integrals and their derivatives
//
void AGM(
double k, // input
double *K, double *E, double *dKdk, double *dEdk){ // output

    double a,b,a0,b0,tn,S;

    // special cases
    if(k==0){ *K = M_PI/2; *E = M_PI/2;
              *dKdk = 0; *dEdk = 0; }
    if(k==1){ *K = INFINITY; *E = 1;
              *dKdk = INFINITY; *dEdk = -INFINITY; }

    if(0 < k && k < 1){
        a0 = 1; b0 = sqrt(1-k*k);
        tn = 1; S  = k*k;
        while(1){
            a = (a0+b0)/2; b = sqrt(a0*b0);
            tn = 2*tn; S = S + tn*(a0-b0)*(a0-b0)/4;
            if(a0==a) break;
            a0 = a; b0 = b;
        }
        *K = M_PI/(2*a); *E = *K *(1-S/2);
        *dKdk = (*E-(1-k*k)*(*K))/(k*(1-k*k));
        *dEdk = ((*E)-(*K))/k;
    }
}

//
// simple and accurate series
//
double K(double k){
    //double n = (log(4) - log(M_PI))/(M_PI/2 - log(4));
    //double b = exp(n*M_PI/2) - pow(4,n);
    // --- pre-computed for speed ---
    double n = 1.3092785997521463;
    double b = 1.678061276031407;

    // avoid numerical divide-by-zero
    if(k==1) return INFINITY;

    return 1/n*log(pow(4/sqrt(1-k*k),n)+b);
}

double invK(double K){
    //double n = (log(4) - log(M_PI))/(M_PI/2 - log(4));
    //double b = exp(n*M_PI/2) - pow(4,n);
    // --- pre-computed for speed ---
    double n = 1.3092785997521463;
    double b = 1.678061276031407;

    // avoid numerical divide-by-zero
    if(K==M_PI/2) return 0;

    return sqrt(1-16/pow(exp(n*K)-b,2/n));
}

double E(double k){
    //double n = log(3*M_PI/2 - 4)/(log(4) - M_PI + 3./2);
    //double b = exp(n*(M_PI-2)) - pow(4/sqrt(M_E),n);
    // --- pre-computed for speed ---
    double n = 1.3283723627880784;
    double b = 1.3103755722411727;

    // avoid numerical divide-by-zero
    if(k==1) return 1;

    return 1 + 1./2*(1-k*k)/n*log( 
        pow(4/sqrt(M_E)/sqrt(1-k*k),n) + b );
}

//
// Compute exact inverse of K
//
void invK_Newton(
double K, // input
double *k, int *nIter){ // output

    double dk, Kn, En, dKdk_n, dEdk_n;
    int MAX=100;

    // special case
    if(K==M_PI/2){ *k = 0.0; return; }

    // initial guess
    *k = invK(K);

    for((*nIter)=1; (*nIter)<MAX; (*nIter)++){
        AGM(*k, &Kn, &En, &dKdk_n, &dEdk_n);
        dk = -(Kn-K)/dKdk_n;
        *k = *k + dk;
        if(fabs(dk)<1e-14) break;
    }
}

int main(){
    double k;
    double K_exact,E_exact,dKdk_exact,dEdk_exact;
    int nIter;

    // test the AGM method
    k = 0.5;
    AGM(k, &K_exact, &E_exact, &dKdk_exact, &dEdk_exact);
    printf("k = %15.10f\n",k);
    printf("Exact K, E            : %15.10f %15.10f\n", K_exact, E_exact);
    printf("Simple & accurate K, E: %15.10f %15.10f\n", K(k), E(k));
    printf("Exact dKdk, dEdk      : %15.10f %15.10f\n", dKdk_exact, dEdk_exact);

    // test Newton root finding
    printf("\n");
    invK_Newton(K_exact, &k, &nIter);
    printf("Exact inverse of K    : %15.10f\n", k);
    printf("Simple & accurate invK: %15.10f\n", invK(K_exact));
}
