/** Code auto-generated by cOde 1.1.1 **/
#include <R.h> 
 #include <math.h> 

static double parms[53];
static double forc[0];
static double cons[0];
static double range[2];

#define nGridpoints 2 
#define nSplines 0 
#define precision 1e-05 

#define k2 parms[0] 
 #define Km2 parms[1] 
 #define k3 parms[2] 
 #define Km3 parms[3] 
 #define k4 parms[4] 
 #define Km4 parms[5] 
 #define k5 parms[6] 
 #define Km5 parms[7] 
 #define deg_InsP6 parms[8] 
 #define prod_InsP6 parms[9] 
 #define k7 parms[10] 
 #define Km7 parms[11] 
 #define k9 parms[12] 
 #define Km9 parms[13] 
 #define k6 parms[14] 
 #define Km6 parms[15] 
 #define k8 parms[16] 
 #define Km8 parms[17] 
 #define phospho0 parms[18] 
 #define phospho1 parms[19] 
 #define phospho2 parms[20] 
 #define phospho_0 parms[21] 
 #define phospho_1 parms[22] 
 #define phospho_2 parms[23] 
 #define y0_0 parms[24] 
 #define y1_0 parms[25] 
 #define y2_0 parms[26] 
 #define y3_0 parms[27] 
 #define y4_0 parms[28] 
 #define y5_0 parms[29] 
 #define y6_0 parms[30] 
 #define y7_0 parms[31] 
 #define y8_0 parms[32] 
 #define y9_0 parms[33] 
 #define y10_0 parms[34] 
 #define y11_0 parms[35] 
 #define y12_0 parms[36] 
 #define y13_0 parms[37] 
 #define y14_0 parms[38] 
 #define y15_0 parms[39] 
 #define y16_0 parms[40] 
 #define y17_0 parms[41] 
 #define y18_0 parms[42] 
 #define y19_0 parms[43] 
 #define y20_0 parms[44] 
 #define y21_0 parms[45] 
 #define y22_0 parms[46] 
 #define y23_0 parms[47] 
 #define y24_0 parms[48] 
 #define y25_0 parms[49] 
 #define y26_0 parms[50] 
 #define y27_0 parms[51] 
 #define y28_0 parms[52] 
#define tmin range[0]
#define tmax range[1]


void odemodel_HCT116_low_to_high_Pi_initmod(void (* odeparms)(int *, double *)) {
	 int N=53;
	 odeparms(&N, parms);
}

void odemodel_HCT116_low_to_high_Pi_initforc(void (* odeforcs)(int *, double *)) {
	 int N=0;
	 odeforcs(&N, forc);
}

/** Derivatives (ODE system) **/
void odemodel_HCT116_low_to_high_Pi_derivs (int *n, double *t, double *y, double *ydot, double *RPAR, int *IPAR) {

	 double time = *t;

	 ydot[0] = -1.0*((k2/(Km2+y[0]))*y[25]*y[0])-1.0*((k2/(Km2+y[0]))*y[26]*y[0])-1.0*((k2/(Km2+y[0]))*y[27]*y[0])-1.0*((k2/(Km2+y[0]))*y[28]*y[0])-1.0*((k3/(Km3+y[0]))*y[25]*y[0])-1.0*((k3/(Km3+y[0]))*y[26]*y[0])-1.0*((k3/(Km3+y[0]))*y[27]*y[0])-1.0*((k3/(Km3+y[0]))*y[28]*y[0])+1.0*((k4/(Km4+y[1]))*y[1])+1.0*((k4/(Km4+y[2]))*y[2])+1.0*((k4/(Km4+y[3]))*y[3])+1.0*((k4/(Km4+y[4]))*y[4])+1.0*((k5/(Km5+y[5]))*y[5])+1.0*((k5/(Km5+y[6]))*y[6])+1.0*((k5/(Km5+y[7]))*y[7])+1.0*((k5/(Km5+y[8]))*y[8])-1.0*(deg_InsP6*y[0])+1.0*(prod_InsP6);
 	 ydot[1] = 1.0*((k2/(Km2+y[0]))*y[25]*y[0])-1.0*((k4/(Km4+y[1]))*y[1])-1.0*((k7/(Km7+y[1]))*y[25]*y[1])-1.0*((k7/(Km7+y[1]))*y[26]*y[1])-1.0*((k7/(Km7+y[1]))*y[27]*y[1])-1.0*((k7/(Km7+y[1]))*y[28]*y[1])+1.0*((k9/(Km9+y[9]))*y[9])+1.0*((k9/(Km9+y[13]))*y[13])+1.0*((k9/(Km9+y[17]))*y[17])+1.0*((k9/(Km9+y[21]))*y[21]);
 	 ydot[2] = 1.0*((k2/(Km2+y[0]))*y[26]*y[0])-1.0*((k4/(Km4+y[2]))*y[2])-1.0*((k7/(Km7+y[2]))*y[25]*y[2])-1.0*((k7/(Km7+y[2]))*y[26]*y[2])-1.0*((k7/(Km7+y[2]))*y[27]*y[2])-1.0*((k7/(Km7+y[2]))*y[28]*y[2])+1.0*((k9/(Km9+y[10]))*y[10])+1.0*((k9/(Km9+y[14]))*y[14])+1.0*((k9/(Km9+y[18]))*y[18])+1.0*((k9/(Km9+y[22]))*y[22]);
 	 ydot[3] = 1.0*((k2/(Km2+y[0]))*y[27]*y[0])-1.0*((k4/(Km4+y[3]))*y[3])-1.0*((k7/(Km7+y[3]))*y[25]*y[3])-1.0*((k7/(Km7+y[3]))*y[26]*y[3])-1.0*((k7/(Km7+y[3]))*y[27]*y[3])-1.0*((k7/(Km7+y[3]))*y[28]*y[3])+1.0*((k9/(Km9+y[11]))*y[11])+1.0*((k9/(Km9+y[15]))*y[15])+1.0*((k9/(Km9+y[19]))*y[19])+1.0*((k9/(Km9+y[23]))*y[23]);
 	 ydot[4] = 1.0*((k2/(Km2+y[0]))*y[28]*y[0])-1.0*((k4/(Km4+y[4]))*y[4])-1.0*((k7/(Km7+y[4]))*y[25]*y[4])-1.0*((k7/(Km7+y[4]))*y[26]*y[4])-1.0*((k7/(Km7+y[4]))*y[27]*y[4])-1.0*((k7/(Km7+y[4]))*y[28]*y[4])+1.0*((k9/(Km9+y[12]))*y[12])+1.0*((k9/(Km9+y[16]))*y[16])+1.0*((k9/(Km9+y[20]))*y[20])+1.0*((k9/(Km9+y[24]))*y[24]);
 	 ydot[5] = 1.0*((k3/(Km3+y[0]))*y[25]*y[0])-1.0*((k5/(Km5+y[5]))*y[5])-1.0*((k6/(Km6+y[5]))*y[25]*y[5])-1.0*((k6/(Km6+y[5]))*y[26]*y[5])-1.0*((k6/(Km6+y[5]))*y[27]*y[5])-1.0*((k6/(Km6+y[5]))*y[28]*y[5])+1.0*((k8/(Km8+y[9]))*y[9])+1.0*((k8/(Km8+y[10]))*y[10])+1.0*((k8/(Km8+y[11]))*y[11])+1.0*((k8/(Km8+y[12]))*y[12]);
 	 ydot[6] = 1.0*((k3/(Km3+y[0]))*y[26]*y[0])-1.0*((k5/(Km5+y[6]))*y[6])-1.0*((k6/(Km6+y[6]))*y[25]*y[6])-1.0*((k6/(Km6+y[6]))*y[26]*y[6])-1.0*((k6/(Km6+y[6]))*y[27]*y[6])-1.0*((k6/(Km6+y[6]))*y[28]*y[6])+1.0*((k8/(Km8+y[13]))*y[13])+1.0*((k8/(Km8+y[14]))*y[14])+1.0*((k8/(Km8+y[15]))*y[15])+1.0*((k8/(Km8+y[16]))*y[16]);
 	 ydot[7] = 1.0*((k3/(Km3+y[0]))*y[27]*y[0])-1.0*((k5/(Km5+y[7]))*y[7])-1.0*((k6/(Km6+y[7]))*y[25]*y[7])-1.0*((k6/(Km6+y[7]))*y[26]*y[7])-1.0*((k6/(Km6+y[7]))*y[27]*y[7])-1.0*((k6/(Km6+y[7]))*y[28]*y[7])+1.0*((k8/(Km8+y[17]))*y[17])+1.0*((k8/(Km8+y[18]))*y[18])+1.0*((k8/(Km8+y[19]))*y[19])+1.0*((k8/(Km8+y[20]))*y[20]);
 	 ydot[8] = 1.0*((k3/(Km3+y[0]))*y[28]*y[0])-1.0*((k5/(Km5+y[8]))*y[8])-1.0*((k6/(Km6+y[8]))*y[25]*y[8])-1.0*((k6/(Km6+y[8]))*y[26]*y[8])-1.0*((k6/(Km6+y[8]))*y[27]*y[8])-1.0*((k6/(Km6+y[8]))*y[28]*y[8])+1.0*((k8/(Km8+y[21]))*y[21])+1.0*((k8/(Km8+y[22]))*y[22])+1.0*((k8/(Km8+y[23]))*y[23])+1.0*((k8/(Km8+y[24]))*y[24]);
 	 ydot[9] = 1.0*((k6/(Km6+y[5]))*y[25]*y[5])+1.0*((k7/(Km7+y[1]))*y[25]*y[1])-1.0*((k8/(Km8+y[9]))*y[9])-1.0*((k9/(Km9+y[9]))*y[9]);
 	 ydot[10] = 1.0*((k6/(Km6+y[5]))*y[26]*y[5])+1.0*((k7/(Km7+y[2]))*y[25]*y[2])-1.0*((k8/(Km8+y[10]))*y[10])-1.0*((k9/(Km9+y[10]))*y[10]);
 	 ydot[11] = 1.0*((k6/(Km6+y[5]))*y[27]*y[5])+1.0*((k7/(Km7+y[3]))*y[25]*y[3])-1.0*((k8/(Km8+y[11]))*y[11])-1.0*((k9/(Km9+y[11]))*y[11]);
 	 ydot[12] = 1.0*((k6/(Km6+y[5]))*y[28]*y[5])+1.0*((k7/(Km7+y[4]))*y[25]*y[4])-1.0*((k8/(Km8+y[12]))*y[12])-1.0*((k9/(Km9+y[12]))*y[12]);
 	 ydot[13] = 1.0*((k6/(Km6+y[6]))*y[25]*y[6])+1.0*((k7/(Km7+y[1]))*y[26]*y[1])-1.0*((k8/(Km8+y[13]))*y[13])-1.0*((k9/(Km9+y[13]))*y[13]);
 	 ydot[14] = 1.0*((k6/(Km6+y[6]))*y[26]*y[6])+1.0*((k7/(Km7+y[2]))*y[26]*y[2])-1.0*((k8/(Km8+y[14]))*y[14])-1.0*((k9/(Km9+y[14]))*y[14]);
 	 ydot[15] = 1.0*((k6/(Km6+y[6]))*y[27]*y[6])+1.0*((k7/(Km7+y[3]))*y[26]*y[3])-1.0*((k8/(Km8+y[15]))*y[15])-1.0*((k9/(Km9+y[15]))*y[15]);
 	 ydot[16] = 1.0*((k6/(Km6+y[6]))*y[28]*y[6])+1.0*((k7/(Km7+y[4]))*y[26]*y[4])-1.0*((k8/(Km8+y[16]))*y[16])-1.0*((k9/(Km9+y[16]))*y[16]);
 	 ydot[17] = 1.0*((k6/(Km6+y[7]))*y[25]*y[7])+1.0*((k7/(Km7+y[1]))*y[27]*y[1])-1.0*((k8/(Km8+y[17]))*y[17])-1.0*((k9/(Km9+y[17]))*y[17]);
 	 ydot[18] = 1.0*((k6/(Km6+y[7]))*y[26]*y[7])+1.0*((k7/(Km7+y[2]))*y[27]*y[2])-1.0*((k8/(Km8+y[18]))*y[18])-1.0*((k9/(Km9+y[18]))*y[18]);
 	 ydot[19] = 1.0*((k6/(Km6+y[7]))*y[27]*y[7])+1.0*((k7/(Km7+y[3]))*y[27]*y[3])-1.0*((k8/(Km8+y[19]))*y[19])-1.0*((k9/(Km9+y[19]))*y[19]);
 	 ydot[20] = 1.0*((k6/(Km6+y[7]))*y[28]*y[7])+1.0*((k7/(Km7+y[4]))*y[27]*y[4])-1.0*((k8/(Km8+y[20]))*y[20])-1.0*((k9/(Km9+y[20]))*y[20]);
 	 ydot[21] = 1.0*((k6/(Km6+y[8]))*y[25]*y[8])+1.0*((k7/(Km7+y[1]))*y[28]*y[1])-1.0*((k8/(Km8+y[21]))*y[21])-1.0*((k9/(Km9+y[21]))*y[21]);
 	 ydot[22] = 1.0*((k6/(Km6+y[8]))*y[26]*y[8])+1.0*((k7/(Km7+y[2]))*y[28]*y[2])-1.0*((k8/(Km8+y[22]))*y[22])-1.0*((k9/(Km9+y[22]))*y[22]);
 	 ydot[23] = 1.0*((k6/(Km6+y[8]))*y[27]*y[8])+1.0*((k7/(Km7+y[3]))*y[28]*y[3])-1.0*((k8/(Km8+y[23]))*y[23])-1.0*((k9/(Km9+y[23]))*y[23]);
 	 ydot[24] = 1.0*((k6/(Km6+y[8]))*y[28]*y[8])+1.0*((k7/(Km7+y[4]))*y[28]*y[4])-1.0*((k8/(Km8+y[24]))*y[24])-1.0*((k9/(Km9+y[24]))*y[24]);
 	 ydot[25] = -1.0*(phospho0*y[25])-1.0*(phospho1*y[25])-1.0*(phospho2*y[25])+1.0*(phospho_0*y[26])+1.0*(phospho_1*y[27])+1.0*(phospho_2*y[28]);
 	 ydot[26] = 1.0*(phospho0*y[25])-1.0*(phospho_0*y[26]);
 	 ydot[27] = 1.0*(phospho1*y[25])-1.0*(phospho_1*y[27]);
 	 ydot[28] = 1.0*(phospho2*y[25])-1.0*(phospho_2*y[28]);

}

