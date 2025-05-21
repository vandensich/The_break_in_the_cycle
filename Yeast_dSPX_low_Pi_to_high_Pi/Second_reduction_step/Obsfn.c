#include <R.h>
 #include <math.h>
 void Obsfn_0o9zhn6z ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = x[0+i**k] ;
y[1+i**l] = p[16]*x[25+i**k] ;
y[2+i**l] = p[17]*x[26+i**k] ;
y[3+i**l] = p[18]*x[27+i**k] ;
y[4+i**l] = p[19]*x[28+i**k] ;
y[5+i**l] = x[1+i**k] ;
y[6+i**l] = x[2+i**k] ;
y[7+i**l] = x[3+i**k] ;
y[8+i**l] = x[5+i**k] ;
y[9+i**l] = x[10+i**k] ;
y[10+i**l] = x[6+i**k] ;
y[11+i**l] = x[7+i**k]+x[11+i**k] ;
y[12+i**l] = x[8+i**k]+x[16+i**k]+x[12+i**k] ; 
}
}