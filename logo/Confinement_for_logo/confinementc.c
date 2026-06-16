#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int floorfunc(double a, double b){
 double zahl;
 int wert,c;
 zahl=a/b;
 if(zahl>=0){
   for(wert=0;wert<=zahl;wert++);
   c=(wert-1);
	    }
	    
 if(zahl<0){
   for(wert=0;wert>=zahl;wert--);
   c=wert;
	   }
 return(c);
				  }




main(){

double x,xm,yue,L,m,b,R,Lp,N,i,f;

FILE *outp;
outp=fopen("confinementc_R_0.05.dat", "w");

N=4000;
m=0.9;
L=1;
b=0.1;
R=0.5*b;
Lp=L-(R*m)/(sqrt(1+m*m));

// printf("Lp=%lf\n",Lp);

for(i=-N/2;i<=N;i++){

x=i*0.001;
xm=x-floorfunc(x,L)*L;

if (xm < R){
    yue=b-sqrt(R*R-xm*xm);
	      } 	    
    if ((R < xm) && (xm < Lp)) {
     yue=b+m*(L-xm)-R*sqrt(1+m*m);
				} 
      if ((Lp <= xm) && (xm <= L)) {
	yue=b-sqrt(R*R-(xm-L)*(xm-L));
				    } 

f=-m*xm+m*L+b;

	fprintf(outp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",x,xm,yue,-yue,f,-f);	 
		 
		 }








fclose(outp);
return 0;
}
