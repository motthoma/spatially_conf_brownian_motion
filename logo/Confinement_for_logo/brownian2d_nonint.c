/*Modulation of brownian motion in 2-d in overdamped limit with confinement with external force*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>



 
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

double positionxfunc(double a, double b, double c,  double f, double g, double h, double i){
 double x;
  
  x=a + sqrt(2*(g/h)*i*b)*c + f*b;
   
  return(x);
					}

double positionyfunc(double a, double b, double c, double d, double e, double f, double g){
 double y;
  
  y=a + sqrt(2*(g/e)*f*b)*c + d*c;
   
  return(y);
					}


 double xmf(double a, double b){
 double c; 

 c=(a-(floorfunc(a,b)*b));            
 
  return(c);
			      }


double yuef(double a, double b, double c, double d, double e, double f){
double g;

if (a <= c){
    g=b-sqrt(c*c-a*a);
	      } 	    
    if ((c <= a) && (a < e)) {
     g=b+d*(f-a)-c*sqrt(1+d*d);
				} 
      if ((e <= a) && (a < f)) {
	g=b-sqrt(c*c-(a-f)*(a-f));
				    }   
 return (g);
 
			}

void specifications(int b,  double c, double d, double f, double g){
unsigned int a;
double e;

    FILE *outpspecs;
    outpspecs = fopen("extbrownian2dconfspecs.dat", "w");		/*erzeugt File,das dt, die Gesamtdauer und die externen Kräfte enthält*/
    fprintf(outpspecs, "dt=%lf\nTemp.=%lf\nR=%lf\n", c,d,g);
       
    for(a=0; a < b; a++){
      e=f*a-2;
      fprintf(outpspecs, "%d.F=%lf\n", a+1, -e);
     }
    fclose(outpspecs);
}
 
int main (){

    
  double x,y,xm,yue,yuew,yuep,dt,t,yo,xo;
  double L,m,B,R,Lp,sigma,T,f,f_y;
  unsigned int i,j,n,N,P,Pn;
  float c;
  double **positionx, **positiony;

  FILE *outp;
  outp=fopen("extbrownian2dconf.dat" ,"w");  


  gsl_rng *r=gsl_rng_alloc (gsl_rng_mt19937);


  printf("Anzahl der Stuetzstellen:");
  scanf("%d", &n); 
  dt=1.7*pow(10,-3);
  printf("dt=%.9lf\n", dt);

  printf("Anzahl der zu plottenden Trajektorien:");
  scanf("%d", &N);

  printf("Faktor um den die ext. Kraft mit jeder Trajektorie erhoeht wird:");
  scanf("%f",&c);

  /*vector<double> p(1000);  dynamischer Speicher geordnet*/
  /*map<int,double> p; dynamischer Speicher, ungeordnet*/
  /*p[1]=100;
  p.resize(N);*/

  sigma=1;
  T=0.1;
  m=0.9;
  L=1;
  B=0.1;
  R=0.5*B;
  Lp=L-(R*m)/(sqrt(1+m*m));
  
  specifications(N, dt, T, c, R);

  positionx = malloc(n * sizeof(double *));
  positiony = malloc(n * sizeof(double *));

  for(i=0; i<n; i++){
    positionx[i] = malloc(N * sizeof(double)); 
    positiony[i] = malloc(N * sizeof(double));         
  }



  for(j=0;j<N;j++){  
  
    t=0;

    if (j == 0){
	    xo = 1.45;
	    yo = 0.05;
	    R = 0.5*B;
	    f = -1.5;
            f_y = 1.2;
    }
    else {
	    xo = 1.5;
	    yo = 0;
	    R = 0.5*B;
	    f = 2.0; 
            f_y = 0;
    }	
    
    //printf("yue start=%lf\t",yue);
    
    
    gsl_rng_set(r, j+100);  
    for (i=0; i<n; i++){       

      
      positionx[i][j]=xo; 
      positiony[i][j]=yo;
      //printf("yin=%lf xin=%lf ",y,x); 
      P=floorfunc(xo, L);
      
      do{
	double u = gsl_ran_gaussian_ziggurat(r,sigma);
	double v = gsl_ran_gaussian_ziggurat(r,sigma);     
	
	x=positionxfunc(xo, dt, u, f, B, R, T);
	y=positionxfunc(yo, dt, v, f_y, B, R, T);	

	xm=xmf(x,L);  
	yue=yuef(xm,B,R,m,Lp,L); 
      }while(fabs(y) > yue);

      //printf("xm=%lf yue=%lf xout=%lf yout=%lf \n",xm,yue,x,y);      

      Pn=floorfunc(x, L);
      //printf("P=%d, Pn=%d, xo=%lf, x=%lf, yo=%lf, y=%lf ",P,Pn,xo,x,yo,y);

      if ((P != Pn) && ((fabs(yo) >= B-R) || (fabs(y) >= B-R))){
	
	x=xo;
	y=yo;  
	
	do{
	  double u = gsl_ran_gaussian_ziggurat(r,sigma);
	  double v = gsl_ran_gaussian_ziggurat(r,sigma); 
	  x=positionxfunc(xo, dt, u, f, B, R, T);
	  y=positionxfunc(yo, dt, v, f_y, B, R, T);
	  Pn=floorfunc(x, L);
	  xm=xmf(x,L);
	  yuep=yuef(xm,B,R,m,Lp,L);
	}while (((Pn != P) && ((fabs(y) >= B-R) || (fabs(yo) >= B-R))) || ((Pn == P) && (fabs(y) > yuep)));
	
      }
  //printf("xP=%lf yP=%lf\n",x,y);   
    
      yo=y;
      xo=x;
    }
      /* printf("\n");*/              
  }
      

  for(i = 0; i < n; i++){	

    t+=dt;
  
    for(j = 0; j < N; j++){ 
      fprintf(outp, "%lf\t%lf\t", positionx[i][j], positiony[i][j]);
    }
			  
    fprintf(outp, "\n");
  }   

    
  for(i = 0; i < n; i++){
    free (positionx[i]);
  }
  free(positionx);      

  for(i = 0; i < n; i++){
    free (positiony[i]);
  } 
  free (positiony);

    gsl_rng_free(r);  
    r=NULL;
    fclose(outp);    
    return 0;
  }
