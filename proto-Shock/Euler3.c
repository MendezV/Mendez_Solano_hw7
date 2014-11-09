//
//  Difusion.c
//  
//
//  Created by Juan Felipe Méndez on 23/10/14.
// resuelve las ecuaciones de euler
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double nasty_half_vect(double uione,double ui,double fione,double fi,double dx,double dt);
double nasty_one_vect(double ui,double fihalf,double fiminhalf,double dx,double dt);
double fvect1(double u2);
double fvect2(double u1,double u2,double u3,double gamma);
double fvect3(double u1,double u2,double u3,double gamma);





int main(int argc, char **argv)
{
    int nx,nt,i,n;
    double dx, dt,xi,xf,ti,tf,gamma,et;
    double rho, p, u,x;
    double *u1;
    double *uhalf1;
    double *u2;
    double *uhalf2;
    double *u3;
    double *uhalf3;
    double *f1;
    double *f2;
    double *f3;
    
    
    FILE *in0;
    char filename[100];
    sprintf(filename, "estado_%s.dat", argv[1]);
    in0 = fopen(filename,"w");
    
    if(!in0){
        printf("problems opening the file %s\n", filename);
        exit(1);
    }
  
    nx=1000;
    nt=1000;
    xi=-10.0;
    xf=10.0;
    ti=0.0;
    tf=atof(argv[1]);
    dx=(xf-xi)/nx;
    dt=(tf-ti)/nt;
    gamma =1.4;
    
    
    u1= malloc(nx*sizeof(double));
    uhalf1= malloc(nx*sizeof(double));
    u2= malloc(nx*sizeof(double));
    uhalf2= malloc(nx*sizeof(double));
    u3= malloc(nx*sizeof(double));
    uhalf3= malloc(nx*sizeof(double));
    f1= malloc(nx*sizeof(double));
    f2= malloc(nx*sizeof(double));
    f3= malloc(nx*sizeof(double));


    /*condiciones iniciales*/
    rho=1;
    p=100000;
    u=0;
    et= p/((gamma-1.0)*rho) + pow(u,2)/2.0;
    for(i=0;i<(nx/2.0);i++){
        
        u1[i]= rho;
        uhalf1[i]=0.0;
        u2[i]= rho*u;
        uhalf2[i]= 0.0;
        u3[i]= rho*et;
        uhalf3[i]=0.0;
        f1[i]= u2[i];
        f2[i]= u2[i]*u2[i]/u1[i] + (gamma - 1.0)*(u3[i] - 0.5*u2[i]*u2[i]/u1[i] );
        f3[i]= (u2[i]/u1[i])*(u3[i] + (gamma - 1.0)*(u3[i] - (0.5*u2[i]*u2[i]/u1[i]))  );
        
         //printf("%f %f %f \n",uhalf1[i], uhalf2[i],uhalf3[i]);
       // printf("%f %f %f \n",u1[i], u2[i],u3[i]);
        
    }
    
    rho=0.125;
    p=10000;
    u=0;
    et= p/((gamma-1.0)*rho) + pow(u,2)/2.0;
    for(i=(nx/2.0);i<nx;i++){
        
        u1[i]= rho;
        uhalf1[i]=0.0;
        u2[i]= rho*u;
        uhalf2[i]= 0.0;
        u3[i]= rho*et;
        uhalf3[i]= 0.0;
        f1[i]= u2[i];
        f2[i]= u2[i]*u2[i]/u1[i] + (gamma - 1.0)*(u3[i] - 0.5*u2[i]*u2[i]/u1[i] );
        f3[i]= (u2[i]/u1[i])*(u3[i] + (gamma - 1.0)*(u3[i] - (0.5*u2[i]*u2[i]/u1[i]))  );
        
        
        printf("%f %f %f \n",uhalf1[i], uhalf2[i],uhalf3[i]);
        printf("%f %f %f \n",u1[i], u2[i],u3[i]);


    }
    
    
    for(n=0;n<nt;n++){

        for(i=1;i<nx;i++){
            uhalf1[i]=(u1[i+1]+u1[i])/2.0-(f1[i+1]-f1[i])*dt/(2.0*dx);
            uhalf2[i]=(u2[i+1]+u2[i])/2.0-(f2[i+1]-f2[i])*dt/(2.0*dx);
            uhalf3[i]=(u3[i+1]+u3[i])/2.0-(f3[i+1]-f3[i])*dt/(2.0*dx);
           //  printf("%f %f %f \n",uhalf1[i], uhalf2[i],uhalf3[i]);
        }
    
        for(i=0;i<nx;i++){
            f1[i]=uhalf2[i];
            f2[i]=uhalf2[i]*uhalf2[i]/uhalf1[i]+(gamma-1)*(uhalf3[i]-uhalf2[i]*uhalf2[i]/(2.0*uhalf1[i]));
            f3[i]=( ( uhalf2[i] )/( uhalf1[i] ) )*( uhalf3[i] + ( gamma - 1.0 )*( uhalf3[i] - 0.5*( pow(uhalf2[i],2) )/( uhalf1[i] ) ) );;
             // printf("%f %f %f \n",f1[i], f2[i],f3[i]);
        }
    
        for(i=1;i<nx;i++)
        {
            u1[i]=u1[i]-(f1[i]-f1[i-1])*dt/(dx);
            u2[i]=u2[i]-(f2[i]-f2[i-1])*dt/(dx);
            u3[i]=u3[i]-(f3[i]-f3[i-1])*dt/(dx);
            //printf("%f %f %f \n",u1[i], u2[i],u3[i]);
        }
        
        
        for(i=0;i<nx;i++){
                    
            f1[i]= u2[i];
            f2[i]= u2[i]*u2[i]/u1[i] + (gamma - 1.0)*(u3[i] - 0.5*u2[i]*u2[i]/u1[i] );
            f3[i]=( ( u2[i] )/( u1[i] ) )*( u3[i] + ( gamma - 1.0 )*( u3[i] - 0.5*( pow(u2[i],2) )/( u1[i] ) ) );
            //        printf("%f %f %f \n",f1[i], f2[i],f3[i]);
                    
                }

    }
                
    /*escrbir el archivo*/
    
    for(i=0;i<nx;i++){
       
        x = xi + i*dx;
        rho = u1[i];
        u = u2[i]/rho;
        p=(gamma - 1.0)*(u3[i] - 0.5*(pow(u2[i],2))/(u1[i]));
    
        fprintf(in0,"%f %f %f %f\n",x,u,p,rho);

    }
    
}





