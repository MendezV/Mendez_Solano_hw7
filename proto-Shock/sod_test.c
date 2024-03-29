

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
    


    /*condiciones iniciales*/
  
    for(i=0;i<(nx);i++){
        
        if(i<nx*0.5){
        rho=1;
        p=100000.0;
        u=0;
        et= p/((gamma-1.0)*rho) + (u*u)/2.0;
        }
        else{
        rho=0.125;
        p=10000.0;
        u=0;
        et= p/((gamma-1.0)*rho) + (u*u)/2.0;
        }
        
        u1[i]= rho;
        uhalf1[i]=0.0;
        u2[i]= rho*u;
        uhalf2[i]= 0.0;
        u3[i]= rho*et;
        uhalf3[i]= 0.0;
        
    }
    
   
    /*evolucion temporal*/
        
    for(n=0;n<nt;n++){
        uhalf1[0] = nasty_half_vect(u1[1],u1[0],fvect1(u2[1]),fvect1(u2[0]),dx,dt);
        uhalf2[0] = nasty_half_vect(u2[1],u2[0],fvect2(u1[1],u2[1],u3[1],gamma),fvect2(u1[0],u2[0],u3[0],gamma),dx,dt);
        uhalf3[0] = nasty_half_vect(u3[1],u3[0],fvect3(u1[1],u2[1],u3[1],gamma),fvect3(u1[0],u2[0],u3[0],gamma),dx,dt);
      
          for(i=1;i<(nx-1.0);i++){
 
              
              uhalf1[i] = nasty_half_vect(u1[i+1],u1[i],fvect1(u2[i+1]),fvect1(u2[i]),dx,dt);
              uhalf2[i] = nasty_half_vect(u2[i+1],u2[i],fvect2(u1[i+1],u2[i+1],u3[i+1],gamma),fvect2(u1[i],u2[i],u3[i],gamma),dx,dt);
              uhalf3[i] = nasty_half_vect(u3[i+1],u3[i],fvect3(u1[i+1],u2[i+1],u3[i+1],gamma),fvect3(u1[i],u2[i],u3[i],gamma),dx,dt);
              
              
              u1[i] = nasty_one_vect(u1[i],fvect1(uhalf2[i]),fvect1(uhalf2[i-1]),dx,dt);
              u2[i] = nasty_one_vect(u2[i],fvect2(uhalf1[i],uhalf2[i],uhalf3[i],gamma),fvect2(uhalf1[i-1],uhalf2[i-1],uhalf3[i-1],gamma),dx,dt);
              u3[i] = nasty_one_vect(u3[i],fvect3(uhalf1[i],uhalf2[i],uhalf3[i],gamma),fvect3(uhalf1[i-1],uhalf2[i-1],uhalf3[i-1],gamma),dx,dt);
              
          
        }
    
    }
        


    /*escrbir el archivo*/
    
    for(i=0;i<nx;i++){
       
        x = xi + i*dx;
        rho = u1[i];
        u = u2[i]/rho;
        p=(gamma-1)*(u3[i]-u2[i]*u2[i]/(2.0*u1[i]));
    
        fprintf(in0,"%f %f %f %f\n",x,u,p,rho);

    }
    
}

double nasty_half_vect(double uione,double ui,double fione,double fi,double dx,double dt){
    
   double uithalf= 0.5*((uione + ui) - (fione - fi)*(dt/dx)) ;
    return uithalf;
}


double nasty_one_vect(double ui,double fihalf,double fiminhalf,double dx,double dt){
    
   double utione= ui - (fihalf - fiminhalf)*(dt/dx) ;
    return utione;
}

double fvect1(double u2){
    
    return u2;
}

double fvect2(double u1,double u2,double u3,double gamma){
    
   double f2=  (u2*u2)/u1 + (gamma - 1.0)*(u3 - ((u2*u2)/(2.0*u1)) );
    return f2;
    
}

double fvect3(double u1,double u2,double u3,double gamma){

   double f3= (u2/u1)*(u3 + (gamma - 1.0)*(u3 - (((u2*u2)/(2.0*u1))))  );
    return f3;
    
}



