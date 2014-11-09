#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//double e_T(double e,double u);
double energia(double p,double rho,double gamma);
double presion(double u_1,double u_2,double u_3,double gamma);
double paso_8(double u_1,double u_2,double f_1,double f_2,double dx,double dt);
double paso_9(double u_1,double f_1,double f_2,double dx,double dt);
double f_1(double u_2);
double f_2(double u_1,double u_2,double u_3,double gamma);
double f_3(double u_1,double u_2,double u_3,double gamma);



int main(int argc, char **argv){
    
    
    // condiciones iniciales
    double p_I = 100000.0,p_D = 10000.0,rho_I = 1.0,rho_D = 0.125,gamma = 1.4;
    
    //tiempo final
    char *t_string = argv[1];
    double t = atof(argv[1]);
    
    int N_t,N_x,i,j;
    double dt,dx, *U_n_x, *U_nh_x, *U_n_y, *U_nh_y, *U_n_z, *U_nh_z;
    //double  *U_n1_x, *U_n1_z, *U_n1_y;
    
    N_x = 1000;
    N_t = 1000;
    dt = t/N_t;
    dx = 20.0/(N_x);
    U_n_x = malloc(N_x * sizeof(double));
    U_nh_x = malloc(N_x * sizeof(double));
    //U_n1_x = malloc(N_x * sizeof(double));
    U_n_y= malloc(N_x * sizeof(double));
    U_nh_y = malloc(N_x * sizeof(double));
    //U_n1_y = malloc(N_x * sizeof(double));
    U_n_z = malloc(N_x * sizeof(double));
    U_nh_z = malloc(N_x * sizeof(double));
    //U_n1_z = malloc(N_x * sizeof(double));
    
    //Se ponen las condiciones iniciales.
    for(i=0; i<N_x/2; i++){
        U_n_y[i] = 0.0; U_nh_y[i] = 0.0; //U_n1_y[i] = 0.0;
        U_n_x[i] = rho_I; U_nh_x[i] = 0.0; //U_n1_x[i] = 0.0;
        U_n_z[i] = (rho_I)*( energia(p_I,rho_I,gamma) ); U_nh_z[i] = 0.0; //U_n1_z[i] = 0.0;
        printf("%f %f %f\n", U_nh_x[i],U_nh_y[i],U_nh_z[i] );
        printf("%f %f %f\n", U_n_x[i-1],U_n_y[i-1],U_n_z[i-1] );
        
    }
    
    for(i=N_x/2; i<N_x; i++){
        U_n_y[i] = 0.0; U_nh_y[i] = 0.0;             //U_n1_y[i] = 0.0;
        U_n_x[i] = rho_D; U_nh_x[i] = 0.0;          //U_nh_x[i] = 0.0;
        U_n_z[i] = (rho_D)*( energia(p_D,rho_D,gamma) ); U_nh_z[i] = 0.0; //U_nh_z[i] = 0.0;
        printf("%f %f %f\n", U_nh_x[i],U_nh_y[i],U_nh_z[i] );
        printf("%f %f %f\n", U_n_x[i-1],U_n_y[i-1],U_n_z[i-1] );
    }
    
    
    // Metodo finito
    //for(j=0;j<N_t;j++){
        
        for(i=1;i<N_x;i++){
            
            U_nh_x[i] = paso_8(U_n_x[i],U_n_x[i-1],f_1(U_n_y[i]),f_1(U_n_y[i-1]),dx,dt);
            U_nh_y[i] = paso_8(U_n_y[i],U_n_y[i-1],f_2(U_n_x[i],U_n_y[i],U_n_z[i],gamma),f_2(U_n_x[i-1],U_n_y[i-1],U_n_z[i-1],gamma),dx,dt);
            U_nh_z[i] = paso_8(U_n_z[i],U_n_z[i-1],f_3(U_n_x[i],U_n_y[i],U_n_z[i],gamma),f_3(U_n_x[i-1],U_n_y[i-1],U_n_z[i-1],gamma),dx,dt);
            
            //printf("%f %f %f\n", U_nh_x[i-1],U_nh_y[i-1],U_nh_z[i-1] );
            
            
            
        }
        
        for(i=1;i<N_x-1;i++){
            U_n_x[i] = paso_9(U_n_x[i],f_1(U_nh_y[i+1]),f_1(U_nh_y[i]),dx,dt);
            U_n_y[i] = paso_9(U_n_y[i],f_2(U_nh_x[i+1],U_nh_y[i+1],U_nh_z[i+1],gamma),f_2(U_nh_x[i],U_nh_y[i],U_nh_z[i],gamma),dx,dt);
            U_n_z[i] = paso_9(U_n_z[i],f_3(U_nh_x[i+1],U_nh_y[i+1],U_nh_z[i+1],gamma),f_3(U_nh_x[i],U_nh_y[i],U_nh_z[i],gamma),dx,dt);
            
            //printf("%f %f %f\n", U_n_x[i-1],U_n_y[i-1],U_n_z[i-1] );
            
            U_n_x[0] = U_n_x[1];
            U_n_y[0] = U_n_y[1];
            U_n_z[0] = U_n_z[1];
            
            U_n_x[N_x-1] = U_n_x[N_x-2];
            U_n_y[N_x-1] = U_n_y[N_x-2];
            U_n_z[N_x-1] = U_n_z[N_x-2];
            
        }
        
   // }
    
    
    
    double p,rho_f,u,x;
    FILE *fileout;
    
    char filename[100];
    sprintf(filename, "estado_%s.dat", t_string);
    fileout= fopen(filename, "w");
    
    for(i=0;i<N_x;i++){
        
        x = -10 + i*dx;
        
        p = presion(U_n_x[i],U_n_y[i],U_n_z[i],gamma);
        rho_f = U_n_x[i];
        u = ( U_n_y[i] )/( rho_f );
        
        fprintf(fileout,"%f %f %f %f\n",x,u,p,rho_f);
        
        
        //fprintf(fileout,"%f %f %f %f\n",x,U_n_x[i],U_n_y[i],U_n_z[i]);
        
    }
    
    fclose(fileout);
    
    
    return 0;
    
}



//energia total
/*double e_T(double e,double u){
 double respuesta = e + 0.5*( pow(u,2) );
 return respuesta;
 }*/

// energia

double energia(double p,double rho,double gamma){
    
    double respuesta = p/( ( gamma - 1.0 )*rho );
    
    return respuesta;
}
// calcula la presion en un puntos
double presion(double u_1,double u_2,double u_3,double gamma){
    
    double respuesta = ( gamma - 1.0 )*( u_3 - 0.5*( pow(u_2,2) )/(u_1) );
    
    return respuesta;
}


//PASO 8: mucho cuidado aca, u_1 no es la primera componente, solo hace referencia a la primera velocidad que sale en la ecuación 8 del enlace recomendado para hacer el ejercicio. Lo mismo aplica para el resto. (ej: u_1 = u^n_{i+1} )
double paso_8(double u_1,double u_2,double f_1,double f_2,double dx,double dt){
    
    double respuesta = 0.5*( u_1 + u_2 ) - 0.5*( f_1 - f_2 )*(dt/dx) ;
    
    return respuesta;
}

//PASO 9 :mismo cuidado que el paso_8
double paso_9(double u_1,double f_1,double f_2,double dx,double dt){
    
    double respuesta =  u_1 - ( f_1 - f_2 )*(dt/dx) ;
    
    return respuesta;
}

// primera componente del vector f
double f_1(double u_2){
    
    return u_2;
}

// segunda componente del vector f
double f_2(double u_1,double u_2,double u_3,double gamma){
    
    double respuesta = ( pow(u_2,2) )/( u_1 ) + ( gamma - 1.0 )*( u_3 - 0.5*( pow(u_2,2) )/(u_1) );
    
    return respuesta;
    
}

// tercera componente del vector f
double f_3(double u_1,double u_2,double u_3,double gamma){
    
    double respuesta = ( ( u_2 )/( u_1 ) )*( u_3 + ( gamma - 1.0 )*( u_3 - 0.5*( pow(u_2,2) )/( u_1 ) ) );
    
    return respuesta;
    
}

