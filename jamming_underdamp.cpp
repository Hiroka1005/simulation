#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cfloat>
#include <algorithm>
#include "BM.h"

#define a1 1.0
#define a2 1.4
#define N1 162
#define N2 162
#define phi0 0.80
#define dtbd 0.01
#define dim 2
#define polydispersity 0.0

void ini_coord(double (*x)[dim],double L){
  for(int i=0;i<N1+N2;i++){
    x[i][0] = (double)rand()/RAND_MAX*L;
    x[i][1] = (double)rand()/RAND_MAX*L;
    // std::cout << x[i][0] << "\t" << x[i][1] << std::endl;
  }
}

void set_diameter(double *a){
  for(int i=0;i<N1;i++){
    a[i]=a1+polydispersity*gaussian_rand();
  }
  for(int i=0;i<N2;i++){
    a[i+N1]=a2+polydispersity*gaussian_rand();
  }
}

void p_boundary(double (*x)[dim],double L){
  for(int i=0;i<N1+N2;i++)
    for(int j=0;j<dim;j++)
      x[i][j]-=L*floor(x[i][j]/L);
}

void ini_array(double (*x)[dim]){
  for(int i=0;i<N1+N2;i++)
    for(int j=0;j<dim;j++)
      x[i][j]=0.0;
}

void affine_transformation(double (*x)[dim],double *L,double phi){
  int i,j=0;
  *L*=sqrt(phi/(phi+0.01));
  for(i=0;i<N1+N2;i++){
    for(j=0;j<dim;j++){
      x[i][j] *= sqrt(phi/(phi+0.01));
    }
  }
}

void calc_force(int (*list)[N1+N2],double (*x)[dim],double (*f)[dim],double *a,double *U,double *P,double *Fmax,double L){
  double dx,dy,dr,dr2,dUr,aij,F;
  int i,j;
  *U=0.0;
  *Fmax=0.0;
  ini_array(f);
  for(i=0;i<N1+N2;i++){
    for(j=0;j<N1+N2;j++){
      list[i][j]=0;
    }
  }
  
  for(i=0;i<N1+N2;i++){
    for(j=0;j<N1+N2;j++){
      if(i<j){
        dx=x[i][0]-x[j][0];
        dy=x[i][1]-x[j][1];
        dx-=L*floor((dx+0.5*L)/L);  // boundary condition
        dy-=L*floor((dy+0.5*L)/L);
        dr2=dx*dx+dy*dy;
        dr=sqrt(dr2);
        aij=0.5*(a[i]+a[j]);

        if(dr<aij){
          dUr=(dr-aij)/(aij*aij);
          f[i][0]-=dUr*dx/dr;
          f[list[i][j]][0]+=dUr*dx/dr;
          f[i][1]-=dUr*dy/dr;
          f[list[i][j]][1]+=dUr*dy/dr;
          *U+=(dr-aij)*(dr-aij)/(2.0*aij*aij);  // alpha=2 at Hertzian potential
          *P-=dUr/(2.0*L*L);
          //std::cout << i << "\t" << j << "\t" << sqrt(dx*dx+dy*dy)-(0.5*(a[i]+a[j])) << "\n";
        }else{
          //std::cout << i << "\t" << j << "\t" << 0.0 << "\n";
        }

      }
    }
  }

  for(int i=0;i<N1+N2;i++){
    F=sqrt(f[i][0]*f[i][0]+f[i][1]*f[i][1]);
    if(F>*Fmax){
      *Fmax=F;
    }
  }
}

void eom_underdamp(int (*list)[N1+N2],double (*v)[dim],double (*x)[dim],double (*f)[dim],double *a,double *U,double *P,double *Fmax,double dt,double L){
  double zeta=1.0;
  
  calc_force(list,x,f,a,&(*U),&(*P),Fmax,L);
  //std::cout << *U << std::endl;  //cout
  for(int i=0;i<N1+N2;i++)
    for(int j=0;j<dim;j++){
      v[i][j]=(1.0-dt)*v[i][j]+f[i][j]*dt;
      x[i][j]+=v[i][j]*dt;
    }
  p_boundary(x,L);
}

void output_coord(double (*x)[dim],double *a,double phi){
  char filename[128];
  std::ofstream file;
  sprintf(filename,"coord_ud_phi%.2f.dat",phi);
  file.open(filename);
  for(int i=0;i<N1+N2;i++)
    file <<x[i][0]<<"\t"<<x[i][1]<<"\t"<<a[i]<<std::endl;
  file.close();
}

void output_potential(double U,double phi){
  char filename[128];
  std::ofstream file;
  sprintf(filename,"potential_ud.dat");
  file.open(filename,std::ios::app);
  file << phi << "\t" << U/(N1+N2) << std::endl;
  file.close();
}

void output_pressure(double P,double phi){
  char filename[128];
  std::ofstream file;
  sprintf(filename,"pressure_ud.dat");
  file.open(filename,std::ios::app);
  file << phi << "\t" << P << std::endl;
  file.close();
}

int main(){
  double x[N1+N2][dim],v[N1+N2][dim],f[N1+N2][dim],a[N1+N2],U,P,phi,F,Fmax;
  double L=sqrt((N1*a1/2.0*a1/2.0+N2*a2/2.0*a2/2.0)*M_PI/phi0);
  int list[N1+N2][N1+N2];
  char filename[128];
  set_diameter(a);
  ini_array(x);
  ini_array(v);
  ini_coord(x,L);
  
  phi=phi0;
  output_coord(x,a,phi+10.0);    // "+10": avoid duplication of filename
  calc_force(list,x,f,a,&U,&P,&Fmax,L);
  
  while(phi<0.9){
    while(Fmax>1.0e-8){
      eom_underdamp(list,v,x,f,a,&U,&P,&Fmax,dtbd,L);
    }
    //std::cout << "after calculating U=" << U << std::endl;  //cout

    output_coord(x,a,phi);
    output_potential(U,phi);
    output_pressure(P,phi);
    affine_transformation(x,&L,phi);
    phi=(N1*a1/2.0*a1/2.0+N2*a2/2.0*a2/2.0)*M_PI/(L*L);
    }

  return 0;
}
