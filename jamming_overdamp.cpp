#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
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
#define dphi_a 0.002
#define dtbd 0.01
#define dim 2
#define polydispersity 0.0
#define skin 2.0

void ini_coord(double (*x)[dim],double L){
  //srand((unsigned int)time(NULL));
  for(int i=0;i<N1+N2;i++){
    x[i][0] = (double)rand()/RAND_MAX*L;
    x[i][1] = (double)rand()/RAND_MAX*L;
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
  *L*=sqrt(phi/(phi+dphi_a));
  for(i=0;i<N1+N2;i++){
    for(j=0;j<dim;j++){
      x[i][j] *= sqrt(phi/(phi+dphi_a));
    }
  }
}

void list_verlet(int (*list)[N1+N2],double (*x)[dim],double L){
  double dx,dy,dr2;
  int i,j;
  for(i=0;i<N1+N2;i++){
    for(j=0;j<N1+N2;j++){
      list[i][j]=0;
    }
  }
  
  for(i=0;i<N1+N2;i++){
    for(j=0;j<N1+N2;j++){
      if(i!=j){
        dx=x[i][0]-x[j][0];   // tagged particle is i
        dy=x[i][1]-x[j][1];
        dx-=L*floor((dx+0.5*L)/L);  // boundary condition
        dy-=L*floor((dy+0.5*L)/L);
        dr2=dx*dx+dy*dy;
        if(dr2<skin*skin){
          list[i][0]++;
          list[i][(int)list[i][0]]=j;
        }
      }
    }
  }
}

void update(double (*x_update)[dim],double (*x)[dim]){
  for(int i=0;i<N1+N2;i++)
    for(int j=0;j<dim;j++)
      x_update[i][j]=x[i][j];
}

void calc_disp_max(double *disp_max,double (*x)[dim],double (*x_update)[dim],double L){
  double dx,dy;
  double disp;
  for(int i=0;i<N1+N2;i++){
    dx=x[i][0]-x_update[i][0];
    dy=x[i][1]-x_update[i][1];
    dx-=L*floor((dx+0.5*L)/L);
    dy-=L*floor((dy+0.5*L)/L);
    disp = dx*dx+dy*dy;
    if(disp > *disp_max){
      *disp_max =disp;
    }
  }
  //std::cout<<"disp_max is"<<"\t"<<*disp_max<<"\n";
}

void auto_list_update(double *disp_max,double (*x)[dim],double (*x_update)[dim],int (*list)[N1+N2],double L){
  static int count=0;
  count++;
  calc_disp_max(&(*disp_max),x,x_update,L);
  if(*disp_max > skin*skin*0.25){
    list_verlet(list,x,L);
    update(x_update,x);
    *disp_max=0.0;
        //std::cout<<"update"<<*disp_max<<" "<<count<<std::endl;
    count=0;
  }
} 

void calc_force(int (*list)[N1+N2],double (*x)[dim],double (*f)[dim],double (*fij)[dim],double *a,double *U,double *P,double *Fmax,double L){
  double dx,dy,dr,dr2,dUr,aij,F,overlap;
  int i,j;

  *U=0.0;
  *P=0.0;
  *Fmax=0.0;
  overlap=0.0;
  ini_array(f);
  
  // calculate force
  for(i=0;i<N1+N2;i++){
    for(j=1;j<=list[i][0];j++){
      //if(i!=j){
        dx=x[i][0]-x[list[i][j]][0];   // tagged particle is i
        dy=x[i][1]-x[list[i][j]][1];
        dx-=L*floor((dx+0.5*L)/L);  // boundary condition
        dy-=L*floor((dy+0.5*L)/L);
        dr2=dx*dx+dy*dy;
        dr=sqrt(dr2);
        aij=0.5*(a[i]+a[list[i][j]]);
        //std::cout << i << "\t" << j << "\t" << dx/dr << "\t" << dy/dr << "\t" << dUr << "\n";  // cout
        if(dr<aij){
          dUr=-(1.0-dr/aij)/aij;  // <0
          f[i][0]-=dUr*dx/dr;
          fij[list[i][j]][0]=-dUr*dx/dr;
          f[i][1]-=dUr*dy/dr;
          fij[list[i][j]][1]=-dUr*dy/dr;
          if(aij-dr>overlap){
            overlap=aij-dr;
          }
          if(i<j){
            *U+=(1.0-dr/aij)*(1.0-dr/aij)/2.0/(double)(N1+N2);  // alpha=2 at Hertzian potential
            *P-=dr*dUr/(2.0*L*L);
          }
        std::cout << i << "\t" << list[i][j] << "\t" << dx/dr << "\t" << dy/dr << "\t" << dUr << "\n";  // cout
        }
      //}
    }
  }

  // evaluate power
  for(i=0;i<N1+N2;i++){
    //std::cout<<i<<"\t"<<x[i][0]<<"\t"<<x[i][1]<<"\n";
    F=sqrt(f[i][0]*f[i][0]+f[i][1]*f[i][1]);
    //std::cout<<i<<"\t"<<F<<"\n";
    if(F>*Fmax){
      *Fmax=F;
    }
  }
  std::cout << *Fmax << "\t" << overlap << "\t" <<  *U << "\t" << *P << "\n";  // cout
}

void eom_overdamp(int (*list)[N1+N2],double (*x)[dim],double (*x_update)[dim],double (*f)[dim],double (*fij)[dim],double *a,double *U,double *P,double *Fmax,double *disp_max,double dt,double L){
  calc_force(list,x,f,fij,a,&(*U),&(*P),&(*Fmax),L);
  //std::cout << *U << std::endl;  //cout
  for(int i=0;i<N1+N2;i++){
    for(int j=0;j<dim;j++){
      x[i][j]+=f[i][j]*dt;
    }
  }
  p_boundary(x,L);
  auto_list_update(disp_max,x,x_update,list,L);
}

void output_coord(double (*x)[dim],double *a,double phi){
  char filename[128];
  std::ofstream file;
  sprintf(filename,"coord_od_phi%.3f_debug.dat",phi);
  file.open(filename);
  for(int i=0;i<N1+N2;i++)
    file <<x[i][0]<<"\t"<<x[i][1]<<"\t"<<a[i]<<std::endl;
  file.close();
}

void output_dr(double (*x)[dim],double *a,double L,double phi){
  double dx,dy,dr,aij;
  int i,j;
  int k=0;
  char filename[128];
  std::ofstream file;
  sprintf(filename,"dr_od_phi%.3f_debug.dat",phi);
  file.open(filename);
  for(i=0;i<N1+N2;i++){
    for(j=0;j<i;j++){
      dx=x[i][0]-x[j][0];
      dy=x[i][1]-x[j][1];
      dx-=L*floor((dx+0.5*L)/L);  // boundary condition
      dy-=L*floor((dy+0.5*L)/L);
      dr=sqrt(dx*dx+dy*dy);
      aij=0.5*(a[i]+a[j]);
      if(dr<aij){
        file << i << "\t" << j << "\t" << aij-dr <<std::endl;
        k++;
      }
    }
  }
  file <<"Number of overlapping particles is"<<"\t"<<k<<std::endl;
  file.close();
}

void output_potential(double U,double phi){
  char filename[128];
  std::ofstream file;
  sprintf(filename,"potential_od_debug.dat");
  file.open(filename,std::ios::app);
  file << phi << "\t" << U << std::endl;
  file.close();
}

void output_pressure(double P,double phi){
  char filename[128];
  std::ofstream file;
  sprintf(filename,"pressure_od_debug.dat");
  file.open(filename,std::ios::app);
  file << phi << "\t" << P << std::endl;
  file.close();
}

void initialize(int (*list)[N1+N2],double (*x)[dim],double (*f)[dim],double (*fij)[dim],double *a,double *U,double *P,double *Fmax,double *phi,double L){
  ini_array(x);
  ini_coord(x,L);
  *phi=phi0;
  list_verlet(list,x,L);
  calc_force(list,x,f,fij,a,&(*U),&(*P),&(*Fmax),L);
  *P=0.0;
  *U=0.0;
}

int main(){
  double x[N1+N2][dim],x_update[N1+N2][dim],f[N1+N2][dim],fij[N1*N2][dim],a[N1+N2],U,P,phi,F,Fmax;
  double L=sqrt((N1*a1/2.0*a1/2.0+N2*a2/2.0*a2/2.0)*M_PI/phi0),disp_max=0.0;
  char filename[128];
  int list[N1+N2][N1+N2];

  set_diameter(a);
  initialize(list,x,f,fij,a,&U,&P,&Fmax,&phi,L);
  //std::cout<<"initialized"<<"\t"<<"\t"<<"P="<<P<<std::endl;
  while(phi<0.9){
    while(Fmax>1.0e-8){
      eom_overdamp(list,x,x_update,f,fij,a,&U,&P,&Fmax,&disp_max,dtbd,L);
    }
    //output_coord(x,a,phi);
    //output_dr(x,a,L,phi);
    //output_potential(U,phi);
    //output_pressure(P,phi);
    //std::cout<<"phi="<<phi<<"\t"<<"U="<<U<<"\t"<<"P="<<P<<std::endl;
    affine_transformation(x,&L,phi);
    phi+=dphi_a;
    calc_force(list,x,f,fij,a,&U,&P,&Fmax,L);
  }

  return 0;
}
