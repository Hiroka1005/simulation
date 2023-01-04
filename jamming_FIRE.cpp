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
#define phi_max 0.86
#define dphi_fixed 1.0e-2
#define dphi_a 2.0e-3
#define dphi_g 1.0e-5
#define dt0 0.01
#define dt_max 0.1
#define zeta 1.0
#define dim 2
#define polydispersity 0.0
#define cut 1.4
#define skin 1.0

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

void affine_transformation(double (*x)[dim],double *L,double *phi,double dphi){
  int i,j=0;
  // update phi
  *phi+=dphi;
  // transform L and x by using new phi
  *L*=sqrt((*phi-dphi)/(*phi));
  for(i=0;i<N1+N2;i++){
    for(j=0;j<dim;j++){
      x[i][j] *= sqrt((*phi-dphi)/(*phi));
    }
  }
}

void affine_geometric(double (*x)[dim],double *L,double *phi){
  int i,j=0;
  // update phi
  *phi*=(1.0+dphi_g);
  // transform L and x by using old phi
  *L*=sqrt(1.0/(1.0+dphi_g));
  for(i=0;i<N1+N2;i++){
    for(j=0;j<dim;j++){
      x[i][j]*=sqrt(1.0/(1.0+dphi_g));
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
        dx=x[i][0]-x[j][0];
        dy=x[i][1]-x[j][1];
        dx-=L*floor((dx+0.5*L)/L);
        dy-=L*floor((dy+0.5*L)/L);
        dr2=dx*dx+dy*dy;
        if(dr2<(cut+skin)*(cut+skin)){
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
    //std::cout<<"update"<<*disp_max<<" "<<count<<std::endl;
    *disp_max=0.0;
    count=0;
  }
} 

void calc_force(int (*list)[N1+N2],double (*x)[dim],double (*f)[dim],double *a,double *U,double *p,double L){
  double dx,dy,dr,dr2,dUr,aij,overlap;
  int i,j;

  *U=0.0;
  *p=0.0;
  overlap=0.0;
  ini_array(f);
  
  // calculate force
  for(i=0;i<N1+N2;i++){
    for(j=1;j<=list[i][0];j++){
      dx=x[i][0]-x[list[i][j]][0];
      dy=x[i][1]-x[list[i][j]][1];
      dx-=L*floor((dx+0.5*L)/L);
      dy-=L*floor((dy+0.5*L)/L);
      dr2=dx*dx+dy*dy;
      dr=sqrt(dr2);
      aij=0.5*(a[i]+a[list[i][j]]);
      if(dr<aij){
        dUr=-(1.0-dr/aij)/aij;
        f[i][0]-=dUr*dx/dr;
        f[list[i][j]][0]+=dUr*dx/dr;
        f[i][1]-=dUr*dy/dr;
        f[list[i][j]][1]+=dUr*dy/dr;
        if(aij-dr>overlap){
          overlap=aij-dr;
        }
        if(i<list[i][j]){
          *U+=(1.0-dr/aij)*(1.0-dr/aij)/2.0/(double)(N1+N2);  // alpha=2 at Hertzian potential
          *p-=dr*dUr/(2.0*L*L);
        }
      }
    }
  }
}

void FIRE(int (*list)[N1+N2],double (*x)[dim],double (*x_update)[dim],double (*f)[dim],double *a,double *U,double *p,double *f_tot,double *disp_max,double L){
  double v[N1+N2][dim],P,v_norm,f_norm;
  double dt=dt0,A=0.1;
  int count=0;
  
  ini_array(v);
  list_verlet(list,x,L);
  calc_force(list,x,f,a,&(*U),&(*p),L);

  for(;;){
    // initialize
    P=0.0;
    *f_tot=0.0;
    // 1. MD integration with Velocity Verlet
    for(int i=0;i<N1+N2;i++){
      for(int j=0;j<dim;j++){
        x[i][j]+=v[i][j]*dt+0.5*f[i][j]*dt*dt;
        v[i][j]+=0.5*f[i][j]*dt;
      }
    }
    calc_force(list,x,f,a,&(*U),&(*p),L);
    for(int i=0;i<N1+N2;i++){
      for(int j=0;j<dim;j++){
        v[i][j]+=0.5*f[i][j]*dt;
      }
    }
    p_boundary(x,L);
    auto_list_update(&(*disp_max),x,x_update,list,L);

    // 3. calculate power P
    for(int i=0;i<N1+N2;i++){
      v_norm=sqrt(v[i][0]*v[i][0]+v[i][1]*v[i][1]);
      f_norm=sqrt(f[i][0]*f[i][0]+f[i][1]*f[i][1]);
      *f_tot+=f_norm/(N1+N2);
      P+=f[i][0]*v[i][0]+f[i][1]*v[i][1];
      // 4. velocity modification
      for(int j=0;j<dim;j++){
        v[i][j]=(1.0-A)*v[i][j]+A*f[i][j]*v_norm/(f_norm+DBL_EPSILON);
      }
    }

    // converge criterion
    if(*f_tot<1.0e-12){
      break;
    }

    // 5. evaluate power P
    if(P>=0.0){
      // 6. P>0
      count++;
      if(count>5){
        dt=std::min(1.1*dt,dt_max);
        A*=0.99;
      }
    }else{
      // 7. P<0
      count=0;
      ini_array(v);
      dt*=0.5;
      A=0.1;
    }

  }
}

void output_coord(int step,double (*x)[dim],double *a,double phi){
  char filename[128];
  std::ofstream file;
  sprintf(filename,"coord_FIRE_step%d_phi%.6f.dat",step,phi);
  file.open(filename);
  for(int i=0;i<N1+N2;i++)
    file <<x[i][0]<<"\t"<<x[i][1]<<"\t"<<a[i]<<std::endl;
  file.close();
}

void output_potential(int step,double U,double phi){
  char filename[128];
  std::ofstream file;
  sprintf(filename,"potential_FIRE_step%d.dat",step);
  file.open(filename,std::ios::app);
  file << phi << "\t" << U << std::endl;
  file.close();
}

void output_pressure(int step,double p,double phi){
  char filename[128];
  std::ofstream file;
  sprintf(filename,"pressure_FIRE_step%d.dat",step);
  file.open(filename,std::ios::app);
  file << phi << "\t" << p << std::endl;
  file.close();
}

void initialize(int (*list)[N1+N2],double (*x)[dim],double (*f)[dim],double *a,double *U,double *p,double *phi,double L){
  *phi=phi0;
  ini_array(x);
  ini_coord(x,L);
  list_verlet(list,x,L);
  calc_force(list,x,f,a,&(*U),&(*p),L);
}

void training(int (*list)[N1+N2],int step,double (*x)[dim],double (*x_update)[dim],double (*f)[dim],double *a,double *U,double *p,double *f_tot,double *phi,double dphi,double *disp_max,double *L){
  affine_transformation(x,&(*L),&(*phi),dphi);
  FIRE(list,x,x_update,f,a,&(*U),&(*p),&(*disp_max),&(*f_tot),*L);

  //output_coord(step,x,a,*phi);
  output_potential(step,*U,*phi);
  output_pressure(step,*p,*phi);
}

int main(){
  double x[N1+N2][dim],x_update[N1+N2][dim],v[N1+N2][dim],f[N1+N2][dim],a[N1+N2],U,p,phi,phi_J,f_tot;
  double L=sqrt((N1*a1/2.0*a1/2.0+N2*a2/2.0*a2/2.0)*M_PI/phi0),disp_max=0.0;
  char filename[128];
  int list[N1+N2][N1+N2],step;

  set_diameter(a);
  initialize(list,x,f,a,&U,&p,&phi,L);

  step=1;
  while(phi<=phi_max){
    training(list,step,x,x_update,f,a,&U,&p,&f_tot,&phi,dphi_a,&disp_max,&L);
  }

  step=2;
  while(p>1.0e-3){
    training(list,step,x,x_update,f,a,&U,&p,&f_tot,&phi,-dphi_a*5.0e-1,&disp_max,&L);
  }
  while(p>1.0e-4){
    training(list,step,x,x_update,f,a,&U,&p,&f_tot,&phi,-dphi_a*5.0e-2,&disp_max,&L);
  }
  while(p>1.0e-5){
    training(list,step,x,x_update,f,a,&U,&p,&f_tot,&phi,-dphi_a*5.0e-3,&disp_max,&L);
  }
  while(p>1.0e-6){
    training(list,step,x,x_update,f,a,&U,&p,&f_tot,&phi,-dphi_a*5.0e-4,&disp_max,&L);
  }
  while(p>1.0e-7){
    training(list,step,x,x_update,f,a,&U,&p,&f_tot,&phi,-dphi_a*5.0e-5,&disp_max,&L);
  }
  while(p>1.0e-8){
    training(list,step,x,x_update,f,a,&U,&p,&f_tot,&phi,-dphi_a*5.0e-6,&disp_max,&L);
  }
  while(p>1.0e-9){
    training(list,step,x,x_update,f,a,&U,&p,&f_tot,&phi,-dphi_a*5.0e-7,&disp_max,&L);
  }
  phi_J=phi;

  step=3;
  while(phi<=phi_J+dphi_fixed){
    affine_geometric(x,&L,&phi);
    FIRE(list,x,x_update,f,a,&U,&p,&f_tot,&disp_max,L);

    //output_coord(step,x,a,phi);
    output_potential(step,U,phi);
    output_pressure(step,p,phi);
  }

  return 0;
}
