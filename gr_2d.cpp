#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cfloat>
#include <vector>
//using namespace std;

#define a1 1.0
#define a2 1.4
#define N1 512
#define N2 512
#define dphi 1.0e-2
#define dr 1.0e-2
#define rmax 5.00
#define thresh 1.0e-3
#define nx_max 128
#define dim 2
#define sample 300

void histogram(double (*x)[dim+1],int (*Z)[2],double *Z0,double *Z1,double *Z2,double *Z3,double *Z4,double *Z5,double *Z6,double *Zave,double *gr,double *gr_SS,double *gr_LL,double *gr_SL,double L){
  double rij,aij,dx,dy,sum;
  double r=0.0;
  int i,j,k;

  //initialize
  //*Zave=0.0;
  //*Z1=0.0;
  for(i=0;i<N1+N2;i++){
    Z[i][0]=i;
    Z[i][1]=0;
  }

  for(i=0;i<N1+N2;i++){
    for(j=0;j<N1+N2;j++){
      if(i!=j){
        //calculate rij
        dx=x[i][0]-x[j][0];
        dy=x[i][1]-x[j][1];
        //std::cout << i << "\t" << j << "\t" << x[i][1] << "\t" << x[j][1] << "\t" << x[i][1]-x[j][1] << std::endl;
        dx-=L*floor((dx+0.5*L)/L);
        dy-=L*floor((dy+0.5*L)/L);
        rij=sqrt(dx*dx+dy*dy);
        aij=0.5*(x[i][2]+x[j][2]);

        // make histogram of \delta Z
        if(rij-aij<=thresh){
          Z[i][1]++;
          //std::cout << i << "\t" << j << "\t" << Z[i][1] << std::endl;
          *Zave+=1.0/(N1+N2)/sample;
        }
        
        //make histogram of g(r)
        if(rij<=rmax){
          k=(int)floor(rij*rmax/dr);
          r=k*dr;
          gr[k]+=1.0/(N1+N2)/(2.0*M_PI*(r+dr*1.0e-4)*dr*((N1+N2)/(L*L)))/sample;
          if(i<N1 && j<N1){
            gr_SS[k]+=1.0/(N1+N2)/(2.0*M_PI*(r+dr*1.0e-4)*dr*((N1+N2)/(L*L)))/sample;
          }else if(i>=N1 && j>=N1){
            gr_LL[k]+=1.0/(N1+N2)/(2.0*M_PI*(r+dr*1.0e-4)*dr*((N1+N2)/(L*L)))/sample;
          }else if(i<N1 && j>=N1){
            gr_SL[k]+=1.0/(N1+N2)/(2.0*M_PI*(r+dr*1.0e-4)*dr*((N1+N2)/(L*L)))/sample;
          }else if(i>=N1 && j<N1){
            gr_SL[k]+=1.0/(N1+N2)/(2.0*M_PI*(r+dr*1.0e-4)*dr*((N1+N2)/(L*L)))/sample;
          }
          goto next;
        }
        next:
        //std::cout<<r<<std::endl;
        r=0.0;
      }
    }
  }

  //calculate rattlers
  for(i=0;i<N1+N2;i++){
    if(Z[i][1]==0){
      *Z0+=1.0/sample;
      //std::cout << i << "\t" << Z[i][1] << "\t" << "Z0" << std::endl;
    }else if(Z[i][1]==1){
      *Z1+=1.0/sample;
      //std::cout << i << "\t" << Z[i][1] << "\t" << "Z1" << std::endl;
    }else if(Z[i][1]==2){
      *Z2+=1.0/sample;
      //std::cout << i << "\t" << Z[i][1] << "\t" << "Z2" << std::endl;
    }else if(Z[i][1]==3){
      *Z3+=1.0/sample;
      //std::cout << i << "\t" << Z[i][1] << "\t" << "Z3" << std::endl;
    }else if(Z[i][1]==4){
      *Z4+=1.0/sample;
      //std::cout << i << "\t" << Z[i][1] << "\t" << "Z4" << std::endl;
    }else if(Z[i][1]==5){
      *Z5+=1.0/sample;
      //std::cout << i << "\t" << Z[i][1] << "\t" << "Z5" << std::endl;
    }else if(Z[i][1]==6){
      *Z6+=1.0/sample;
      //std::cout << i << "\t" << Z[i][1] << "\t" << "Z6" << std::endl;
    }
    
  }
}

void output_Zeach(int (*Z)[2]){
  char filename[128];
  std::ofstream file;
  sprintf(filename,"Zeach_FIRE_dphi%1.1e_many.dat",dphi);
  file.open(filename,std::ios::app);
  for(int i=0;i<N1+N2;i++){
    file<<Z[i][0]<<"\t"<<Z[i][1]<<std::endl;
    //cout<<dr*(k+0.5)<<"\t"<<hist[k]<<endl;
  }
  file.close();
}

void output_Zaverage(double phi,double Zave){
  char filename[128];
  std::ofstream file;
  sprintf(filename,"Zaverage_FIRE_dphi1.0e-06_debug.dat");
  file.open(filename,std::ios::app);
  file << std::setprecision(10) <<thresh<<"\t"<<Zave<<std::endl;
  file.close();
}

void output_rattlers(double phi,double Z0){
  char filename[128];
  std::ofstream file;
  sprintf(filename,"rattlers_FIRE_dphi1.0e-06_debug.dat");
  file.open(filename,std::ios::app);
  file<<phi<<"\t"<<Z0<<std::endl;
  file.close();
}

void output_gr(double *gr,double *gr_SS,double *gr_LL,double *gr_SL){
  char filename[128];
  std::ofstream file;
  sprintf(filename,"gr_FIRE_dphi1.0e-06_debug.dat");
  file.open(filename);
  for(int k=0;k<(int)rmax/dr;k++){
    file<<dr*(k+0.5)<<"\t"<<gr[k]<<std::endl;
  }
  file.close();
  std::ofstream file_SS;
  sprintf(filename,"grSS_FIRE_dphi1.0e-06_debug.dat");
  file_SS.open(filename);
  for(int k=0;k<(int)rmax/dr;k++){
    file_SS<<dr*(k+0.5)<<"\t"<<gr_SS[k]<<std::endl;
  }
  file_SS.close();
  std::ofstream file_LL;
  sprintf(filename,"grLL_FIRE_dphi1.0e-06_debug.dat");
  file_LL.open(filename);
  for(int k=0;k<(int)rmax/dr;k++){
    file_LL<<dr*(k+0.5)<<"\t"<<gr_LL[k]<<std::endl;
  }
  file_LL.close();
  std::ofstream file_SL;
  sprintf(filename,"grSL_FIRE_dphi1.0e-06_debug.dat");
  file_SL.open(filename);
  for(int k=0;k<(int)rmax/dr;k++){
    file_SL<<dr*(k+0.5)<<"\t"<<gr_SL[k]<<std::endl;
  }
  file_SL.close();
}

int main(){
  double phiJ[sample],L[sample],x[N1+N2][3];
  double Zave,Z0,Z1,Z2,Z3,Z4,Z5,Z6,phi;
  double gr[500],gr_SS[500],gr_LL[500],gr_SL[500];
  int i,j,k,Z[N1+N2][2];
  
  // initialize
  for(k=0;k<(int)rmax/dr;k++){
    gr[k]=0.0;
    gr_SS[k]=0.0;
    gr_LL[k]=0.0;
    gr_SL[k]=0.0;
  }
  //std::cout << "initiaized" << std::endl;

  std::ifstream file_J("phiJ_FIRE_many.dat");
  std::ifstream file_c("coord_FIRE_many_dphi1.0e-06.dat");

  for(j=0;j<sample;j++){
    // input datfiles
    file_J >> phiJ[j];
    L[j]=sqrt((N1*a1/2.0*a1/2.0+N2*a2/2.0*a2/2.0)*M_PI/phiJ[j]);
    //std::cout << phiJ[j] << "\t" << L[j] << std::endl;

    for(i=0;i<N1+N2;i++){
      file_c >> x[i][0] >> x[i][1] >> x[i][2];
    } 
    histogram(x,Z,&Z0,&Z1,&Z2,&Z3,&Z4,&Z5,&Z6,&Zave,gr,gr_SS,gr_LL,gr_SL,L[j]);
    //output_Zeach(Z);
  }

  file_J.close();
  file_c.close();
  output_Zaverage(dphi,Zave);
  //output_rattlers(dphi,Z0+Z1+Z2);
  //output_gr(gr,gr_SS,gr_LL,gr_SL);
  std::cout<<"Z_average is "<<Zave<<std::endl;
  std::cout<<"Z0="<<Z0<<"\t"<<"Z1="<<Z1<<"\t"<<"Z2="<<Z2<<"\t"<<"Z3="<<Z3<<"\t"<<"Z4="<<Z4<<"\t"<<"Z5="<<Z5<<"\t"<<"Z6="<<Z6<<std::endl;

  return 0;
}