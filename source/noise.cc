#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <complex>
#include "wigner3j.hh"
#define PI 3.14159265358979323846
#define l_max 200

/*************************** Eq 33***************************/
double f(int l1, int l2, int L, std::vector<double> C_l1){
  return std::sqrt((2.*l1+1.)*(2.*l2+1.)/(4.*PI*(2.*L+1.)))*wigner3j(l2,l1,L,0,0,0)
          *(C_l1[l1]+pow_minus1(l1+l2+L)*C_l1[l2]);
}

/******************** weighting factor***********************/
std::complex<double> g(int l1, int l2, int L){
  return std::complex<double>(1.,0.);
}

/*************************** Eq 36***************************/
std::complex<double> A(int L, int M){
  std::complex<double> sum (0.,0.);
  std::complex<double> sum0 (0.,0.);
  int l1,l2;
  for(l1=0;l1<=l_max;l1++){
    sum0=std::complex<double>(0.,0.);
    for(l2=0;l2<=l_max;l2++){
      sum0 = sum0 + g(l1,l2,L);
    }
    sum = sum + sum0;
  }
  return std::complex<double>(1.,0.)/sum;
}

/******************** diagonal variance***********************/
std::complex<double> diag_variance(int L, int M , std::vector<double> C_l, std::vector<double> C_l1){
  std::complex<double> sum3 (0.,0.);
  std::complex<double> sum2 (0.,0.);
  std::complex<double> sum1 (0.,0.);
  std::complex<double> sum0 (0.,0.);
  int l1,l2,m1, m2;

  for(l1=0; l1<=l_max;l1++){
    sum2 = std::complex<double>(0.,0.);
    for(l2=0; l2<=l_max;l2++){
      sum1 = std::complex<double>(0.,0.);
      if((L>=std::abs(l1-l2))&&(L<=(l1+l2))&&((l1+l2+L)%2==0)){
      for(m1=-l1;m1<=l1;m1++){
        sum0 = std::complex<double>(0.,0.);
        if((L%2==0)){
          if(M==0){
            for(m2=-l2;m2<=l2;m2++){
              sum0 = sum0 +g(l1,l1,L)*std::conj(wigner3j(l1,l1,L,m1,-m1,0)
                      *wigner3j(l2,l2,L,m2,-m2,0)/(f(l1,l1,L,C_l1)*f(l2,l2,L,C_l1)))*std::conj(g(l2,l2,L));
            }
            sum0 = sum0 + g(l1,l2,L)*std::conj(wigner3j(l2,l1,L,m1,-m1,0)*
                    wigner3j(l2,l1,L,m1,-m1,0)/(f(l1,l2,L,C_l1)))*(
                      std::conj(g(l1,l2,L))*std::conj(1./f(l1,l2,L,C_l1))+std::conj(pow_minus1(l1+l2+L))*std::conj(g(l2,l1,L))*std::conj(1./f(l2,l1,L,C_l1)));
          }
          else{
            sum0 = sum0 + g(l1,l2,L)*std::conj(wigner3j(l2,l1,L,m1-M,-m1,M)*
                    wigner3j(l2,l1,L,m1-M,-m1,M)/(f(l1,l2,L,C_l1)))*(
                      std::conj(g(l1,l2,L))*std::conj(1./f(l1,l2,L,C_l1))+std::conj(pow_minus1(l1+l2+L))*std::conj(g(l2,l1,L))*std::conj(1./f(l2,l1,L,C_l1)));
          }
        }
        else{
          sum0 = sum0 + g(l1,l2,L)*std::conj(wigner3j(l2,l1,L,m1-M,-m1,M)*
                    wigner3j(l2,l1,L,m1-M,-m1,M)/(f(l1,l2,L,C_l1)))*(
                      std::conj(g(l1,l2,L))*std::conj(1./f(l1,l2,L,C_l1))+std::conj(pow_minus1(l1+l2+L))*std::conj(g(l2,l1,L))*std::conj(1./f(l2,l1,L,C_l1)));
        }

        sum1 = sum1 +sum0;
      }}
      sum2 = sum2 + sum1*std::conj(C_l[l2]);
    }
    sum3 = sum3 + sum2*std::conj(C_l[l1]);
  }
  return A(L,M)*A(L,M)*sum3;
}

int
main(int argc, char* argv[])
{
  std::ifstream file1("data/Psi_N.txt");
  std::vector<double> Psi_N;
  double temp_read;
  
  while(!file1.eof())
  {
  file1 >> temp_read;
  Psi_N.push_back(temp_read);
  }
  
  std::ifstream file2("data/monopole.txt");
  std::vector<double> monopole;

  while(!file2.eof())
  {
  file2 >> temp_read;
  monopole.push_back(temp_read);
  }

  std::ifstream file3("data/klin_cmb.txt");
  std::vector<double> klin_cmb;

  while(!file3.eof())
  {
  file3 >> temp_read;
  klin_cmb.push_back(temp_read);
  }

  std::ifstream file4("data/tt.txt");
  std::vector<double> C_CMB = {0.,0.};

  while(!file4.eof())
  {
  file4 >> temp_read;
  C_CMB.push_back(temp_read);
  }

  std::ifstream file5("data/C_l1.txt");
  std::vector<double> C_l1 ;

  while(!file5.eof())
  {
  file5 >> temp_read;
  C_l1.push_back(temp_read);
  }

  Psi_N.pop_back();
  monopole.pop_back();
  klin_cmb.pop_back();
  C_CMB.pop_back();
  C_l1.pop_back();

  int nklin_cmb = klin_cmb.size();
  int L;
  std::ofstream fout("/global/homes/p/peikail/TimeDelay/source/diag_variance/L.txt");
  for(L=0;L<=50;L++){
    std::cout << L << "\n";
    fout <<diag_variance(L,0,C_CMB,C_l1)<<"\n";
  }
  fout.close();
  return 0;
}
