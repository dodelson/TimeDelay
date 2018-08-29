#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <assert.h>
#include <stdio.h>
#include <vector>
#include <complex>
#define PI 3.14159265358979323846
#define l_min 1000
#define l_max 4000
#include "../wignerpy/src/_wignerpy/include/wignerSymbols/commonFunctions.h"
#include "../wignerpy/src/_wignerpy/include/wignerSymbols/wignerSymbols-cpp.h"

/*The symbol of these functions are the same as Okamoto & Hu's paper */
double F(int l1, int L,int l2){
	return (L*(L+1.)+l2*(l2+1.)-l1*(l1+1.))*std::sqrt((2.*l1+1.)*(2.*l2+1.)*(2.*L+1.)/(16.*PI))*WignerSymbols::wigner3j(l1,L,l2,0.,0.,0.);
}

double f(int l1,int L,int l2,std::vector<double> C_l){
	return C_l[l2]*F(l1,L,l2)+C_l[l1]*F(l2,L,l1);
}

double g(int l1, int L,int l2,std::vector<double> C_l,std::vector<double> C_l_noise){
	return f(l1,L,l2,C_l)/(2.*C_l_noise[l1]*C_l_noise[l2]);
}

double A(int L, std::vector<double> C_l,std::vector<double> C_l_noise){
	double sum0=0.;
	double sum1=0.;
	int l1, l2;
	for(l1=l_min;l1<=l_max;l1++){
		sum0=sum0+sum1;
		for(l2=l_min;l2<=l_max;l2++){
			sum1 = sum1 + g(l1,L,l2,C_l,C_l_noise)*f(l1,L,l2,C_l);
		}
		std::cout << l1 << " " << L*(L+1.)*(2.*L+1.)/sum0 << "\n";
	}
	return L*(L+1.)*(2.*L+1.)/sum0;
}

double diag_variance(int L , std::vector<double> C_l,std::vector<double> C_l_noise){
	return A(L,C_l,C_l_noise);
}

int
main(int argc, char* argv[])
{
  double temp_read;
  
  /*TT spectrum I restored is C_{l}/(2*PI), so I timed a 2*PI here*/
  std::ifstream file0("data/tt.txt");
  std::vector<double> C_TT={0.,0.};

  while(!file0.eof())
  {
  file0 >> temp_read;
  C_TT.push_back(temp_read*2.*PI);
  }

  std::ifstream file1("data/C_TT_1.txt");
  std::vector<double> C_TT1={0.,0.};

  while(!file1.eof())
  {
  file1 >> temp_read;
  C_TT1.push_back(temp_read*2.*PI);
  }

  C_TT.pop_back();
  C_TT1.pop_back();

  /* NOISE term for TT power spectrum, I ploted this and checked it with "MASS RECONSTRUCTION WITH CMB POLARIZATION", Hu & Okamoto
    ,so this term should be correct...*/

  std::vector<double> C_TT_noise;

  double coeff_noise = (PI/10800.)*(PI/10800.);
  double fwhm2_noise = 16.*(PI/10800.)*(PI/10800.);

  for(int i=0;i<=5000;i++){
    C_TT_noise.push_back(C_TT[i]+coeff*std::exp(i*(i+1.)*fwhm2_noise/(8.*std::log(2))));
  }

  int L=2;
  std::cout << L << " " << L*(L+1.)*diag_variance(L,C_TT,C_TT_noise)/(2.*PI) << "\n";

  return 0;
}
