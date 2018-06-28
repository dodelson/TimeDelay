#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "wigner3j.hh"

double pow_minus1(int n){
	if((n%2)==0){
		return 1.;
	}
	else{
		return -1.;
	}
}
long double factorial(int n)
{	
	long double temp=1.;
    if(n > 1)
        {for(;n>1;n--)
        {
        	temp=temp*double(n);
        }
        return temp;}
    else
        return 1.;
}

double triangle_coefficient(int j1, int j2, int j3){
	return factorial(j1+j2-j3)*factorial(j1-j2+j3)*factorial(-j1+j2+j3)
	 		/factorial(j1+j2+j3+1);
}

long double sum_coefficient(int t,int j1, int j2, int j3, int m1 ,int m2){
	return factorial(t)*factorial(j3-j2+t+m1)*factorial(j3-j1+t-m2)*factorial(j1+j2-j3-t)
			*factorial(j1-t-m1)*factorial(j2-t+m2);
}

long double sum_t(int j1, int j2, int j3, int m1, int m2){
	int min1 = 0;
	int min2 = j2 - j3 - m1;
	int min3 = j1 - j3 + m2;
	int min12 = std::max(min1, min2);
	int min =  std::max(min12, min3);
	int max1 = j1 + j2 - j3;
	int max2 = j1 - m1;
	int max3 = j2 + m2;
	int max12 = std::min(max1, max2);
	int max = std::min(max12, max3);

	long double temp = 0.;
	int t;
	for(t=min;t<=max;t++){
		temp = temp + pow_minus1(t)/sum_coefficient(t, j1, j2, j3, m1, m2);
	}
	return temp;
}

double wigner3j(int j1, int j2, int j3, int m1, int m2, int m3){
    if((m1>=-j1)&&(m1<=j1)&&(m2>=-j2)&&(m2<=j2)&&(m3>=-j3)&&(m3<=j3)
        &&(m1+m2+m3==0)&&(j3>=std::abs(j1-j2))&&(j3<=(j1+j2)))
    {
        if((m1==0)&&(m2==0)&&(m3==0)&&((j1+j2+j3)%2==1)){
        	return 0.;
        }
        else{
        	return pow_minus1(j1-j2-m3)*std::sqrt(triangle_coefficient(j1, j2, j3)*
        			factorial(j1+m1)*factorial(j1-m1)*factorial(j2+m2)*
        			factorial(j2-m2)*factorial(j3+m3)*factorial(j3-m3))*
        			sum_t(j1,j2,j3,m1,m2);
        }
    }
    else{
        return 0.;
    }
}