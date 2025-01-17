#ifndef __COMPLEX_H__
#define __COMPLEX_H__
#pragma once
#include "exportimport.h"

#include "afxtempl.h"
#include "math.h"

#include <string>
#include <cmath>

#include "ComplexArray.h"

class CLASS_CMPLXLIB_DLL complex
{
public:
	//complex(void);
	//complex(double xInit=0, double yInit=0);
	complex(double xInit = 0, double yInit = 0)
		: x(xInit), y(yInit), FlagError(false) {}
	~complex(void) {};
	double x,y;
	std::string TextError;
	bool FlagError;

};

class CLASS_CMPLXLIB_DLL Pixel : public complex
{
public:
	//Pixel(void){};
	Pixel(double xInit=0, double yInit=0){x=xInit; y=yInit;Break=false;};
	~Pixel(void){};
	bool Break;
};

class CLASS_CMPLXLIB_DLL CFraction
{
public:
	int num,dem, part;
	CFraction();
	CFraction(int x, int y){num=x, dem=y;};
	CFraction(double value, double EPS);
};

extern "C" CLASS_CMPLXLIB_DLL complex operator+(complex A, complex &B);
extern "C" CLASS_CMPLXLIB_DLL complex operator-(complex A, complex &B);
extern "C" CLASS_CMPLXLIB_DLL complex operator+=(complex A, complex B);
extern "C" CLASS_CMPLXLIB_DLL complex operator*(complex &A, complex &B);
extern "C" CLASS_CMPLXLIB_DLL complex operator*(double A, complex B);
extern "C" CLASS_CMPLXLIB_DLL complex operator*(complex A, double B);

complex CLASS_CMPLXLIB_DLL operator/(complex &A,complex &B);
complex CLASS_CMPLXLIB_DLL operator/(complex A,double B);
complex CLASS_CMPLXLIB_DLL operator/(double A,complex &B);
complex CLASS_CMPLXLIB_DLL operator+(double A,complex B);
complex CLASS_CMPLXLIB_DLL operator+(complex A,double B);
complex CLASS_CMPLXLIB_DLL operator-(double A,complex B);
complex CLASS_CMPLXLIB_DLL operator-(complex A,double B);
complex CLASS_CMPLXLIB_DLL operator+(complex &A);
complex CLASS_CMPLXLIB_DLL operator-(complex A);
complex CLASS_CMPLXLIB_DLL operator+=(complex A, complex B);
bool CLASS_CMPLXLIB_DLL operator!=(complex &A, complex &B);
int CLASS_CMPLXLIB_DLL operator==(complex &A,complex &B);
complex CLASS_CMPLXLIB_DLL conj(complex A);
complex CLASS_CMPLXLIB_DLL exp(complex &A);
complex CLASS_CMPLXLIB_DLL log(complex &A);
complex CLASS_CMPLXLIB_DLL pow(complex &A,double B);
double CLASS_CMPLXLIB_DLL arg(complex &A, bool bDegree=TRUE);
double CLASS_CMPLXLIB_DLL abs(complex &A);
double CLASS_CMPLXLIB_DLL abs_db(complex &A);
complex CLASS_CMPLXLIB_DLL polar(double mag, double angle, bool bDegree=TRUE);
complex CLASS_CMPLXLIB_DLL sqrt(complex &A);
double CLASS_CMPLXLIB_DLL real(complex &A);
double CLASS_CMPLXLIB_DLL imag(complex &A);
//double abs(double A);
double CLASS_CMPLXLIB_DLL norm(complex &A);
complex CLASS_CMPLXLIB_DLL sin(complex A);
complex CLASS_CMPLXLIB_DLL cos(complex A);
complex CLASS_CMPLXLIB_DLL tan(complex A);
complex CLASS_CMPLXLIB_DLL lg(complex A);
complex CLASS_CMPLXLIB_DLL acos(complex A);
complex CLASS_CMPLXLIB_DLL acosh(complex A);
complex CLASS_CMPLXLIB_DLL asin(complex A);
complex CLASS_CMPLXLIB_DLL asinh(complex A);
complex CLASS_CMPLXLIB_DLL atan(complex A);
complex CLASS_CMPLXLIB_DLL atanh(complex A);
complex CLASS_CMPLXLIB_DLL acot(complex A);
complex CLASS_CMPLXLIB_DLL ceil(complex A, bool& f);
complex CLASS_CMPLXLIB_DLL floor(complex A, bool& f);
complex CLASS_CMPLXLIB_DLL cosh(complex A);
complex CLASS_CMPLXLIB_DLL coth(complex A);
complex CLASS_CMPLXLIB_DLL sinh(complex A);
complex CLASS_CMPLXLIB_DLL csc(complex A);
complex CLASS_CMPLXLIB_DLL csch(complex A);
complex CLASS_CMPLXLIB_DLL Im(complex A);
complex CLASS_CMPLXLIB_DLL Re(complex A);
complex CLASS_CMPLXLIB_DLL sec(complex A);	
complex CLASS_CMPLXLIB_DLL asec(complex A);
complex CLASS_CMPLXLIB_DLL acsc(complex A);
complex CLASS_CMPLXLIB_DLL sech(complex A);
complex CLASS_CMPLXLIB_DLL tanh(complex A);
complex CLASS_CMPLXLIB_DLL asech(complex A);
complex CLASS_CMPLXLIB_DLL acsch(complex A);
complex CLASS_CMPLXLIB_DLL modul(complex A);
complex CLASS_CMPLXLIB_DLL pow(complex A, complex B);
complex CLASS_CMPLXLIB_DLL pow(complex A, complex B, bool& P);
double CLASS_CMPLXLIB_DLL Atan2(double y, double x);
complex CLASS_CMPLXLIB_DLL SQRT(complex B, complex A);
complex CLASS_CMPLXLIB_DLL SQRT1(complex A, complex B);
complex CLASS_CMPLXLIB_DLL log(complex A, complex B);
complex CLASS_CMPLXLIB_DLL signum(complex A);
double CLASS_CMPLXLIB_DLL modul1(complex A);
complex	CLASS_CMPLXLIB_DLL ln(complex A);
complex CLASS_CMPLXLIB_DLL sqrt2(double a, double n);
complex CLASS_CMPLXLIB_DLL acoth(complex A);


// Определение оператора сложения
//extern "C" CLASS_CMPLXLIB_DLL complex operator+(const complex A, const complex &B) {
//	return complex(A.x + B.x, A.y + B.y);
//}


//complex CLASS_CMPLXLIB_DLL operator+(complex A, const complex &B) {
//	A.x += B.x;
//	A.y += B.y;
//	return A;
//}



//typedef CArray <complex, complex> ComplexArray;
//	int _matherr(struct _exception* pExp);
//	complex CLASS_CMPLXLIB_DLL operator+(complex A,complex &B);
//	complex CLASS_CMPLXLIB_DLL operator-(complex A,complex &B);
//	complex CLASS_CMPLXLIB_DLL operator*(complex &A,complex &B);
//	complex CLASS_CMPLXLIB_DLL operator*(double A,complex B);
//	complex CLASS_CMPLXLIB_DLL operator*(complex A,double B);
//	complex CLASS_CMPLXLIB_DLL operator/(complex &A,complex &B);
//	complex CLASS_CMPLXLIB_DLL operator/(complex A,double B);
//	complex CLASS_CMPLXLIB_DLL operator/(double A,complex &B);
//	complex CLASS_CMPLXLIB_DLL operator+(double A,complex B);
//	complex CLASS_CMPLXLIB_DLL operator+(complex A,double B);
//	complex CLASS_CMPLXLIB_DLL operator-(double A,complex B);
//	complex CLASS_CMPLXLIB_DLL operator-(complex A,double B);
//	complex CLASS_CMPLXLIB_DLL operator+(complex &A);
//	complex CLASS_CMPLXLIB_DLL operator-(complex A);
//	complex CLASS_CMPLXLIB_DLL operator+=(complex A, complex B);
//	bool CLASS_CMPLXLIB_DLL operator!=(complex &A, complex &B);
// 	int CLASS_CMPLXLIB_DLL operator==(complex &A,complex &B);
//	complex CLASS_CMPLXLIB_DLL conj(complex A);
//	complex CLASS_CMPLXLIB_DLL exp(complex &A);
//	complex CLASS_CMPLXLIB_DLL log(complex &A);
//	complex CLASS_CMPLXLIB_DLL pow(complex &A,double B);
//	double CLASS_CMPLXLIB_DLL arg(complex &A, bool bDegree=TRUE);
//	double CLASS_CMPLXLIB_DLL abs(complex &A);
//	double CLASS_CMPLXLIB_DLL abs_db(complex &A);
//	complex CLASS_CMPLXLIB_DLL polar(double mag, double angle, bool bDegree=TRUE);
//	complex CLASS_CMPLXLIB_DLL sqrt(complex &A);
//	double CLASS_CMPLXLIB_DLL real(complex &A);
//	double CLASS_CMPLXLIB_DLL imag(complex &A);
//	//double abs(double A);
//	double CLASS_CMPLXLIB_DLL norm(complex &A);
//	complex CLASS_CMPLXLIB_DLL sin(complex A);
//	complex CLASS_CMPLXLIB_DLL cos(complex A);
//	complex CLASS_CMPLXLIB_DLL tan(complex A);
//	complex CLASS_CMPLXLIB_DLL lg(complex A);
//	complex CLASS_CMPLXLIB_DLL acos(complex A);
//	complex CLASS_CMPLXLIB_DLL acosh(complex A);
//	complex CLASS_CMPLXLIB_DLL asin(complex A);
//	complex CLASS_CMPLXLIB_DLL asinh(complex A);
//	complex CLASS_CMPLXLIB_DLL atan(complex A);
//	complex CLASS_CMPLXLIB_DLL atanh(complex A);
//	complex CLASS_CMPLXLIB_DLL acot(complex A);
//	complex CLASS_CMPLXLIB_DLL ceil(complex A, bool& f);
//	complex CLASS_CMPLXLIB_DLL floor(complex A, bool& f);
//	complex CLASS_CMPLXLIB_DLL cosh(complex A);
//	complex CLASS_CMPLXLIB_DLL coth(complex A);
//	complex CLASS_CMPLXLIB_DLL sinh(complex A);
//	complex CLASS_CMPLXLIB_DLL csc(complex A);
//	complex CLASS_CMPLXLIB_DLL csch(complex A);
//	complex CLASS_CMPLXLIB_DLL Im(complex A);
//	complex CLASS_CMPLXLIB_DLL Re(complex A);
//    complex CLASS_CMPLXLIB_DLL sec(complex A);	
//	complex CLASS_CMPLXLIB_DLL asec(complex A);
//	complex CLASS_CMPLXLIB_DLL acsc(complex A);
//	complex CLASS_CMPLXLIB_DLL sech(complex A);
//	complex CLASS_CMPLXLIB_DLL tanh(complex A);
//	complex CLASS_CMPLXLIB_DLL asech(complex A);
//	complex CLASS_CMPLXLIB_DLL acsch(complex A);
//	complex CLASS_CMPLXLIB_DLL modul(complex A);
//	complex CLASS_CMPLXLIB_DLL pow(complex A, complex B);
//	complex CLASS_CMPLXLIB_DLL pow(complex A, complex B, bool& P);
//	double CLASS_CMPLXLIB_DLL Atan2(double y, double x);
//	complex CLASS_CMPLXLIB_DLL SQRT(complex B, complex A);
//	complex CLASS_CMPLXLIB_DLL SQRT1(complex A, complex B);
//	complex CLASS_CMPLXLIB_DLL log(complex A, complex B);
//	complex CLASS_CMPLXLIB_DLL signum(complex A);
//	double CLASS_CMPLXLIB_DLL modul1(complex A);
//	complex	CLASS_CMPLXLIB_DLL ln(complex A);
//	complex CLASS_CMPLXLIB_DLL sqrt2(double a, double n);
//	complex CLASS_CMPLXLIB_DLL acoth(complex A);
//	complex Fraction(double a, double EPS);

#endif