#include "StdAfx.h"
#include "complex.h"
#include "MathConstants.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

#define EPS_EQ 10e-15

complex::complex(double xInit, double yInit)
{
	x=xInit; y=yInit;
	FlagError=false;
}

complex::~complex(void)
{
}

complex operator+(complex A,complex &B)
{
 A.x+=B.x;
 A.y+=B.y;
return A;
};

complex operator-(complex A,complex &B)
{
 A.x-=B.x;
 A.y-=B.y;
return A;
};

complex operator+=(complex A, complex B)
{
	A.x+=B.x;
	A.y+=B.y;
	return A;
}

complex operator*(complex &A,complex &B)
{
 complex C(0,0);
 C.x=A.x*B.x-A.y*B.y;
 C.y=A.y*B.x+B.y*A.x;
return C;
};

complex operator*(double A,complex B)
{
 B.x*=A;
 B.y*=A;
return B;
};

complex operator*(complex A,double B)
{
 A.x*=B;
 A.y*=B;
return A;
};

complex operator/(double A,complex &B)
{
	complex C;
	if(B==complex())
	{
		C.TextError="ƒеление на ноль";
		C.FlagError=true;
		return C;
	}
	C=A * conj(B) / norm(B);
	return C;
};

complex operator/(complex &A,complex &B)
{complex C;
  C=A * conj(B) / norm(B);
return C;
};

complex operator/(complex A,double B)
{
 A.x/=B;
 A.y/=B;
return A;
};

complex operator+(double A,complex B)
{
 B.x+=A;
return B;
};


complex operator+(complex A,double B)
{
 A.x+=B;
return A;
};

complex operator-(double A,complex B)
{
 B.x=A-B.x;
 B.y=-B.y;
return B;
};

complex operator-(complex A,double B)
{
 A.x-=B;
return A;
};

complex operator+(complex &A)
{
return A;
}

complex operator-(complex A)
{
A.x=-A.x;
A.y=-A.y;
return A;
}

bool operator!=(complex &A, complex &B)
{
	return !(A==B);
}

//ѕоиск сопр€женного
complex conj(complex A)
{
A.y=-A.y;
return A;
}

complex exp(complex &A)
{
	complex C;
	//ASSERT((A.x<7.097827e+002)&&(A.x>-7.083964e+002));
	//ASSERT((A.y<7.097827e+002)&&(A.y>-7.083964e+002));

	double a=exp(A.x);
	C.x=a*cos(A.y);
	C.y=a*sin(A.y);
	return C;
}

complex log(complex &A)
{
	complex C;
	complex M = modul(A);
	complex S(0,0);
	if(M==S)
		return S;
	C.x=log(abs(A));
	C.y=arg(A,false);
return C;
}

complex pow(complex &A,double B)
{
return exp(B*log(A));
}

double arg(complex &A, bool bDegree)
{
	if (bDegree)
		return Atan2(A.x,A.y)*180/M_PI;
	else
		return Atan2(A.x,A.y);
}

double abs(complex &A)
{
return sqrt(A.x * A.x + A.y * A.y);
}


double abs_db(complex &A)
{
return 10*log10(sqrt(A.x * A.x + A.y * A.y));
}


complex polar(double mag, double angle, bool bDegree)
{
	if (bDegree) angle=angle*M_PI/180;
	complex c;
	c.x=mag * cos(angle);
	c.y=mag * sin(angle);
	return c;
}

complex sqrt(complex &A)
{
return polar(sqrt(abs(A)), arg(A) / 2);
}

double real(complex &A)
{
return A.x;
}

double imag(complex &A)
{
	return A.y;
}

/*double abs(double A)
{
	if (A>0) return A;
	return -A;
}*/

int operator==(complex &A,complex &B)
{
	return (fabs(A.x-B.x)<EPS_EQ) && (fabs(A.y-B.y))<EPS_EQ;
}

double norm(complex &A)
{
return A.x*A.x+A.y*A.y;
}

complex sin(complex A)
{
	complex b;
	complex i;//мнима€ еденица
	i.x=0;i.y=1;
	/*if(A.y==0)
		b.x=sin(A.x);
	else
	{*/
		b=(exp(A*i)-exp(-A*i))/(2*i);
	//}
	return b;
}

complex cos(complex A)
{
	complex b;
	complex i;
	i.x=0;
	i.y=1;	
	if(A.y==0)
		b.x=cos(A.x);
	else
		b=(exp(A*i)+exp(-A*i))/(2);
	//return cos(A.x);
	return b;
}

complex tan(complex A)
{
	complex B;
	if(abs(cos(A))!=0)
	{
		complex C1=sin(A);
		complex C2=cos(A);
		//B=sin(A)/cos(A);
		B=C1/C2;
	}
	return B;
}

complex lg(complex A)
{
	complex B;
	//if(A.y==0)
	//	B.x=log10(A.x);
	//else
	//{
		B.x=log10(abs(A));
		B.y=arg(A,false)/log(10.);			
	//}
	return B;
}

complex acos(complex A)
{
	complex B;
	complex S = pow(A*A-1,0.5);
	complex L=log(A+S);
	B=-i1*L;
	return B;
}

complex acosh(complex A)
{
	return log(A+sqrt(A*A-1));
}

complex asin(complex A)
{
	complex i(0,1);
	complex sin=-i1*log(i*A+pow(1-A*A,0.5));
	/*if((A.y==0)&&(A.x<0))
		sin=M_PI_2-sin;*/
	return sin;	
}

complex asinh(complex A)
{
	return log(A+sqrt(A*A+1));
}

complex atan(complex A)
{
	complex i(0,1);
	if(A==i)
	{
		A.TextError=ER_SPEZIAL_POINT;
		A.FlagError = true;
		return A;
	}
	complex b=(i-A)/(i+A);
	complex l=log(b);
	complex a=(-i/2*l);
	return a;
}

complex ln(complex A)
{
	complex i(0,1);
	double d=((complex)modul(A)).x;
	double b=arg(A,false);
	complex r=log(d)+i*b;
	return r;
}

complex atanh(complex A)
{
	if(A==complex(1,0))
	{
		A.TextError=ER_SPEZIAL_POINT;
		A.FlagError = true;
		return A;
	}
	return (log((1+A)/(1-A))/2);
}

complex ceil(complex A, bool& f)
{
	complex B(0,0);
	if(A.y!=0)
		f=false;
	else
		B.x=ceil(A.x);
	return B;
}

complex floor(complex A, bool& f)
{
	complex B(0,0);
	if(A.y!=0)
		f=false;
	else
		B.x=floor(A.x);
	return B;
}

complex cosh(complex A)
{
	return ((exp(A)+exp(-A))/2);
}

complex coth(complex A)
{
	if(sinh(A)==complex())
	{
		A.TextError=ER_SPEZIAL_POINT;
		A.FlagError = true;
		return A;
	}
	return (cosh(A)/(sinh(A)));
}

complex sinh(complex A)
{
	return ((exp(A)-exp(-A))/2);
}

complex csc(complex A)
{
	if(sin(A)==complex())
	{
		A.TextError=ER_SPEZIAL_POINT;
		A.FlagError = true;
		return A;
	}
	return (1/sin(A));
}

complex csch(complex A)
{
	if(sinh(A)==complex())
	{
		A.TextError=ER_SPEZIAL_POINT;
		A.FlagError = true;
		return A;
	}
	return (1/sinh(A));
}

complex Im(complex A)
{
	return complex(A.y,0);
}

complex Re(complex A)
{
	return complex(A.x,0);
}

complex sec(complex A)
{
	return (1/cos(A));
}

complex sech(complex A)
{
	complex cosh1=cosh(A);
	if(cosh1==complex())
	{
		A.TextError=ER_SPEZIAL_POINT;
		A.FlagError = true;
		return A;
	}
	return(1/cosh1);
}

complex tanh(complex A)
{
	complex cosh1=cosh(A);
	if(cosh1==complex())
	{
		A.TextError=ER_SPEZIAL_POINT;
		A.FlagError = true;
		return A;
	}
	return (sinh(A)/cosh1);
}

complex modul(complex A)
{
	return complex(sqrt(A.x*A.x+A.y*A.y),0);
}

double modul1(complex a)
{
	return (double)sqrt(a.x*a.x+a.y*a.y);
}

complex pow(complex A, complex B)
{
	complex x(0,0);
	if(fabs(A.y)<10e-8)
		A.y=0;
	if((A.x==0)&&(A.y==0))
	     return x;
	double p=Atan2(A.x,A.y);
	double m=modul(A).x;
	double r=log(m);
	double v=B.x*p+B.y*r;
	//v=v/B.x;
	double w=0;
	w=(double)exp(B.x*r-B.y*p);
	 x.x=w*cos(v);
	x.y=w*sin(v);
	return x;
}

complex pow(complex A, complex B, bool& P)
{
	complex x(1,0);
	/*if(fabs(A.y)<EPS_EQ)
	{
		return (B==0)?complex(1,0):complex();
	}*/
	if((A.x==0)&&(A.y==0))
	{
		if((B.y==0)&&(B.x<0))
		{
			P=true;
			return complex();
		}
	    return complex(0,0);
	}
	double p=Atan2(A.x,A.y);
	double m=modul(A).x;
	double r=log(m);
	double v=B.x*p+B.y*r;
	//v=v/B.x;
	double w=0;
	try
	{
		double f=B.x*r-B.y*p;
		if(f>709.782712893)
			throw "ѕереполнение при возведении в степень";
		w=(double)exp(f);
	}
	catch(char* str)
	{
		x.TextError = str;
		P=true;
		return x;
	}
	P=false;
	x.x=w*cos(v);
	x.y=w*sin(v);
	return x;
}

int _matherr(struct _exception* pExp)
{
	AfxMessageBox("произошла ошибка");
	TRACE(pExp->name);
	return 0;
}

double Atan2(double x, double y)
{
	double Arctan=0;
	if(y==0)
		if(x>0)
			return 0;
		else
			return M_PI;
	if(x==0) 	
	{
		if(y==0) return 0;
		else
		{
			if(y>0) return (M_PI/2);
			else return (-M_PI/2);
		}
	}
	else
	{
		if((x>0)&&(y>0)) Arctan=atan(fabs(y)/fabs(x));
		if((x<0)&&(y>0)) Arctan=M_PI-atan(fabs(y)/fabs(x));//втора€ четверть
		if((x<0)&&(y<0)) Arctan=-M_PI+atan(fabs(y)/fabs(x));	//треть€ четверть	
		if((x>0)&&(y<0)) Arctan=-atan(fabs(y)/fabs(x));//четверта€ четверть
		
	}
	//while(Arctan<-M_PI/2) Arctan+=M_PI;
	//while(Arctan>M_PI/2) Arctan-=M_PI;
	return Arctan;
}


#define SQRT_EPS 10e-6

complex SQRT(complex B,complex A)
{
	//return pow(A,1/B);
	ComplexArray* Mass = new ComplexArray;
	if(A.y!=0) //if B=complex = then 
	{
		Mass->Add(pow(A,1/B));
	}
	double n=B.x;
	double M=modul(A).x;
	double Ro=pow(M,1/n);
	double fi=arg(A,false);//gradus or radian
	for(int i=0;i<n;i++)
	{
		double psi=(fi+2*M_PI*i)/n;
		Mass->Add(complex(Ro*cos(psi),Ro*sin(psi)));
	}
	for(int i=0;i<Mass->GetCount();i++)
	{
		complex a=Mass->GetAt(i);
		if(fabs(a.y)<SQRT_EPS)
		{
			delete Mass;       
			return a;
		}
		else
			if(fabs(a.x)<SQRT_EPS)
			{
				delete Mass;
				return a;
			}
	}
	complex c = Mass->GetAt(0);
	delete Mass;
	return c;
}

complex SQRT1(complex A, complex B)
{
	return pow(A,1/B);
}

complex log(complex A, complex B)
{
	complex a=log(A);
	complex b=log(B);
	return (a/b);
}


complex signum(complex A)
{
	return (A/modul(A));
}

CFraction::CFraction(double value, double EPS)
{
	//part = value
	num=0; dem=1;
	double c;
	part = (int)value;
    value-=part;
	ASSERT(value>0);
	if((value<0)&&(part!=0))
		value*=-1;
	if(value)
	{
		do
		{
			num++;
			double s=(double)num/dem;
			c=s-value;
			if(num==dem)
			{
				dem++;
				if(part==0)
					num=-dem;
				else
					num=0;
			}
		}while(fabs(c)>EPS);
	}
	//if(num>dem
}

complex sqrt2(double a, double n)
{
	complex x(0,0);
	double p=0;
	double m=fabs(a);
	double r=log(m);
	double v=n*p;
	double w=0;
	//double p=Atan2(A.x,A.y);
	w=(double)exp(n*r);
	if(a>0)
	{
		x.x=w;
		x.y=0;
	}
	else
	{
		int b=(int)fabs(1/n);
		int f = b % 2;
		if(f==0)
		{
			x.x=0;
			x.y=w;
		}
		else
		{
			x.x=-w;
			x.y=0;
		}
	}
	return x;
}

complex acot(complex A)
{
	if(A==complex(0,1))
	{
		A.TextError=ER_SPEZIAL_POINT;
		A.FlagError = true;
		return A;
	}
	complex B=((i1/2)*ln((A-i1)/(A+i1)));
	if((fabs(B.y)<10e-10)&&(B.x<0))
		B.x+=M_PI;
	return B;
}


complex asec(complex A)
{
	//return atan(A/pow(pow(A,2)-1,2));
	if(A==complex())
	{
		A.TextError=ER_SPEZIAL_POINT;
		A.FlagError = true;
		return A;
	}
	return -i1*ln(1/A*(1+pow((1-pow(A,2)),0.5)));
}

complex acsc(complex A)	
{
	//return -i1*ln((i1+pow((A-1)*(A+1),0.5))/A);
	return asin(pow(A,-1));
}

complex acoth(complex A)
{
	if(A==complex(1,0))
	{
		A.TextError=ER_SPEZIAL_POINT;
		A.FlagError = true;
		return A;
	}
	return ((A.x>0)?ln(pow(pow(A,2)-1,0.5)/(A-1)):ln(-pow(pow(A,2)-1,0.5)/(A-1)));
}

complex asech(complex A)
{
	return ln((1+pow(1-pow(A,2),0.5))/A);
}
complex acsch(complex A)
{
	return (A.x>0)?ln((1+pow(1+pow(A,2),0.5))/A)
					:ln((1-pow(1+pow(A,2),0.5))/A);
}