#include "nBinEqvCod.h"


nBinEqvCod::nBinEqvCod(int n, int w, int q) : n(n), w(w), q(q)
{
	qw = qPow((q - 1), w);
	M = qw * (fact(n)) / (fact(w)*fact(n - w));

	code = new nBinEqvVec[M];
	
	qDebug() << "Creating non_binary equivalent codes";
	qDebug() << "n:" << n <<"\tw: " << w <<"\tq:"<<q <<"\tM: " <<M;

	for (int i = 0; i < M; ++i)
	{
		calc_eVec(i, code[i]);
	}
}

int nBinEqvCod::fact(const int& n)
{
	if (n >= 0)
		return n ? (n * fact(n - 1)) : 1;
	else
		return n;
}		

int nBinEqvCod::getM() const
{
	return M;
}

QString nBinEqvCod::getEqvVec(int A)
{
	if (0 <= A && A < M)
	{
		return code[A].ToStr();
	}
	else return "A not in 0<=A<M";
}

void nBinEqvCod::calc_eVec(int A, nBinEqvVec& v)
{
	v.A = A;
	v.n = n;
	v.w = w;

	int Ab = A / qw;
	int Ap = A % qw;

	v.ab = new int[n];
	v.a = new int[w];


	/*иру 3*/
	int x = Ab;
	int l = 0;

	for (int i = 0; i < n; ++i)
	{
		int b = fact(n - i - 1) / (fact(w - l)*fact((n - i - 1) - (w - l)));
		if (b > x) v.ab[n - i - 1] = 0;
		else
		{
			v.ab[n - i - 1] = 1;
			x -= b;
			l++;
		}
	}

	/*иру 4*/
	x = Ap;
	l = 0;
	while (l < w)
	{
		v.a[l] = (x % (q - 1)) + 1;
		x = x / (q - 1);
		l++;
	}

	/*иру 5*/
	for (int i = 0, l = 0; i < n; ++i)
	{
		if (v.ab[i])
		{
			v.Ca += QString::number(v.a[l]);
			l++;
		}
		else v.Ca += "0";
	}
}


nBinEqvCod::~nBinEqvCod()
{
	if (code) delete[] code;
}
