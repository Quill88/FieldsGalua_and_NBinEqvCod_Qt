#include "nBinEqvCod.h"


nBinEqvCod::nBinEqvCod(int n, int w, int q) : n(n), w(w), q(q)
{
	qw = qPow((q - 1), w);
	M = qw * (fact(n)) / (fact(w)*fact(n - w));
	
	qDebug() << "Creating non_binary equivalent codes";
	qDebug() << "n:" << n <<"\tw: " << w <<"\tq:"<<q <<"\tM: " <<M;
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
        return calc_eVec(A).ToStr();
	}
	else return "A not in 0<=A<M";
}

nBinEqvVec nBinEqvCod::calc_eVec(int A)
{
    nBinEqvVec v;
	v.A = A;
		
	int Ab = A / qw;
	int Ap = A % qw;

	int* ab = new int[n];
	int* a = new int[w];


	/*иру 3*/
	int x = Ab;
	int l = 0;

	for (int i = 0; i < n; ++i)
	{
		int b = fact(n - i - 1) / (fact(w - l)*fact((n - i - 1) - (w - l)));
		if (b > x) ab[n - i - 1] = 0;
		else
		{
			ab[n - i - 1] = 1;
			x -= b;
			l++;
		}
	}

	/*иру 4*/
	x = Ap;
	l = 0;
	while (l < w)
	{
		a[l] = (x % (q - 1)) + 1;
		x = x / (q - 1);
		l++;
	}

	/*иру 5*/
	for (int i = 0, l = 0; i < n; ++i)
	{
		if (ab[i])
		{
			v.Ca += QString::number(a[l]);
			l++;
		}
		else v.Ca += "0";
	}
	
	delete ab;
	delete a;

    return v;
}

nBinEqvCod::~nBinEqvCod()
{
}
