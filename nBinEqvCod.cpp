#include "nBinEqvCod.h"


nBinEqvCod::nBinEqvCod(int n, int w, int q) : n(n), w(w), q(q)
{
	M = qPow(q - 1, w) * (fact(n)) / (fact(w)*fact(n - w));
	
	qDebug() << "M: " << M;
}

int nBinEqvCod::fact(const int& n)
{
	return n ? (n * fact(n - 1)) : 1;
}

QString nBinEqvCod::getEqvVec(int A)
{
	if (0 <= A && A < M)
	{
		return calc_eVec(A).ToStr();
	}
	else return "A not in 0<=A<M";
}

nBinEqvCod::~nBinEqvCod()
{
}
