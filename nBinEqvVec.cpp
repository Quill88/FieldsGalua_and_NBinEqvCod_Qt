#include "nBinEqvVec.h"


nBinEqvVec::nBinEqvVec() : Ab(0), Ap(0)
{
	Ca = "";
	ab = nullptr;
	a = nullptr;
	CaInt = nullptr;
}

nBinEqvVec::nBinEqvVec(QString Ca) : Ca(Ca), Ab(0), Ap(0)
{
	ab = nullptr;
	a = nullptr;
	CaInt = nullptr;
}

nBinEqvVec::~nBinEqvVec()
{
	if (ab) delete[] ab;
	if (a) delete[] a;
	if (CaInt) delete[] CaInt;
}

QString nBinEqvVec::ToStr()
{
    QString s = A.get_str().c_str();
    s += "   ";
    s += A.get_str(2).c_str();
    s += "  {";
    s += Ab.get_str().c_str();
    s += "} {";

	for (int i = 0; i < n; ++i)
	{
		s += QString::number(ab[i]); 
		if (i != (n - 1)) s += ",";
	}

	s += "} {";
	s += Ap.get_str().c_str(); s += "} {";

	for (int i = 0; i < w; ++i)
	{
		s += QString::number(a[i]);
		if (i != (w - 1)) s += ",";
	}

	s += "}  ";
	s += Ca;

    return s;
}
