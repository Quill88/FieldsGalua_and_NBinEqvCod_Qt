#include "nBinEqvVec.h"


nBinEqvVec::nBinEqvVec()
{
	Ca = "";
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
	QString s = QString::number(A) + "   ";

	if (A >= 64 && A <= 127) s += "0";
	else if (A >= 32 && A <= 63) s += "00";
	else if (A >= 16 && A <= 31) s += "000";
	else if (A >= 8 && A <= 15) s += "0000";
	else if (A >= 4 && A <= 7) s += "00000";
	else if (A >= 2 && A <= 3) s += "000000";
	else if (A >= 0 && A <= 1) s += "0000000";

	s += QString::number(A, 2) + "  {";

	for (int i = 0; i < n; ++i)
	{
		s += QString::number(ab[i]); 
		if (i != (n - 1)) s += ", ";
	}
	s += "} {";

	for (int i = 0; i < w; ++i)
	{
		s += QString::number(a[i]);
		if (i != (w - 1)) s += ", ";
	}

	s += "}  ";
	s += Ca;

    return s;
}
