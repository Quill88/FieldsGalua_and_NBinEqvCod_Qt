#include "nBinEqvVec.h"


nBinEqvVec::nBinEqvVec()
{
	Ca = "";
	ab = nullptr;
	a = nullptr;
}


nBinEqvVec::~nBinEqvVec()
{
	if (ab) delete[] ab;
	if (a) delete[] a;
}

QString nBinEqvVec::ToStr()
{
	QString s = QString::number(A) + "\t";
	
	s += QString::number(A, 2) + "\t{";

	for (int i = 0; i < n; ++i)
	{
		s += QString::number(ab[i]); 
		if (i != (n - 1)) s += ", ";
	}
	s += "}\t{";

	for (int i = 0; i < w; ++i)
	{
		s += QString::number(a[i]);
		if (i != (w - 1)) s += ", ";
	}

	s += "}\t";
	s += Ca;

    return s;
}
