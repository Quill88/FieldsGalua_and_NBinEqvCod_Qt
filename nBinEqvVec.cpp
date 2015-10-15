#include "nBinEqvVec.h"


nBinEqvVec::nBinEqvVec()
{
	Ca = "";
}


nBinEqvVec::~nBinEqvVec()
{
}

QString nBinEqvVec::ToStr()
{
    return QString::number(A)+"\t"+Ca;
}
