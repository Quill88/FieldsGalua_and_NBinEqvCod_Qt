#pragma once
#ifndef NBINEQVVEC_H
#define NBINEQVVEC_H

#include <QString>

/*nonЦbinary equivalent vector
недвоичный равновесный вектор*/
class nBinEqvVec
{
public:
	nBinEqvVec();
	~nBinEqvVec();

    QString ToStr();
	
	int A;
	int n;
	int w;

	int* ab;
	int* a;

    QString Ca;
};

#endif /* NBINEQVVEC_H */
