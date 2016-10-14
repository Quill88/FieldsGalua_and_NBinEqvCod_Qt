#pragma once
#ifndef NBINEQVVEC_H
#define NBINEQVVEC_H

#include <gmpxx.h>
#include <QString>

/*nonЦbinary equivalent vector
недвоичный равновесный вектор*/
class nBinEqvVec
{
public:
	nBinEqvVec();
	nBinEqvVec(QString);
	~nBinEqvVec();

    QString ToStr();
	
	mpz_class A;
	int n;
	int w;

    QString Ca;

	mpz_class Ab;
	mpz_class Ap;
    int* ab;
    int* a;
    int* CaInt;
};

#endif /* NBINEQVVEC_H */
