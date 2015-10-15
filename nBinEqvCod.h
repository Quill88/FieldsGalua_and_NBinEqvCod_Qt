#pragma once
#ifndef NBINEQVCOD_H
#define NBINEQVCOD_H
#include <qmath.h>
#include <QDebug>
#include <QString>
#include "nBinEqvVec.h"

/*nonЦbinary equivalent codes
недвоичный равновесный код*/
class nBinEqvCod
{
public:
	nBinEqvCod(int n, int w, int q);
	~nBinEqvCod();

	QString getEqvVec(int A);

private:
	int M;
	int n;
	int w;
	int q;
	
	nBinEqvVec calc_eVec(int A);

	int fact(const int& n);
};

#endif /* NBINEQVCOD_H */