#pragma once
#ifndef NBINEQVCOD_H
#define NBINEQVCOD_H
#include <qmath.h>
#include <QDebug>
#include <QString>
#include "nBinEqvVec.h"

/*non�binary equivalent codes
���������� ����������� ���*/
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