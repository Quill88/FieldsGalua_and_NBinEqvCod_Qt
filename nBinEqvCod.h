#pragma once
#ifndef NBINEQVCOD_H
#define NBINEQVCOD_H
#include <qmath.h>
#include <QDebug>
#include <QString>
#include <QMap>
#include <QException>
#include "nBinEqvVec.h"
#include <omp.h>
#include <QGlobal.h>
#include <QTime>

/*nonЦbinary equivalent codes
недвоичный равновесный код*/
class nBinEqvCod
{
public:
	nBinEqvCod(int n, int w, int q);
	~nBinEqvCod();

	QString getStringEqvVec(int A);
	int getM() const;
	nBinEqvVec* getEqvVecByNum(int A);
	nBinEqvVec* getEqvVecByCa(QString Ca);
	
	void calc_eVec(int i, nBinEqvVec* v);

private:
	int M;
	int n;
	int w;
	int q;
	int qw;

	nBinEqvVec* code;
	QMap<QString, nBinEqvVec*> mapNBEV;

	int fact(const int& n);
};

#endif /* NBINEQVCOD_H */