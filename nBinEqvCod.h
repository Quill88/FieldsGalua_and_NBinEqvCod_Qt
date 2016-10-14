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
#include <QtGlobal>
#include <QTime>
#include "smart_fact.h"

/*nonbinary equivalent codes
недвоичный равновесный код*/
class nBinEqvCod
{
public:
	nBinEqvCod(int n, int w, int q);
	~nBinEqvCod();

	QString getStringEqvVec(mpz_class);
	mpz_class getM() const;
	//long get_long_M() const;
	nBinEqvVec* getEqvVecByNum(mpz_class);
	nBinEqvVec* getEqvVecByCa(QString);

	void test();

private:
	mpz_class M;
    int n;
    int w;
    int q;
	mpz_class qw;

	nBinEqvVec* code;
	QMap<QString, nBinEqvVec*> mapNBEV;

	//int from_bigint_to_int(Bigint) const;

	void calc_eVec(mpz_class, nBinEqvVec*);
	void calc_eVec_byStr(nBinEqvVec*);

    mpz_class comb(long, long);
};

#endif /* NBINEQVCOD_H */
