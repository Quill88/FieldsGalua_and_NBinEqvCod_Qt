#pragma once
#ifndef KODER_H
#define KODER_H

#include <QString>
#include <QStringList>
#include <QDebug>
#include <QBitArray>
#include <QVector>
#include <QVector3D>
#include "galuafield.h"
#include "galuarow.h"
#include <QGlobal.h>
#include "nBinEqvCod.h"

class Koder
{
private:
    bool MakeField();
	bool MakeHdots();
	bool MakeHpoly();
	bool MakeH();
	bool MakeX();
	bool MakeP();
	bool MakeD();
	bool MakeHx();

	bool checkPoint(int x, int y, int z);
	int* gauss(int** Array, int row, int col);
	
public:
    int** H;
	int** Hdots;
	int** Hpoly;
	int** Hx;
	int** X;
	int** inverseX;
	int*  P;
	int*  inverseP;
	int*  D;
	int*  inverseD;
	
	int* Sx;

	QBitArray fixedMch;
    GaluaField *f;
	nBinEqvCod *nBC;
    int n, k, d, m, t, t2, degF;
    Koder(int t, int k);
    void Init();
    ~Koder();

	void test();

    void encode(QString I, int**);
	QString decode(int*);

    GaluaRow getRow(int);
};

#endif // KODER_H
