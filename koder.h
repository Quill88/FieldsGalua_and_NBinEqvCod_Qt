#pragma once
#ifndef KODER_H
#define KODER_H

#include <QString>
#include <QDebug>
#include <QBitArray>
#include <QVector>
#include <QVector3D>
#include "galuafield.h"
#include "galuarow.h"
#include <QGlobal.h>

class Koder
{
private:
    bool MakeField();
    bool MakeH();
    bool MakeX();
	bool MakeHx();

	bool checkPoint(int x, int y, int z);

    int* CStar;
    int* S;

public:
    int** H;
	int** Hx;
    int** X;
	int* P;
	int* D;

	QBitArray fixedMch;
    GaluaField *f;
    int n, k, d, m, t, degF;
    Koder(int t, int k);
    void Init();
    ~Koder();

    void encode(QString I);
    QString decode();

    GaluaRow getRow(int);
};

#endif // KODER_H
