#pragma once
#ifndef GALUAFIELD_H
#include "galuarow.h"
#include <QMap>
#include <QString>
#include <QBitArray>
#include <QtGlobal>
#include <QTime>
#include <qmath.h>
#include <omp.h>
#include <QDebug>
#define GALUAFIELD_H


class GaluaField
{
public:
    GaluaField(int &m);
    ~GaluaField();
    QString GetStrPoly();
    QBitArray GetBitPoly();

    QString BitArrToString(const QBitArray &a);

    QMap<int,GaluaRow> Field;
    QMap<QString, GaluaRow> BitField;

    int** pl;
    int** um;

    int plus(int&, int&);
    int umnozh(int&, int&);
    int del(int&, int&);
    double determ(int**, int);

    int size;

    int powerNum(int num, int power);
    int randIntInField();
    int** RandMatrixSxS(int s);
	int** inverseMatrix(int** matrix, int s);
	int** MatrixMult(int **a, int m, int n, int** b, int n2, int q);

private:
    int m;
    QMap<int,QString> GP;
    QMap<int,QBitArray> GPB;

    void FillField();
    void FillOper();
};

#endif // GALUAFIELD_H
