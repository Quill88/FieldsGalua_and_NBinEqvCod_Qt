#pragma once
#ifndef GALUAFIELD_H
#include "galuarow.h"
#include <QMap>
#include <QString>
#include <QBitArray>
#include <QGlobal.h>
#include <QTime>
#include <qmath.h>
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
    int** RandMatrixSize();
    double** inverseMatrixSize(int** matrix);
    void MatrixMult(int**, double);

    //QBitArray plus(const QBitArray &a1, const QBitArray &a2);
private:
    int m;\
    QMap<int,QString> GP;
    QMap<int,QBitArray> GPB;

    void FillField();
    void FillOper();
};

#endif // GALUAFIELD_H
