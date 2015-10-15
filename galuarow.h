#pragma once
#ifndef GALUAROW_H
#include <QBitArray>
#include <QString>
#include <QDataStream>
#include <QTextStream>
#define GALUAROW_H

class GaluaRow
{
public:
    GaluaRow();
    GaluaRow(const QBitArray &binary);
    GaluaRow(const int &alpha, const QBitArray &binary, const QString &poly);
    bool operator==(GaluaRow const &a);
    //inline bool operator==(const GaluaRow& mc1, const GaluaRow& mc2);
    QString toString();
    //int m;
    bool Bitwise(); //сдвиг битового представления на 1, возвращает произошло переполнение или нет
    //QTextStream operator<<(const GaluaRow &row);
    //QDataStream &operator>>(QDataStream &in, GaluaRow &row);
    //QDataStream &operator<<(const GaluaRow &row);   //(QDataStream &out, const GaluaRow &row);
    //QDataStream &operator>>(GaluaRow &row);

    QBitArray binary; //битовое представление
    QString poly; //полином
    int alpha;
    //int N;
    bool overflow;

    void binaryTOpoly();
};

#endif // GALUAROW_H

