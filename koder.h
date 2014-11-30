#ifndef KODER_H
#define KODER_H

#include <QString>
#include <QDebug>
#include "galuafield.h"
#include "galuarow.h"

class Koder
{
private:
    bool MakeField();
    bool MakeH();
    bool MakeG();
    //Koder(int, int, int);
public:
    int** H;
    int** G;
    GaluaField *f;
    int n, k, d, m, t;
    Koder(int t, int k);
    void Init();
    ~Koder();

    GaluaRow getRow(int);
};

#endif // KODER_H
