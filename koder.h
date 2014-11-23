#ifndef KODER_H
#define KODER_H

#include "galuafield.h"
#include "galuarow.h"

class Koder
{
private:
    GaluaField *f;
    int MakeField();
    //Koder(int, int, int);
public:
    int n, k, d, m;
    Koder(int, int);
    void Init();
    ~Koder();

    GaluaRow getRow(int);
};

#endif // KODER_H
