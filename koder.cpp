#include "koder.h"

Koder::Koder(int t, int k)
{
    this->k = k;
    d = k + (2*t + 1);
    int i = 1;
    while (true)
    {
        n = static_cast<int>(pow(2, static_cast<double>(i)));
        if(n > d) break;
        ++i;
    }
    d-=k;
    this->m = i;
    this->MakeField();
}

//Koder::Koder(int n, int k, int d) : n(n), k(k), d(d)
//{ }
Koder::~Koder()
{
    delete f;
}

int Koder::MakeField()
{
    if(m>15 || m<2)
    {
        return -1;
    }
    f = new GaluaField(m);
    return 0;
}

GaluaRow Koder::getRow(int i)
{
    GaluaRow r = f->Field[i];
    return r;
}

void Koder::Init()
{

}
