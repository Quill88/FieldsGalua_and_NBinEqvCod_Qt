#include "koder.h"

Koder::Koder(int t, int k)
{
    this->t = t;
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
    Init();
}

//Koder::Koder(int n, int k, int d) : n(n), k(k), d(d)
//{ }
Koder::~Koder()
{
    if(H)
    {
        for(int i = 0; i < 2*t; ++i)
        {
            delete[] H[i];
        }
        delete[] H;
    }
    if(G)
    {
        for(int i = 0; i < k; ++i)
        {
            delete[] G[i];
        }
        delete[] G;
    }
    delete f;
}

bool Koder::MakeField()
{
    if(m>15 || m<2)
    {
        return false;
    }
    f = new GaluaField(m);
    return 1;
}

GaluaRow Koder::getRow(int i)
{
    GaluaRow r = f->Field[i];
    return r;
}

void Koder::Init()
{
    if(!MakeH())qDebug()<<"error in building H";
//    if(!MakeG())qDebug()<<"error in building G";
}

bool Koder::MakeH()
{
    H = new int*[2*t];
    for(int i = 0; i < 2*t; ++i)
    {
        H[i] = new int[n-1];
    }
    int c = f->size-1;
    for(int i = 0; i < 2*t; ++i)
    {
        for(int j = 0; j < n-1; ++j)
        {
            H[i][j] = static_cast<int>((i * j)%c) + 1;
            //qDebug()<<H[i][j];
        }
    }
    return true;
}

bool Koder::MakeG()
{
    G = new int*[k];
    for(int i = 0; i < k; ++i)
    {
        G[i] = new int[n-1];
    }

    int*** GS = new int**[k];
    for(int i = 0; i < k; ++i)
    {
        GS[i] = new int*[2*t];
        for(int j = 0; j < 2*t; ++j)
        {
            GS[i][j] = new int[d];
        }
    }

    //Podgotovka
    for(int it = 0; it < k; ++it)
    {
        for(int i = 0; i < 2*t; ++i)
        {
            for(int j = 0; j < (d-1); ++j)
            {
                GS[it][i][j] = H[i][j+k];
            }
            GS[it][i][d-1] = H[i][it];
        }
    }

    int xx[d];
    int N = 2*t;
    int tmp = 0;
    //GAUSS
    for(int it = 0; it < k; ++it)
    {
        for(int i = 0; i < N; ++i)
        {
            tmp = GS[it][i][i];
            for(int j = N; j >= i; j--)
                GS[it][i][j] = f->del(GS[it][i][j], tmp);
            for(int j = i+1; j < N; ++j)
            {
                tmp = GS[it][j][i];
                for(int l = N; l >= i; l--)
                    GS[it][j][l] = f->pl[GS[it][j][l]][f->um[tmp][GS[it][i][k]]];
            }
        }
        xx[N-1] = GS[it][N-1][N];
        for(int i = N-2; i >= 0; i--)
        {
            xx[i] = GS[it][i][N];
            for(int j = i+1; j < N; j++)
                xx[i] = f->pl[xx[i]][f->um[GS[it][i][j]][xx[i]]];
        }
        QString s = "";
        for(int i = 0; i < N; ++i)
        {
            for(int j = 0; j<d; ++j)
                s+=QString::number(GS[it][i][j]) + " ";
            qDebug() << s;
            s = "";
        }
        for(int i = 0; i < d; ++i)
        {
            s+=QString::number(xx[i]) + " ";
        }
        qDebug() << endl << s << endl;
    }

    for(int i = 0; i < k; ++i)
    {
        for(int j = 0; j < 2*t; ++j)
        {
            delete[] GS[i][j];
        }
        delete[] GS[i];
    }
    delete[] GS;
    return true;
}
