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
        if(n >= d) break;
        ++i;
    }
    n--;
    d-=k;
    this->m = i;
    this->MakeField();
    Init();
}

Koder::~Koder()
{
    if(H)
    {
        int I = n-k;
        for(int i = 0; i < I; ++i)
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
    if(CStar)
    {
        delete[] CStar;
    }
    if(S)
    {
        delete[] S;
    }
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
    if(!MakeG())qDebug()<<"error in building G";
}

bool Koder::MakeH()
{
    int h = n - k;
    H = new int*[h];
    for(int i = 0; i < h; ++i)
    {
        H[i] = new int[n];
    }
    int c = f->size-1;
    for(int i = 0; i < h; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            H[i][j] = static_cast<int>((i * j)%c) + 1;            
        }
    }

    qDebug()<<"H"<<endl;
    for(int i = 0; i < h; ++i)
    {
        QString s = "";
        for(int j = 0; j < n; ++j)
        {
            s+=QString::number(H[i][j]);
            s+="\t";
        }
        qDebug()<<s;
    }
    qDebug()<<endl;
    return true;
}

bool Koder::MakeG()
{
    int D = n - k;
    G = new int*[k];
    for(int i = 0; i < k; ++i)
    {
        G[i] = new int[n];
    }

    int*** GS = new int**[k];
    for(int i = 0; i < k; ++i)
    {
        GS[i] = new int*[D];
        for(int j = 0; j < D; ++j)
        {
            GS[i][j] = new int[D+1];
        }
    }

    //Podgotovka
    for(int it = 0; it < k; ++it)
    {
        for(int i = 0; i < D; ++i)
        {
            for(int j = 0; j < (D); ++j)
            {
                GS[it][i][j] = H[i][j+k];
            }
            GS[it][i][D] = H[i][it];
        }
    }

    int tmp = 0;
    //GAUSS+
    for(int it = 0; it < k; ++it)
    {
        //прямой проход
        for(int i = 1; i < D; ++i)
        {
            tmp = GS[it][i][0];
            for(int j = 0; j < D+1; ++j)
                GS[it][i][j] = f->pl[GS[it][i][j]][tmp];
        }
        for(int i = 1; i < D; ++i)
        {
            for(int j = i+1; j < D; ++j)
            {
                tmp = f->del(GS[it][j][i], GS[it][i][i]);
                for(int itD = i; itD < D+1; ++itD)
                {
                    GS[it][j][itD] = f->pl[GS[it][j][itD]][f->um[GS[it][i][itD]][tmp]];
                }
            }
        }
        //обратный проход
        for(int i = D-1; i >= 0; --i)
        {
            for(int j = i-1; j >= 0; --j)
            {
                tmp = f->del(GS[it][j][i], GS[it][i][i]);
                for(int itD = 0; itD < D+1; ++itD)
                {
                    GS[it][j][itD] = f->pl[GS[it][j][itD]][f->um[GS[it][i][itD]][tmp]];
                }
            }
        }
        //результат
        for(int i = 0; i < k; ++i)
        {
            if(it!=i) G[it][i] = 0;
            else G[it][i] = 1;
        }
        for(int i = k, j = 0; i < n; ++i, ++j)
        {
            G[it][i] = f->del(GS[it][j][D], GS[it][j][j]);
        }
    }

    for(int i = 0; i < k; ++i)
    {
        for(int j = 0; j < D; ++j)
        {
            delete[] GS[i][j];
        }
        delete[] GS[i];
    }
    delete[] GS;

    qDebug()<<"G"<<endl;
    for(int i = 0; i < k; ++i)
    {
        QString s = "";
        for(int j = 0; j < n; ++j)
        {
            s+=QString::number(G[i][j]);
            s+="\t";
        }
        qDebug()<<s;
    }
    qDebug()<<endl;
    return true;
}

void Koder::encode(QString I)
{
    int* C = new int[n];
    int* E = new int[n];
    CStar = new int[n];
    int* inf = new int[k];

    for(int i = 0; i < k; ++i)
    {
        inf[i] = static_cast<int>(I[i].digitValue());
    }

    for(int i = 0; i < n; ++i)
    {
        int temp = 0;
        for(int j = 0; j < k; ++j)
        {
            temp = f->pl[temp][f->um[inf[j]][G[j][i]]];
        }
        C[i] = temp;
        E[i] = 0;
    }

    for(int i = 0; i<t; ++i)
    {
        E[f->randIntInField()-1] = f->randIntInField();
    }

    QString temp = "";
    for(int i = 0; i < n; ++i)
    {
        CStar[i] = f->pl[C[i]][E[i]];
        temp += QString::number(E[i])+" ";
    }
    qDebug()<<"e: "<<temp;

    temp = "";
    for(int i = 0; i < n; ++i)
    {
        temp += QString::number(CStar[i])+" ";
    }
    qDebug()<<"C*: "<<temp;

    delete[] inf;
    delete[] C;
    delete[] E;
}

QString Koder::decode()
{
    int* lE = new int[n];
    int h = n-k;
    S = new int[h];
    for(int i = 0; i < h; ++i)
    {
        int temp = 0;
        for(int j = 0; j < n; ++j)
        {
            temp = f->pl[temp][f->um[CStar[j]][H[i][j]]];
            lE[j] = 0;
        }
        S[i] = temp;
    }

    QString temp = "";
    for(int i = 0; i < h; ++i)
    {
        temp += QString::number(S[i]);
    }
    qDebug()<<"S: "<<temp;

    int** A = new int*[t];
    for(int i = 0; i < t; ++i)
    {
        A[i] = new int[t+1];
    }

    //Podgotovka
    for(int i = 0; i < t; ++i)
    {
        for(int j = 0; j < t+1; ++j)
        {
            A[i][j] = S[i+j];
        }
    }
    int tmp = 0;
    //прямой проход
    for(int i = 1; i < t; ++i)
    {
        if(A[i][0])
        {
            tmp = A[i][0];
            for(int j = 0; j < t+1; ++j)
                A[i][j] = f->pl[A[i][j]][tmp];
        }
    }
    for(int i = 1; i < t; ++i)
    {
        for(int j = i+1; j < t; ++j)
        {
            tmp = f->del(A[j][i], A[i][i]);
            for(int itD = i; itD < t+1; ++itD)
            {
                A[j][itD] = f->pl[A[j][itD]][f->um[A[i][itD]][tmp]];
            }
        }
    }
    //обратный проход
    for(int i = t-1; i >= 0; --i)
    {
        for(int j = i-1; j >= 0; --j)
        {
            tmp = f->del(A[j][i], A[i][i]);
            for(int itD = 0; itD < t+1; ++itD)
            {
                A[j][itD] = f->pl[A[j][itD]][f->um[A[i][itD]][tmp]];
            }
        }
    }
    //результат
    for(int i = 0; i < t; ++i)
    {
        A[i][t] = f->del(A[i][t], A[i][i]);
        A[i][i] = f->del(A[i][i], A[i][i]);
    }


    qDebug()<<"A";
    for(int i = 0; i < t; ++i)
    {
        temp="";
        for(int j = 0; j < t+1; ++j)
        {
            temp+=QString::number(A[i][j])+" ";
        }
        qDebug()<<temp;
    }

    for(int i = 1; i < n+1; ++i)
    {
        int tempI = 0;
        for(int j = t; j >= 0; --j)
        {
            int x = 1;
//            switch (j) {
//            case t:
//                x = f->powerNum(i, t);
//                break;
//            case 1:
//                x = f->um[i][A[1][t]];
//                break;
//            case 0:
//                x = A[0][t];
//                break;
//            default:
//                x = f->um[f->powerNum(i, j)][A[j][t]];
//                break;
//            }
            if(j==t) x = f->powerNum(i, t);
            else if(j==1) x = f->um[i][A[1][t]];
            else if(j==0) x = A[0][t];
            else x = f->um[f->powerNum(i, j)][A[j][t]];
            tempI = f->pl[tempI][x];
        }
        if(tempI==0)
        {
            lE[i-1] = 1;
        }
    }

    temp="lE: ";
    for(int i = 0; i < n; ++i)
    {
        temp+=QString::number(lE[i])+" ";
    }
    qDebug()<<temp;

    int** E = new int*[t];
    for(int i = 0; i < t; ++i)
    {
        E[i] = new int[t+1];
    }
    //Podgotovka
    for(int i = 0; i < t; ++i)
    {
        for(int j = 0; j < t; ++j)
        {
            E[i][j] = 1;
        }
        E[i][t] = S[i];
    }
    for(int i = 0, ij = 0; i<n; ++i)
    {
        if(lE[i])
        {
            for(int j = 0; j<t; ++j)
            {
                E[j][ij] = H[j][i];
            }
            ++ij;
        }
    }
    tmp = 0;
    //прямой проход
    for(int i = 0; i < t; ++i)
    {
        for(int j = i+1; j < t; ++j)
        {
            tmp = f->del(E[j][i], E[i][i]);
            for(int ij = 0; ij < t+1; ++ij)
                E[j][ij] = f->pl[E[j][ij]][f->um[E[i][ij]][tmp]];
        }
    }
    for(int i = 1; i < t; ++i)
    {
        for(int j = i+1; j < t; ++j)
        {
            tmp = f->del(E[j][i], E[i][i]);
            for(int itD = i; itD < t+1; ++itD)
            {
                E[j][itD] = f->pl[E[j][itD]][f->um[E[i][itD]][tmp]];
            }
        }
    }
    //обратный проход
    for(int i = t-1; i >= 0; --i)
    {
        for(int j = i-1; j >= 0; --j)
        {
            tmp = f->del(E[j][i], E[i][i]);
            for(int itD = 0; itD < t+1; ++itD)
            {
                E[j][itD] = f->pl[E[j][itD]][f->um[E[i][itD]][tmp]];
            }
        }
    }
    //результат
    for(int i = 0; i < t; ++i)
    {
        E[i][t] = f->del(E[i][t], E[i][i]);
        E[i][i] = f->del(E[i][i], E[i][i]);
    }

    qDebug()<<"E";
    for(int i = 0; i < t; ++i)
    {
        temp="";
        for(int j = 0; j < t+1; ++j)
        {
            temp+=QString::number(E[i][j])+" ";
        }
        qDebug()<<temp;
    }

    temp = "";
    for(int i = 0, j = 0; i<n; ++i)
    {
        if(lE[i])
        {
            lE[i] = E[j][t];
            ++j;
        }
        temp+=QString::number(lE[i])+" ";
    }
    qDebug()<<"e_rebuild: "<<temp;

    int* enC = new int[n];
    temp = "";
    for(int i = 0; i < n; ++i)
    {
        enC[i] = f->pl[CStar[i]][lE[i]];
        temp+=QString::number(enC[i])+" ";
    }
    qDebug()<<"decode C: "<<temp;

    temp = "";
    for(int i = 0; i < k; ++i)
    {
        temp+=QString::number(enC[i]);
    }

    for(int i = 0; i < t; ++i)
    {
        delete[] A[i];
        delete[] E[i];
    }
    delete[] A;
    delete[] E;
    delete[] lE;
    delete[] enC;

    return temp;
}
