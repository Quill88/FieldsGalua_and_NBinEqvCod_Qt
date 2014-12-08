#include "galuafield.h"

GaluaField::GaluaField(int &m)
{
    this->m = m;
    this->size = static_cast<int>(pow(2, static_cast<double>(m)));

    QBitArray a(3);
    switch (m) {
    case (2):
        GP.insert(2, "x2 + x + 1");
        a[0] = true; a[1] = true; a[2] = true;
        GPB.insert(2, a);
        break;
    case (3):
        GP.insert(3, "x3 + x + 1");
        a.resize(4);
        a[0] = true; a[1] = true; a[2] = false; a[3] = true;
        GPB.insert(3, a);
        break;
    case (4):
        GP.insert(4, "x4 + x + 1");
        a.resize(5);
        a[0] = true; a[1] = true; a[2] = false; a[3] = false; a[4] = true;
        GPB.insert(4, a);
        break;
    case (5):
        GP.insert(5, "x5 + x2 + 1");
        a.resize(6);
        a[0] = true; a[1] = false; a[2] = true; a[3] = false; a[4] = false;
        a[5] = true;
        GPB.insert(5, a);
        break;
    case (6):
        GP.insert(6, "x6 + x + 1");
        a.resize(7);
        a[0] = true;  a[1] = true; a[2] = false; a[3] = false; a[4] = false;
        a[5] = false; a[6] = true;
        GPB.insert(6, a);
        break;
    case (7):
        GP.insert(7, "x7 + x3 + 1");
        a.resize(8);
        a[0] = true;    a[1] = false; a[2] = false; a[3] = false; a[4] = true;
        a[5] = false;   a[6] = false; a[7] = true;
        GPB.insert(7, a);
        break;
    case (8):
        GP.insert(8, "x8 + x4 + x3 + x2 + 1");
        a.resize(9);
        a[0] = true;  a[1] = false; a[2] = true;  a[3] = true; a[4] = true;
        a[5] = false; a[6] = false; a[7] = false; a[8] = true;
        GPB.insert(8, a);
        break;
    case (9):
        GP.insert(9, "x9 + x4 + 1");
        a.resize(10);
        a[0] = true;  a[1] = false; a[2] = false; a[3] = false; a[4] = true;
        a[5] = false; a[6] = false; a[7] = false; a[8] = false; a[9] = true;
        GPB.insert(9, a);
        break;
    case (10):
        GP.insert(10, "x10 + x3 + 1");
        a.resize(11);
        a[0] = true;  a[1] = false; a[2] = false; a[3] = true;  a[4] = false;
        a[5] = false; a[6] = false; a[7] = false; a[8] = false; a[9] = false;
        a[10] = true;
        GPB.insert(10, a);
        break;
    case (11):
        GP.insert(11, "x11 + x2 + 1");
        a.resize(12);
        a[0] = true;    a[1] = false;   a[2] = true;    a[3] = false;    a[4] = false;
        a[5] = false;   a[6] = false;   a[7] = false;   a[8] = false;    a[9] = false;
        a[10] = false;   a[11] = true;
        GPB.insert(11, a);
        break;
    case (12):
        GP.insert(12, "x12 + x6 + x4 + x + 1");
        a.resize(13);
        a[0] = true;    a[1] = true;   a[2] = false;   a[3] = false;    a[4] = true;
        a[5] = false;   a[6] = true;   a[7] = false;   a[8] = false;    a[9] = false;
        a[10] = false;   a[11] = false;   a[12] = true;
        GPB.insert(12, a);
        break;
    case (13):
        GP.insert(13, "x13 + x4 + x3 + x + 1");
        a.resize(14);
        a[0] = true;    a[1] = true;   a[2] = false;   a[3] = true;    a[4] = true;
        a[5] = false;   a[6] = false;  a[7] = false;   a[8] = false;   a[9] = false;
        a[10] = false;  a[11] = false; a[12] = false;  a[13] = true;
        GPB.insert(13, a);
        break;
    case (14):
        GP.insert(14, "x14 + x10 + x6 + x + 1");
        a.resize(15);
        a[0] = true;    a[1] = true;   a[2] = false;   a[3] = false;   a[4] = false;
        a[5] = false;   a[6] = true;   a[7] = false;   a[8] = false;   a[9] = false;
        a[10] = true;   a[11] = false; a[12] = false;  a[13] = false;  a[14] = true;
        GPB.insert(14, a);
        break;
    case (15):
        GP.insert(15, "x15 + x + 1");
        a.resize(16);
        a[0] = true;    a[1] = true;   a[2] = false;   a[3] = false;   a[4] = false;
        a[5] = false;   a[6] = false;  a[7] = false;   a[8] = false;   a[9] = false;
        a[10] = false;  a[11] = false; a[12] = false;  a[13] = false;  a[14] = false;
        a[15] = true;
        GPB.insert(15, a);
        break;
    }

    QTime time = QTime::currentTime();
    qsrand((uint)time.msec());

    FillField();
    FillOper();
}

QString GaluaField::GetStrPoly()
{
    return GP.find(m).value();
}

QBitArray GaluaField::GetBitPoly()
{
    return GPB.find(m).value();
}

/*
QBitArray GaluaField::plus(const QBitArray &a1, const QBitArray &a2)
{
    bool p = false;
    int n = 0, s1 = a1.size(), s2 = a2.size();
    QBitArray A1(a1);
    QBitArray A2(a2);
    if(s1>s2)
    {
        n = s1;
        A2.resize(n);
    }
    else
    {
        n = s2;
        A1.resize(n);
    }
    QBitArray result(n);
    for(int i = 0; i < n; ++i)
    {
        result[i] = p;
        p = false;
        if(A1[i] && A2[i]) p = true;
        else if(A1[i] || A2[i])
        {
            if(result[i])
            {
                p = true;
                result[i] = false;
            }
            else result[i] = true;
        }
    }
    if(p)
    {
        result.resize(n+1);
        result[n] = true;
    }
    return result;
}
*/

void GaluaField::FillField()
{   
    int ch = 0;
    QBitArray arr(m);
    GaluaRow gRow(arr);
    gRow.binaryTOpoly();
    gRow.alpha = ch-1;
    Field.insert(ch, gRow);
    ++ch;
    BitField.insert(BitArrToString(gRow.binary), gRow);

    gRow.binary[0] = true;
    gRow.binaryTOpoly();
    gRow.alpha = ch-1;
    Field.insert(ch, gRow);
    ++ch;
    BitField.insert(BitArrToString(gRow.binary), gRow);

    gRow.overflow = false;

    while (ch<size)
    {
        gRow.overflow = gRow.Bitwise();
        if(!gRow.overflow)
        {
            gRow.binaryTOpoly();
            gRow.alpha = ch-1;
            Field.insert(ch, gRow);
            ++ch;
            BitField.insert(BitArrToString(gRow.binary), gRow);
        }
        else
        {
            while(gRow.overflow)
            {
                arr = gRow.binary ^ this->GetBitPoly();
                gRow.overflow = false;
                arr.resize(m);
                gRow.binary = arr;
            }
            gRow.binaryTOpoly();
            gRow.alpha = ch-1;
            Field.insert(ch, gRow);
            ++ch;
            BitField.insert(BitArrToString(gRow.binary), gRow);
        }
    }
}

QString GaluaField::BitArrToString(const QBitArray &a)
{
    QString s = "";
    for(int i = a.size(); i > 0; --i)
    {
        if(a[i-1])s+="1";
        else s+="0";
    }
    return s;
}

int GaluaField::plus(int &a, int &b)
{
    QBitArray barr = Field.find(a).value().binary ^ Field.find(b).value().binary;
    return BitField.find(BitArrToString(barr)).value().alpha+1;
}

int GaluaField::umnozh(int &a, int &b)
{
    int i = 0;
    if(a!=0 && b!=0)
        i = ((Field[a].alpha + Field[b].alpha) % (size-1))+1;
    return i;
}

int GaluaField::del(int &a, int &b)
{
    int obr = 0;
    if(b==1) obr = 1;
    if(b > 1) obr = ((size-1) - (b-1)) + 1;
    return umnozh(a, obr);
}

void GaluaField::FillOper()
{
    pl = new int* [size];
    um = new int* [size];
    for(int i = 0; i < size; ++i)
    {
        pl[i] = new int [size];
        um[i] = new int [size];
    }

    for(int i = 0; i < size; ++i)
    {
        for(int j = 0; j < size; ++j)
        {
            pl[i][j] = plus(i, j);
            um[i][j] = umnozh(i, j);
        }
    }
}

int** GaluaField::RandMatrixSize()
{
    QTime midnight;
    qsrand(midnight.secsTo(QTime::currentTime()));
    int** a = new int* [size];
    for(int i = 0; i < size; ++i)
    {
        a[i] = new int [size];
    }
    for(int i = 0; i < size; ++i)
    {
        for(int j = 0; j < size; ++j)
        {
            a[i][j] = qrand() % size;
        }
    }
    return a;
}

double** GaluaField::inverseMatrixSize(int **matrix)
{
    int i,j,k,j1;
    double **res = new double *[size];
    double **tmp = new double *[size];
    double **tmp2 = new double *[size];
    double S;
    for(i=0;i<size;++i)
    {
        tmp[i] = new double [2*size];
        tmp2[i] = new double [2*size];
        res[i] = new double [size];
    }

    //Inizialization
    for(i=0;i<size;++i)
    {
        for(j=0;j<size;++j)
        {
            tmp[i][j] = matrix [i][j];
        }
        tmp[i][i+size] = 1;
    }

    //Calculation
    for(k=0;k<size;k++)
    {
        if(tmp[k][k]==0)
        {
            for(i=0;i<size;++i)
            {
               delete [] tmp[i];
               delete [] tmp2[i];
               delete [] res[i];
            }
            delete [] tmp;
            delete [] tmp2;
            delete [] res;

            return NULL;
        }
        for(j=k;j<2*size;j++)
            tmp2[k][j]=tmp[k][j]/tmp[k][k];
        for(i=k+1;i<size;i++)
        {
            for(j=k+1;j<2*size;j++)
                tmp[i][j]=tmp[i][j]-tmp[i][k]*tmp2[k][j];
        }
     }

     for(j=0;j<size;j++)
        res[size-1][j]=tmp2[size-1][j+size];
     for(k=size-2;k>=0;k--)
        for(j1=0;j1<size;j1++)
        {
            S=0;
            for(i=k+1;i<size;i++)
                S+=tmp2[k][i]*res[i][j1];
            res[k][j1]=tmp2[k][j1+size]-S;
        }
    for(i=0;i<size;++i)
    {
        delete [] tmp[i];
        delete [] tmp2[i];
    }
    delete [] tmp;
    delete [] tmp2;

    return res;
}

void GaluaField::MatrixMult(int **a, double b)
{
    for(int i = 0; i < size; ++i)
    {
        for(int j = 0; j < size; ++j)
        {
            a[i][j] = (static_cast<int>(a[i][j] * b)%size);
        }
    }
}

GaluaField::~GaluaField()
{
    if(pl && um)
    {
        for(int i = 0; i < size; ++i)
        {
            delete [] pl[i];
            delete [] um[i];
        }
        delete [] pl;
        delete [] um;
    }
}

double GaluaField::determ(int** Arr, int size)
{
        int i,j;
        double det=0;
        int** matr;
        if(size==1)
        {
                det=Arr[0][0];
        }
        else if(size==2)
        {
                det=Arr[0][0]*Arr[1][1]-Arr[0][1]*Arr[1][0];
        }
        else
        {
                matr=new int*[size-1];
                for(i=0;i<size;++i)
                {
                        for(j=0;j<size-1;++j)
                        {
                                if(j<i)
                                        matr[j]=Arr[j];
                                else
                                        matr[j]=Arr[j+1];
                        }
                        det+=pow((double)-1, (i+j))*determ(matr, size-1)*Arr[i][size-1];
                }               
                delete[] matr;
        }
        return det;
}

int GaluaField::randIntInField()
{
    return qrand() % (size - 1) + 1;
}

int GaluaField::powerNum(int num, int power)
{
    int temp = num;
    for(int i = 1; i < power; ++i)
    {
        temp = um[temp][num];
    }
    return temp;
}
