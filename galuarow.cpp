#include "galuarow.h"

GaluaRow::GaluaRow()
{
}

GaluaRow::GaluaRow(const QBitArray &binary)
{
    this->binary = binary;
}

GaluaRow::GaluaRow(const int &alpha, const QBitArray &binary, const QString &poly)
{
    //this->N = N;
    this->alpha = alpha;
    this->binary = binary;
    this->poly = poly;
}

void GaluaRow::binaryTOpoly()
{
    poly = "";
    int m = binary.size();
    for(int i = m-1; i>=0; --i)
    {
        if(binary[i])
        {
            if(poly.size()!=0)
            poly+=" + ";
            if(i==0)
            {
                poly+="1";
            }
            else
            {
                if(i==1) poly+="x";
                else poly+="x"+QString::number(i);
            }
        }
    }
}

bool GaluaRow::Bitwise()
{
    int m = binary.size();
    QBitArray arr(binary);
    binary[0] = false;
    for(int i = 1; i < m; ++i)
    {
        binary[i] = arr[i-1];
    }
    return arr[m-1];
}


bool GaluaRow::operator ==(GaluaRow const & a)
{
    return a.binary == this->binary ? true : false;
}

QString GaluaRow::toString()
{
    QString s = "";
    s+=QString::number(this->alpha+1)+"\t";
    for(int i = this->binary.size(); i > 0; --i)
    {
        if(this->binary[i-1])s+="1";
        else s+="0";
    }
    s+="\t";
    s+="L" + QString::number(this->alpha) + "\t";
    s+=this->poly;
    return s;
}
