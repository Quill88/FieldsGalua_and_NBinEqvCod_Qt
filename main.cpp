#include "galuafield.h"
#include "galuarow.h"
#include <QTextStream>
#include <QString>
#include <QBitArray>
#include <QMap>
#include <QMapIterator>
#include <QFile>
#include <QTime>

int main()
{
    QTime timer;
    QFile file("out.txt");
        if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
            return 0;

    QTextStream cout(stdout);
    QTextStream fout(&file);

    QTextStream cin(stdin);

    timer = QTime::currentTime();
    cout<<timer.toString("hh:mm:ss.zzz")<<endl;

    cout << "Vvedite t: ";
    cout.flush();
    int t;
    try {
        cin >> t;
    } catch (...) {
        cout<<"Error!!!";
    }
    cout << "Vvedite k: ";
    cout.flush();
    int k;
    try {
        cin >> k;
    } catch (...) {
        cout<<"Error!!!";
    }
    cout<<endl;
    int N = k + (2*t + 1);
    int m = 1, n = 0;
    while (true)
    {
        n = static_cast<int>(pow(2, static_cast<double>(m)));
        if(n > N) break;
        ++m;
    }

    //cout << "M (2 <= M <= 15): "<<QString::number(m);
    cout.flush();

    //cin >> m;
    cout<<"(n, k, d)"<<endl;
    cout<<"("<<n<<", "<<k<<", "<<(2*t+1)<<")"<<endl<<endl;

    if(m>15 || m<2)
    {
        return 0;
    }

    GaluaField f(m);

    timer = QTime::currentTime();
    //cout<<"Field build at "<<QString::number(timer.elapsed())<<endl;
    cout<<timer.toString("hh:mm:ss.zzz")<<endl;

    QMapIterator<int, GaluaRow> i(f.Field);
    fout << "no\t" << "bin\t" << "L\t" << "poly\t" << endl;
    while (i.hasNext())
    {
        i.next();
        fout << i.key() << ":\t" << f.BitArrToString(i.value().binary)
             << "\tL" << QString::number(i.value().alpha)
             << "\t" << i.value().poly << endl;
    }


    QString fout2 = "  ";
    fout << "+" << endl;
    fout << "  " << "\t";
    fout2 += " \t";
    for(int i = 0; i < f.size; ++i)
    {
        fout << QString::number(i) << "\t";
        fout2 += QString::number(i) + "\t";
    }
    fout << endl;
    fout2 += "\n";
    for(int i = 0; i < f.size; ++i)
    {
        fout << QString::number(i) + "\t";
        fout2 += QString::number(i) + "\t";
        for(int j = 0; j < f.size; ++j)
        {
            fout << QString::number(f.pl[i][j]) << "\t";
            fout2 += QString::number(f.um[i][j])+ "\t";
        }
        fout << endl;
        fout2 += "\n";
    }
    fout << endl << "*" << endl << fout2;


    timer = QTime::currentTime();
    cout<<timer.toString("hh:mm:ss.zzz")<<endl;

    fout<<endl<<"Rand matrix:"<<endl;
    int** j = f.RandMatrixSize();    

//    cout<<f.determ(j, f.size)<<endl;
//    if(!f.determ(j, f.size))
//    {
//        while(!f.determ(j, f.size))
//        {
//            for(int i = 0; i < f.size; ++i)
//            {
//                delete [] j[i];
//            }
//            delete [] j;

//            j = f.RandMatrixSize();
//            cout<<f.determ(j, f.size)<<endl;
//        }
//    }
//    fout<<"Determ:"<<f.determ(j, f.size)<<endl;

    for(int k = 0; k < f.size; ++k)
    {
        for(int n = 0; n < f.size; ++n)
        {
            fout << QString::number( j[k][n] ) << "\t";
        }
        fout<<endl;
    }

    double** iMatr = f.inverseMatrixSize(j);
    fout<<endl<<"Inverse Matrix"<<endl;
    for(int k = 0; k < f.size; ++k)
    {
        for(int n = 0; n < f.size; ++n)
        {
            fout << QString::number( iMatr[k][n] ) << "\t";
        }
        fout<<endl;
    }

    fout<<endl<<"Matrix mult number "<<QString::number(m)<<endl;
    f.MatrixMult(j, static_cast<double>(m));
    for(int k = 0; k < f.size; ++k)
    {
        for(int n = 0; n < f.size; ++n)
        {
            fout << QString::number(j[k][n]) << "\t";
        }
        fout<<endl;
    }

    for(int i = 0; i < f.size; ++i)
    {
        delete [] j[i];
        delete [] iMatr[i];
    }
    delete [] j;
    delete [] iMatr;

    timer = QTime::currentTime();
    cout<<timer.toString("hh:mm:ss.zzz")<<endl;

    //system("pause");
    return 0;
}
