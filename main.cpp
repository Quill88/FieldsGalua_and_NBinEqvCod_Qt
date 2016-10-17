//#include "galuafield.h"
//#include "galuarow.h"
#include "koder.h"
//#include "nBinEqvCod.h"
#include <QTextStream>
#include <QString>
#include <QBitArray>
#include <QMap>
#include <QMapIterator>
#include <QFile>
#include <QFileInfo>
#include <QTime>
#include <QDebug>

#include <gmpxx.h>

int main()
{
	QTime timer;
	timer = QTime::currentTime();

    QTextStream cout(stdout);
//    QTextStream cin(stdin);

    int t = 4; int k = 6;
    //int t = 20; int k = 15;
	Koder Kdr(t, k);

	cout << endl << "G_field" << endl;
	cout << "(n, k, d)" << endl;
    cout << "(" << Kdr.n << ", " << Kdr.k << ", " << Kdr.d << ")" << endl;

    Kdr.mult_test(1024,false);

//--------------------------Testing encode NBEC--------------------------------
//    Kdr.test_encode_nonBEC();

//------------------------------koder test-------------------------------------
//    Kdr.test(3, 2, 3);

//-------------------------------gmp test--------------------------------------
//    Kdr.Fact_test();


	//QString I = "AbC";
	//Kdr.encode(I, Kdr.Hx);

	//cout<<endl<<"Message: "<<Kdr.decode(Kdr.Sx)<<endl;
	
    //-----------Запись в файл-------------------------------------------------

//	QFile file2("nBinEqvCod.txt");
//	if (!file2.open(QIODevice::WriteOnly | QIODevice::Text))
//		return 0;
//	QTextStream fout3(&file2);
//	fout3 << "A\tbin       Ab\t[ab]\t\tAp\t[a]\tCa" << endl;
//	for (int i = 0; i < Kdr.nBC->getM(); ++i)
//	{
//		nBinEqvVec* v = Kdr.nBC->getEqvVecByNum(i);
//		fout3 << v->ToStr() << endl;
//		delete v;
//	}
//	fout3 << endl << endl << endl;
//	fout3 << "A\tbin       Ab\t[ab]\t\tAp\t[a]\tCa" << endl;
//	for (int i = 0; i < Kdr2.nBC->getM(); ++i)
//	{
//		nBinEqvVec* v = Kdr2.nBC->getEqvVecByNumFromCode(i);
//		fout3 << v->ToStr() << endl;
//	}
//	file2.close();

    QFile file("out.txt");
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
        return 0;

    QTextStream fout(&file);

	QMapIterator<int, GaluaRow> i(Kdr.f->Field);
	fout << "no\t" << "bin     " << "L\t" << "poly\t" << endl;
	while (i.hasNext())
	{
		i.next();
		fout << i.key() << ":\t" << Kdr.f->BitArrToString(i.value().binary)
			<< "\tL" << QString::number(i.value().alpha)
			<< "\t" << i.value().poly << endl;
	}

	QString fout2 = "  ";
	fout << endl << "+" << endl;
	fout << "  " << "\t";
	fout2 += " \t";
	for (int i = 0; i < Kdr.f->size; ++i)
	{
		fout << QString::number(i) << ".\t";
		fout2 += QString::number(i) + ".\t";
	}
	fout << endl;
	fout2 += "\n";
	for (int i = 0; i < Kdr.f->size; ++i)
	{
		fout << QString::number(i) + "|\t";
		fout2 += QString::number(i) + "|\t";
		for (int j = 0; j < Kdr.f->size; ++j)
		{
			fout << QString::number(Kdr.f->pl[i][j]) << "\t";
			fout2 += QString::number(Kdr.f->um[i][j]) + "\t";
		}
		fout << endl;
		fout2 += "\n";
	}
	fout << endl << "*" << endl << fout2;

	fout << endl << "/" << endl << "\t";
	for (int i = 0; i < Kdr.f->size; ++i)
	{
		fout << QString::number(i) << ".\t";
	}
	fout << endl;
	for (int i = 0; i < Kdr.f->size; ++i)
	{
		fout << QString::number(i) + "|\t";
		for (int j = 0; j < Kdr.f->size; ++j)
		{
			fout << QString::number(Kdr.f->del(i, j)) << "\t";
		}
		fout << endl;
    }
    //-------------------------------------------------------------------------

    cout << endl << "Time spend: " << timer.msecsTo(QTime::currentTime())
         << " msecs" << endl;

    cout << endl << "Output file: " << QFileInfo(file).absoluteFilePath()
         << endl;

	file.close();
	return 0;
}
