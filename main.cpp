#include "galuafield.h"
#include "galuarow.h"
#include "koder.h"
#include "nBinEqvCod.h"
#include <QTextStream>
#include <QString>
#include <QBitArray>
#include <QMap>
#include <QMapIterator>
#include <QFile>
#include <QTime>
#include <QDebug>

int main()
{
	QTime timer;

	QFile file("out.txt");
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
		return 0;

	QTextStream cout(stdout);
	QTextStream fout(&file);
	QTextStream cin(stdin);

	int t = 5; int k = 10;
	Koder Kdr(t, k);

	cout << "G_field" << endl;
	cout << "(n, k, d)" << endl;
	cout << "(" << Kdr.n << ", " << Kdr.k << ", " << Kdr.d << ")" << endl << endl;


	int n = 2; int w = 2; int q = Kdr.n;

	timer = QTime::currentTime();

	nBinEqvCod nBC(n, w, q);



	QString I = "5 10 7 2 9 12 29 20 18 9";
	Kdr.encode(I);

	//-----------Запись в файл-----------------------------------------------------------

	fout << "A\tbin\tab\t\ta\tCa" << endl;
	for (int i = 0; i < nBC.getM(); ++i)
	{
		fout << nBC.getStringEqvVec(i) << endl;
	}
	fout << endl;
	QMapIterator<int, GaluaRow> i(Kdr.f->Field);
	fout << "no\t" << "bin\t" << "L\t" << "poly\t" << endl;
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
	//-----------------------------------------------------------------------------------

	cout << endl << "Time spend: " << timer.msecsTo(QTime::currentTime()) << " msecs" << endl;

	return 0;
}

/*
	QTime timer;
	QFile file("out.txt");
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
	return 0;

	QTextStream cout(stdout);
	QTextStream fout(&file);

	QTextStream cin(stdin);

	timer = QTime::currentTime();
	cout<<timer.toString("hh:mm:ss.zzz")<<endl;

	cout << "Enter t: ";
	cout.flush();
	int t;
	try {
	cin >> t;
	} catch (...) {
	cout<<"Error!!!";
	}
	cout << "Enter k: ";
	cout.flush();
	int k;
	try {
	cin >> k;
	} catch (...) {
	cout<<"Error!!!";
	}
	cout<<endl;

	Koder K(t, k);

	cout<<"(n, k, d)"<<endl;
	cout<<"("<<K.n<<", "<<K.k<<", "<<K.d<<")"<<endl<<endl;

	timer = QTime::currentTime();
	cout<<timer.toString("hh:mm:ss.zzz")<<endl;



	cout << "Enter I (lenth="<<k<<"): ";
	cout.flush();
	QString I;
	try {
	cin >> I;
	} catch (...) {
	cout<<"Error!!!";
	}
	cout<<endl;

	K.encode(I);
	QString temp = K.decode();
	qDebug()<<endl<<"I: "<<temp<<endl;


	QMapIterator<int, GaluaRow> i(K.f->Field);
	fout << "no\t" << "bin\t" << "L\t" << "poly\t" << endl;
	while (i.hasNext())
	{
	i.next();
	fout << i.key() << ":\t" <<K.f->BitArrToString(i.value().binary)
	<< "\tL" << QString::number(i.value().alpha)
	<< "\t" << i.value().poly << endl;
	}

	QString fout2 = "  ";
	fout << endl << "+" << endl;
	fout << "  " << "\t";
	fout2 += " \t";
	for(int i = 0; i < K.f->size; ++i)
	{
	fout << QString::number(i) << ".\t";
	fout2 += QString::number(i) + ".\t";
	}
	fout << endl;
	fout2 += "\n";
	for(int i = 0; i < K.f->size; ++i)
	{
	fout << QString::number(i) + "|\t";
	fout2 += QString::number(i) + "|\t";
	for(int j = 0; j < K.f->size; ++j)
	{
	fout << QString::number(K.f->pl[i][j]) << "\t";
	fout2 += QString::number(K.f->um[i][j])+ "\t";
	}
	fout << endl;
	fout2 += "\n";
	}
	fout << endl << "*" << endl << fout2;

	fout << endl << "/" << endl << "\t";
	for(int i = 0; i < K.f->size; ++i)
	{
	fout << QString::number(i) << ".\t";
	}
	fout << endl;
	for(int i = 0; i < K.f->size; ++i)
	{
	fout << QString::number(i) + "|\t";
	for(int j = 0; j < K.f->size; ++j)
	{
	fout << QString::number(K.f->del(i,j)) << "\t";
	}
	fout << endl;
	}

	timer = QTime::currentTime();
	cout<<timer.toString("hh:mm:ss.zzz")<<endl;

	*/