#include "nBinEqvCod.h"


nBinEqvCod::nBinEqvCod(int n, int w, int q) : n(n), w(w), q(q)
{
	//omp_set_dynamic(1);
	//omp_set_num_threads(8);
	//
	//qw = qPow((q - 1), w);

	///*иру 1*/
	//M = qw * (fact(n)) / (fact(w)*fact(n - w));

	//code = new nBinEqvVec[M];
	//
	code = nullptr;
	M = 2027025000;

	qDebug() << "Creating non_binary equivalent codes";
	qDebug() << "n:" << n <<"\tw: " << w <<"\tq:"<<q <<"\tM: " <<M;

	//#pragma omp parallel for
	//for (int i = 0; i < M; ++i)
	//{
	//	calc_eVec(i, code[i]);
	//	mapNBEV.insert(code[i].Ca, &code[i]);
	//}
	//--------------------------------------------------
	//QTime midnight(0, 0, 0);
	//qsrand(midnight.secsTo(QTime::currentTime()));
	//code = new nBinEqvVec[3];
	//calc_eVec(qrand() % M, code[0]);
	//calc_eVec(qrand() % M, code[1]);
	//calc_eVec(qrand() % M, code[2]);
}

int nBinEqvCod::fact(const int& n)
{
	if (n >= 0)
		return n ? (n * fact(n - 1)) : 1;
	else
		return n;
}		

int nBinEqvCod::getM() const
{
	return M;
}

QString nBinEqvCod::getStringEqvVec(int A)
{
	if (0 <= A && A < M)
	{
		return code[A].ToStr();
	}
	else return "A not in 0<=A<M";
}

void nBinEqvCod::calc_eVec(int A, nBinEqvVec* v)
{
	v->A = A;
	v->n = n;
	v->w = w;

	/*иру 2*/
	int Ab = A / qw;
	int Ap = A % qw;

	v->ab = new int[n];
	v->a = new int[w];
	v->CaInt = new int[n];
	
	/*иру 3*/
	int x = Ab;
	int l = 0;

	for (int i = 0; i < n; ++i)
	{
		int b = fact(n - i - 1) / (fact(w - l)*fact((n - i - 1) - (w - l)));
		if (b > x) v->ab[n - i - 1] = 0;
		else
		{
			v->ab[n - i - 1] = 1;
			x -= b;
			l++;
		}
	}

	/*иру 4*/
	x = Ap;
	l = 0;
	while (l < w)
	{
		v->a[l] = (x % (q - 1)) + 1;
		x = x / (q - 1);
		l++;
	}

	/*иру 5*/
	for (int i = 0, l = 0; i < n; ++i)
	{
		if (v->ab[i])
		{
			v->Ca += QString::number(v->a[l]);
			v->CaInt[i] = v->a[l];
			l++;
		}
		else
		{
			v->CaInt[i] = 0;
			v->Ca += "0";
		}
		v->Ca += " ";
	}
}


nBinEqvCod::~nBinEqvCod()
{
	if (code) delete[] code;
}

nBinEqvVec* nBinEqvCod::getEqvVecByNum(int A)
{
	if (A >= 0 && A < M)
	{
		nBinEqvVec* v = new nBinEqvVec();
		calc_eVec(A, v);
		return v;
	}
	else return nullptr;
}


nBinEqvVec* nBinEqvCod::getEqvVecByCa(QString Ca)
{
	if (mapNBEV.contains(Ca)) return mapNBEV[Ca];
	else return nullptr;
}