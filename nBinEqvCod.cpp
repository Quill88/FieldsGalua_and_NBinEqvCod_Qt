#include "nBinEqvCod.h"


nBinEqvCod::nBinEqvCod(int n, int w, int q) : n(n), w(w), q(q)
{
    omp_set_dynamic(1);
    omp_set_num_threads(8);
    mpz_ui_pow_ui(qw.get_mpz_t(), (q - 1), w);
    //qw = qPow((q - 1), w);
	//qw = Bigint(q-1).pow(w);

	M = mpz_class(qw) * comb(n, w); // factorial(n) / (factorial(w) * factorial(n - w));

    //code = new nBinEqvVec[M];
    code = nullptr;

	qDebug() << "Creating non_binary equivalent codes";
	qDebug() << "n:" << n << "\tw: " << w << "\tq:" << q << "\tM: " << M.get_str().c_str();

	//#pragma omp parallel for
	//for (int i = 0; i < M; ++i)
	//{
	//	calc_eVec(i, &code[i]);
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

mpz_class nBinEqvCod::getM() const
{
	return M;
}

QString nBinEqvCod::getStringEqvVec(mpz_class A)
{
    if (0 <= A && A < M)
    {
        if (code)
            return nullptr;//code[A].ToStr();
        else
        {
			QString s = "Err";
			nBinEqvVec* v = getEqvVecByNum(A);
			if (v)
			{
				s = v->ToStr();
				delete v;
			}
			return s;
        }
    }
    else return "A not in 0<=A<M";
}

mpz_class nBinEqvCod::comb(long n, long k)
{
	if (n == k) return 1;
	else if (k == 0) return 1;
	else if (n < k) return 0;
	else
	{
		mpz_class pf, f;
		if (k > (n - k))
		{
			pf = partial_fact(k + 1, n);
			//f = factorial(from_bigint_to_int(n - k));
			f = FactTree(n - k);
		}
		else
		{
			pf = partial_fact(n - k + 1, n);
			//f = factorial(from_bigint_to_int(k));
			f = FactTree(k);
		}

		return pf / f;
	}
}

void nBinEqvCod::test()
{
	//qDebug() << comb(14, 2);
	//qDebug() << comb(13, 2);
}

void nBinEqvCod::calc_eVec(mpz_class A, nBinEqvVec* v)
{
	v->A = A;
	v->n = n;
	v->w = w;

    v->ab = new int[n];
    v->a = new int[w];
    v->CaInt = new int[n];

#pragma omp parallel sections
{
#pragma omp section
	{
		mpz_class Ab = A / qw;
		v->Ab = Ab.get_ui();
		mpz_class x = Ab;
        int l = 0;

		for (int i = 0; i < n; ++i)
		{
			mpz_class b = comb((n - i - 1), (w - l));
			if (b > x)
				v->ab[n - i - 1] = 0;
			else
			{
				v->ab[n - i - 1] = 1;
				x -= b;
				l++;
			}
		}
	}
#pragma omp section
	{
		mpz_class Ap = A % qw;
		v->Ap = Ap;
		mpz_class x = Ap;
        int l = 0;
		while (l < w)
		{
			mpz_class temp_a = (x % (q - 1)) + 1;
			v->a[l] = temp_a.get_ui();
			x = x / (q - 1);
			l++;
		}
	}
}
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

void nBinEqvCod::calc_eVec_byStr(nBinEqvVec* v)
{
	v->n = n;
	v->w = w;

    v->ab = new int[n];
    v->a = new int[w];
    v->CaInt = new int[n];

	QStringList list = v->Ca.split(' ', QString::SkipEmptyParts);
	if (list.count() != n)
		qDebug() << "Split err!";

	for (int i = 0, j = 0; i < n; ++i)
	{
		v->CaInt[i] = list[i].toInt();
		if (v->CaInt[i] != 0)
		{
			v->ab[i] = 1;
			v->a[j] = v->CaInt[i];
			++j;
		}
		else
			v->ab[i] = 0;
	}

	v->Ab = 0;
	for (int i = 0, l = 0; i < n; ++i)
	{
		if (v->ab[n - i - 1] == 1)
		{
            mpz_class j = (comb((n - i - 1), (w - l)));
			if (j >= 0)
			{
				v->Ab += j;
				++l;
			}
		}
	}

	v->Ap = 0;
	for (int i = 0; i < w; ++i)
    {
        mpz_class res;
        mpz_ui_pow_ui(res.get_mpz_t(), (q - 1), i);
        v->Ap += res * (v->a[i] - 1);
	}

    mpz_class temp;
    mpz_ui_pow_ui(temp.get_mpz_t(), (q - 1), w);
    v->A = v->Ab * temp + v->Ap;
}

nBinEqvCod::~nBinEqvCod()
{
	if (code)
		delete[] code;
}

nBinEqvVec* nBinEqvCod::getEqvVecByNum(mpz_class A)
{
    if (A >= 0 && A < M)
    {
        if (code)
            return nullptr; //&code[A];
        else
        {
			nBinEqvVec* v = new nBinEqvVec();
			calc_eVec(A, v);
			return v;
        }
    }
    else return nullptr;
}

nBinEqvVec* nBinEqvCod::getEqvVecByCa(QString Ca)
{
	if (Ca == "")
		return nullptr;
	if (mapNBEV.contains(Ca))
		return mapNBEV[Ca];
	else
	{
		nBinEqvVec* v = new nBinEqvVec(Ca);
		calc_eVec_byStr(v);
		return v ? v : nullptr;
	}
}
