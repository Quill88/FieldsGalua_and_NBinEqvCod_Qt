#include "koder.h"

Koder::Koder(int t, int k) : t(t), k(k)
{
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
	H = nullptr;
	X = nullptr;
	P = nullptr;
	D = nullptr;
	Hx = nullptr;
	CStar = nullptr;
	S = nullptr;
    Init();
}

Koder::~Koder()
{
    if(H)
    {
        delete[] H[0];
		delete[] H[1];
		delete[] H[2];
        delete[] H;
    }
    if(X)
    {
        for(int i = 0; i < k; ++i)
        {
            delete[] X[i];
        }
        delete[] X;
    }
	if (P)
	{
		delete[] P;
	}
	if (D)
	{
		delete[] D;
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
	QTime time = QTime::currentTime();
	qsrand((uint)time.msec());
	
	fixedMch = QBitArray(10);
	fixedMch[0] = true; fixedMch[1] = true; fixedMch[2] = true;
	fixedMch[3] = true;

	degF = 4;

    if(!MakeH())qDebug()<<"error in building H";
    if(!MakeX())qDebug()<<"error in building X";
}

bool Koder::MakeH()
{
    H = new int*[3];
    for(int i = 0; i < 3; ++i)
    {
        H[i] = new int[n];
    }
	
	QVector<QVector3D> allPoints;
	if (checkPoint(1, 0, 0)) allPoints.push_back(QVector3D(1, 0, 0));
	for (int x = 0; x <= n; ++x)
		if (checkPoint(x, 1, 1)) allPoints.push_back(QVector3D(x, 1, 1));
	for (int x = 0; x <= n; ++x)
		for (int y = 0; y <= n; ++y)
			if (checkPoint(x, y, 1)) allPoints.push_back(QVector3D(x, y, 1));

	qDebug() << "allPoints";
	for(int i = 0; i < allPoints.count();++i)
	{
		qDebug() << i << "\t {" << allPoints.at(i).x()
			<< ", " << allPoints.at(i).y()
			<< ", " << allPoints.at(i).z() << "}";
	}
	qDebug() << endl;

    for(int i = 0; i < n; ++i)
    {
		int r = qrand() % (allPoints.count());
		H[0][i] = allPoints.at(r).x();
		H[1][i] = allPoints.at(r).y();
		H[2][i] = allPoints.at(r).z();
    }

    qDebug()<<"H"<<endl;
    for(int i = 0; i < 3; ++i)
    {
        QString s = "";
        for(int j = 0; j < n; ++j)
        {
			if (H[i][j] >= 0 && H[i][j] <= 9) s += " ";
            s+=QString::number(H[i][j]);
            s+=" ";
        }
        qDebug()<<s;
    }
    qDebug()<<endl;
    return true;
}

bool Koder::MakeX()
{
	if (degF<3 || degF>9)
	{
		return false;
	}

	X = new int*[n];
	for (int i = 0; i < n; ++i)
	{
		X[i] = new int[3];
	}

	QVector<QVector3D> allFunctions;

	qDebug() << "Functions";
	for (int i = degF; i >= 0; --i)
		for (int j = degF - i; j >= 0; --j)
		{
			allFunctions.push_back(QVector3D(i, j, degF - i - j));
			qDebug() << i << " " << j << " "<< degF - i - j;
		}


	for (int i = 0; i < n; ++i)
	{
		int r = qrand() % allFunctions.count();
		X[i][0] = allFunctions.at(r).x();
		X[i][1] = allFunctions.at(r).y();
		X[i][2] = allFunctions.at(r).z();
	}

	qDebug() << endl << "X" << endl;
	for (int i = 0; i < n; ++i)
	{
		qDebug() << X[i][0] << " " << X[i][1] << " " << X[i][2];
	}
	qDebug() << endl;

    return true;
}

int dotsPoly[10][3] = { { 3, 0, 0 }, { 2, 1, 0 }, { 2, 0, 1 }, { 1, 1, 1 }, 
{ 1, 2, 0 }, { 1, 0, 2 }, { 0, 3, 0 }, { 0, 2, 1 }, { 0, 1, 2 }, { 0, 0, 3 } };
bool Koder::checkPoint(int x, int y, int z)
{
	int res = 0;
	for (int i = 0; i < 10; ++i)
	{
		if (fixedMch[i])
		{
			int mch = 1;
			for (int j = 0; j < dotsPoly[i][0]; ++j)
				mch = f->um[mch][x];
			for (int j = 0; j < dotsPoly[i][1]; ++j)
				mch = f->um[mch][y];
			for (int j = 0; j < dotsPoly[i][2]; ++j)
				mch = f->um[mch][z];
			res = f->pl[res][mch];
		}
	}
	return !(res);
}

bool Koder::MakeHx()
{
	return true;
}

void Koder::encode(QString I)
{
    //int* C = new int[n];
    //int* E = new int[n];
    //CStar = new int[n];
    //int* inf = new int[k];

    //for(int i = 0; i < k; ++i)
    //{
    //    inf[i] = static_cast<int>(I[i].digitValue());
    //}

    //for(int i = 0; i < n; ++i)
    //{
    //    int temp = 0;
    //    for(int j = 0; j < k; ++j)
    //    {
    //        temp = f->pl[temp][f->um[inf[j]][G[j][i]]];
    //    }
    //    C[i] = temp;
    //    E[i] = 0;
    //}

    //for(int i = 0; i<t; ++i)
    //{
    //    E[f->randIntInField()-1] = f->randIntInField();
    //}

    //QString temp = "";
    //for(int i = 0; i < n; ++i)
    //{
    //    CStar[i] = f->pl[C[i]][E[i]];
    //    temp += QString::number(E[i])+" ";
    //}
    //qDebug()<<"e: "<<temp;

    //temp = "";
    //for(int i = 0; i < n; ++i)
    //{
    //    temp += QString::number(CStar[i])+" ";
    //}
    //qDebug()<<"C*: "<<temp;

    //delete[] inf;
    //delete[] C;
    //delete[] E;
}

QString Koder::decode()
{
//    int* lE = new int[n];
//    int h = n-k;
//    S = new int[h];
//    for(int i = 0; i < h; ++i)
//    {
//        int temp = 0;
//        for(int j = 0; j < n; ++j)
//        {
//            temp = f->pl[temp][f->um[CStar[j]][H[i][j]]];
//            lE[j] = 0;
//        }
//        S[i] = temp;
//    }
//
//    QString temp = "";
//    for(int i = 0; i < h; ++i)
//    {
//        temp += QString::number(S[i]);
//    }
//    qDebug()<<"S: "<<temp;
//
//    int** A = new int*[t];
//    for(int i = 0; i < t; ++i)
//    {
//        A[i] = new int[t+1];
//    }
//
//    //Podgotovka
//    for(int i = 0; i < t; ++i)
//    {
//        for(int j = 0; j < t+1; ++j)
//        {
//            A[i][j] = S[i+j];
//        }
//    }
//    int tmp = 0;
//    //прямой проход
//    for(int i = 1; i < t; ++i)
//    {
//        if(A[i][0])
//        {
//            tmp = A[i][0];
//            for(int j = 0; j < t+1; ++j)
//                A[i][j] = f->pl[A[i][j]][tmp];
//        }
//    }
//    for(int i = 1; i < t; ++i)
//    {
//        for(int j = i+1; j < t; ++j)
//        {
//            if(A[i][i]!=0)
//            {
//                tmp = f->del(A[j][i], A[i][i]);
//                for(int itD = i; itD < t+1; ++itD)
//                {
//                    A[j][itD] = f->pl[A[j][itD]][f->um[A[i][itD]][tmp]];
//                }
//            }
//        }
//    }
//    //обратный проход
//    for(int i = t-1; i >= 0; --i)
//    {
//        for(int j = i-1; j >= 0; --j)
//        {
//            if(A[i][i]!=0)
//            {\
//                tmp = f->del(A[j][i], A[i][i]);
//                for(int itD = 0; itD < t+1; ++itD)
//                {
//                    A[j][itD] = f->pl[A[j][itD]][f->um[A[i][itD]][tmp]];
//                }
//            }
//        }
//    }
//    //результат
//    for(int i = 0; i < t; ++i)
//    {
//        A[i][t] = f->del(A[i][t], A[i][i]);
//        A[i][i] = f->del(A[i][i], A[i][i]);
//    }
//
//
//    qDebug()<<"A";
//    for(int i = 0; i < t; ++i)
//    {
//        temp="";
//        for(int j = 0; j < t+1; ++j)
//        {
//            temp+=QString::number(A[i][j])+" ";
//        }
//        qDebug()<<temp;
//    }
//
//    for(int i = 1; i < n+1; ++i)
//    {
//        int tempI = 0;
//        for(int j = t; j >= 0; --j)
//        {
//            int x = 1;
////            switch (j) {
////            case t:
////                x = f->powerNum(i, t);
////                break;
////            case 1:
////                x = f->um[i][A[1][t]];
////                break;
////            case 0:
////                x = A[0][t];
////                break;
////            default:
////                x = f->um[f->powerNum(i, j)][A[j][t]];
////                break;
////            }
//            if(j==t) x = f->powerNum(i, t);
//            else if(j==1) x = f->um[i][A[1][t]];
//            else if(j==0) x = A[0][t];
//            else x = f->um[f->powerNum(i, j)][A[j][t]];
//            tempI = f->pl[tempI][x];
//        }
//        if(tempI==0)
//        {
//            lE[i-1] = 1;
//        }
//    }
//
//    temp="lE: ";
//    for(int i = 0; i < n; ++i)
//    {
//        temp+=QString::number(lE[i])+" ";
//    }
//    qDebug()<<temp;
//
//    int** E = new int*[t];
//    for(int i = 0; i < t; ++i)
//    {
//        E[i] = new int[t+1];
//    }
//    //Podgotovka
//    for(int i = 0; i < t; ++i)
//    {
//        for(int j = 0; j < t; ++j)
//        {
//            E[i][j] = 1;
//        }
//        E[i][t] = S[i];
//    }
//    for(int i = 0, ij = 0; i<n; ++i)
//    {
//        if(lE[i])
//        {
//            for(int j = 0; j<t; ++j)
//            {
//                E[j][ij] = H[j][i];
//            }
//            ++ij;
//        }
//    }
//    tmp = 0;
//    //прямой проход
//    for(int i = 0; i < t; ++i)
//    {
//        for(int j = i+1; j < t; ++j)
//        {
//            if(E[i][i]!=0)
//            {
//                tmp = f->del(E[j][i], E[i][i]);
//                for(int ij = 0; ij < t+1; ++ij)
//                    E[j][ij] = f->pl[E[j][ij]][f->um[E[i][ij]][tmp]];
//            }
//
//        }
//    }
//    for(int i = 1; i < t; ++i)
//    {
//        for(int j = i+1; j < t; ++j)
//        {
//            if(E[i][i]!=0)
//            {
//                tmp = f->del(E[j][i], E[i][i]);
//                for(int itD = i; itD < t+1; ++itD)
//                {
//                    E[j][itD] = f->pl[E[j][itD]][f->um[E[i][itD]][tmp]];
//                }
//            }
//        }
//    }
//    //обратный проход
//    for(int i = t-1; i >= 0; --i)
//    {
//        for(int j = i-1; j >= 0; --j)
//        {
//            if(E[i][i] != 0)
//            {
//                tmp = f->del(E[j][i], E[i][i]);
//                for(int itD = 0; itD < t+1; ++itD)
//                {
//                    E[j][itD] = f->pl[E[j][itD]][f->um[E[i][itD]][tmp]];
//                }
//            }
//        }
//    }
//    //результат
//    for(int i = 0; i < t; ++i)
//    {
//        E[i][t] = f->del(E[i][t], E[i][i]);
//        E[i][i] = f->del(E[i][i], E[i][i]);
//    }
//
//    qDebug()<<"E";
//    for(int i = 0; i < t; ++i)
//    {
//        temp="";
//        for(int j = 0; j < t+1; ++j)
//        {
//            temp+=QString::number(E[i][j])+" ";
//        }
//        qDebug()<<temp;
//    }
//
//    temp = "";
//    for(int i = 0, j = 0; i<n; ++i)
//    {
//        if(lE[i])
//        {
//            lE[i] = E[j][t];
//            ++j;
//        }
//        temp+=QString::number(lE[i])+" ";
//    }
//    qDebug()<<"e_rebuild: "<<temp;
//
//    int* enC = new int[n];
//    temp = "";
//    for(int i = 0; i < n; ++i)
//    {
//        enC[i] = f->pl[CStar[i]][lE[i]];
//        temp+=QString::number(enC[i])+" ";
//    }
//    qDebug()<<"decode C: "<<temp;
//
//    temp = "";
//    for(int i = 0; i < k; ++i)
//    {
//        temp+=QString::number(enC[i]);
//    }
//
//    for(int i = 0; i < t; ++i)
//    {
//        delete[] A[i];
//        delete[] E[i];
//    }
//    delete[] A;
//    delete[] E;
//    delete[] lE;
//    delete[] enC;
//
    return "dec";
}
