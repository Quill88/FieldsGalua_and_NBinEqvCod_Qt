#include "koder.h"

Koder::Koder(int t, int k) : t(t), k(k)
{
	d = k + (2 * t + 1);
	int i = 1;
	while (true)
	{
		n = static_cast<int>(pow(2, static_cast<double>(i)));
		if (n >= d) break;
		++i;
	}
	n--;
	d -= k;
	this->m = i;
	this->MakeField();
	H = nullptr;
	X = nullptr;
	inverseX = nullptr;
	Hdots = nullptr;
	Hpoly = nullptr;
	P = nullptr;
	D = nullptr;
	Hx = nullptr;
	CStar = nullptr;
	S = nullptr;

	nBC = nullptr;

	Init();
}

Koder::~Koder()
{
	if (Hdots)
	{
		delete[] Hdots[0];
		delete[] Hdots[1];
		delete[] Hdots[2];
		delete[] Hdots;
	}
	if (Hpoly)
	{
		for (int i = 0; i < k; ++i)
		{
			delete[] Hpoly[i];
		}
		delete[] Hpoly;
	}
	if (H)
	{
		for (int i = 0; i < k; ++i)
		{
			delete[] H[i];
		}
		delete[] H;
	}
	if (X)
	{
		for (int i = 0; i < k; ++i)
		{
			delete[] X[i];
			delete[] inverseX[i];
		}
		delete[] X;
		delete[] inverseX;
	}
	if (P) delete[] P;
	if (D) delete[] D;
	if (Hx)
	{
		for (int i = 0; i < k; ++i)
		{
			delete[] Hx[i];
		}
		delete[] Hx;
	}
	delete f;
	if (nBC) delete nBC;
	if (CStar) delete[] CStar;
	if (S) delete[] S;
}

bool Koder::MakeField()
{
	if (m > 15 || m < 2)
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


	nBC = new nBinEqvCod(2, 1, (n + 1));


	if (!MakeHdots())qDebug() << "error in building Hdots";
	if (!MakeHpoly())qDebug() << "error in building Hpoly";
	if (!MakeH())qDebug() << "error in building H";
	if (!MakeX())qDebug() << "error in building X";
	if (!MakeP())qDebug() << "error in building P";
	if (!MakeD())qDebug() << "error in building D";
	if (!MakeHx())qDebug() << "error in building Hx";
}

bool Koder::MakeHdots()
{
	Hdots = new int*[3];
	for (int i = 0; i < 3; ++i)
	{
		Hdots[i] = new int[n];
	}

	QVector<QVector3D> allPoints;
	if (checkPoint(1, 0, 0)) allPoints.push_back(QVector3D(1, 0, 0));
	for (int x = 0; x <= n; ++x)
		if (checkPoint(x, 1, 1)) allPoints.push_back(QVector3D(x, 1, 1));
	for (int x = 0; x <= n; ++x)
		for (int y = 0; y <= n; ++y)
			if (checkPoint(x, y, 1)) allPoints.push_back(QVector3D(x, y, 1));

	qDebug() << "allPoints";
	for (int i = 0; i < allPoints.count(); ++i)
	{
		qDebug() << i << "\t {" << allPoints.at(i).x()
			<< ", " << allPoints.at(i).y()
			<< ", " << allPoints.at(i).z() << "}";
	}
	qDebug() << endl;

	for (int i = 0; i < n; ++i)
	{
		int r = qrand() % (allPoints.count());
		Hdots[0][i] = allPoints.at(r).x();
		Hdots[1][i] = allPoints.at(r).y();
		Hdots[2][i] = allPoints.at(r).z();
	}

	qDebug() << "Hdots" << endl;
	for (int i = 0; i < 3; ++i)
	{
		QString s = "";
		for (int j = 0; j < n; ++j)
		{
			if (Hdots[i][j] >= 0 && Hdots[i][j] <= 9) s += " ";
			s += QString::number(Hdots[i][j]);
			s += " ";
		}
		qDebug() << s;
	}
	qDebug() << endl;
	return true;
}

bool Koder::MakeHpoly()
{
	if (degF < 3 || degF>9)
	{
		return false;
	}

	Hpoly = new int*[k];
	for (int i = 0; i < k; ++i)
	{
		Hpoly[i] = new int[3];
	}

	QVector<QVector3D> allFunctions;

	qDebug() << "Functions";
	for (int i = degF; i >= 0; --i)
		for (int j = degF - i; j >= 0; --j)
		{
			allFunctions.push_back(QVector3D(i, j, degF - i - j));
			qDebug() << i << " " << j << " " << degF - i - j;
		}

	for (int i = 0; i < k; ++i)
	{
		int r = qrand() % allFunctions.count();
		Hpoly[i][0] = allFunctions.at(r).x();
		Hpoly[i][1] = allFunctions.at(r).y();
		Hpoly[i][2] = allFunctions.at(r).z();
	}

	qDebug() << endl << "Hpoly" << endl;
	for (int i = 0; i < k; ++i)
	{
		qDebug() << Hpoly[i][0] << " " << Hpoly[i][1] << " " << Hpoly[i][2];
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

bool Koder::MakeP()
{
	P = new int[n];

	for (int i = 0; i < n; ++i)
		P[i] = f->randIntInField();

	qDebug() << endl << "P";
	QString s = "";
	for (int i = 0; i < n; ++i)
		s += QString::number(P[i]) + " ";
	qDebug() << s;

	return true;
}

bool Koder::MakeD()
{
	D = new int[n];

	for (int i = 0; i < n; ++i)
		D[i] = i + 1;

	std::random_shuffle(D, D + n);

	qDebug() << endl << "D";
	QString s = "";
	for (int i = 0; i < n; ++i)
		s += QString::number(D[i]) + " ";
	qDebug() << s;

	return true;
}

bool Koder::MakeH()
{
	H = new int*[k];
	for (int i = 0; i < k; ++i)
	{
		H[i] = new int[n];
	}
	
	for (int i = 0; i < k; ++i)
		for (int j = 0; j < n; ++j)
		{
			H[i][j] = f->powerNum(Hdots[0][j], Hpoly[i][0]);
			H[i][j] = f->um[H[i][j]][f->powerNum(Hdots[1][j], Hpoly[i][1])];
			H[i][j] = f->um[H[i][j]][f->powerNum(Hdots[2][j], Hpoly[i][2])];
		}
	return true;
}

bool Koder::MakeX()
{
	X = f->RandMatrixSxS(k);
	inverseX = f->inverseMatrix(X, k);

	qDebug() << endl << "X                                                                                    inverseX"<< endl;
	for (int i = 0; i < k; ++i)
	{
		QString s = "";
		QString s2 = "";
		for (int j = 0; j < k; j++)
		{
			if (X[i][j] >= 0 && X[i][j] <= 9) s += "  ";
			if (X[i][j] >= 10 && X[i][j] <= 99) s += " ";
			s += QString::number(X[i][j]) + " ";
			if (inverseX[i][j] >= 0 && inverseX[i][j] <= 9) s2 += "  ";
			if (inverseX[i][j] >= 10 && inverseX[i][j] <= 99) s2 += " ";
			s2 += QString::number(inverseX[i][j]) + " ";
		}
		qDebug() << s <<"|"<<s2;
	}
	qDebug() << endl;

	int** tempX = new int*[k];
	for (int i = 0; i < k; ++i)
	{
		tempX[i] = new int[k];
	}
	tempX = f->MatrixMult(X, k, k, inverseX, k, k);

	qDebug() << endl << "X cheking" << endl;
	for (int i = 0; i < k; ++i)
	{
		QString s = "";
		for (int j = 0; j < k; j++)
		{
			if (tempX[i][j] >= 0 && tempX[i][j] <= 9) s += " ";
			s += QString::number(tempX[i][j]) + " ";
		}
		qDebug() << s;
	}
	qDebug() << endl;

	if (tempX)
	{
		for (int i = 0; i < k; ++i)
		{
			delete[] tempX[i];
		}
		delete[] tempX;
	}
	return true;
}

bool Koder::MakeHx()
{
	Hx = f->MatrixMult(X, k, k, H, k, n);

	for (int i = 0; i < k; ++i)
		for (int j = 0; j < n; ++j)
			Hx[i][j] = f->um[Hx[i][j]][P[j]];

	int** tempHx = new int*[k];
	for (int i = 0; i < k; ++i)
	{
		tempHx[i] = new int[n];
	}

	for (int i = 0; i < k; ++i)
		for (int j = 0; j < n; ++j)
		{
			tempHx[i][j] = Hx[i][D[j] - 1];
		}

	for (int i = 0; i < k; ++i)
	{
		delete[] Hx[i];
	}
	delete[] Hx;

	Hx = tempHx;
	tempHx = nullptr;

	qDebug() << endl << "Hx" << endl;
	for (int i = 0; i < k; ++i)
	{
		QString s = "";
		for (int j = 0; j < n; j++)
		{
			if (Hx[i][j] >= 0 && Hx[i][j] <= 9) s += " ";
			s += QString::number(Hx[i][j]) + " ";
		}
		qDebug() << s;
	}
	qDebug() << endl;

	return true;
}

void Koder::encode(QString I)
{
	//QStringList splited = I.split(" ");
	//if (splited.size() != k) qDebug() << splited.size() << " wrong count!";

	int* tempI = new int[n];
	//for (int i = 0; i < k; ++i)
	//{
	//	tempI[i] = splited[i].toInt();
	//}
	
	if (I.length() != (n/2)) qDebug() << "I wrong length!";

	for (int i = 0, l = 1; i < I.length(); ++i)
	{
		ushort ch = I.at(i).toLatin1();
		tempI[l] = nBC->getEqvVecByNum(static_cast<int>(ch))->CaInt[0];
		tempI[l+1] = nBC->getEqvVecByNum(static_cast<int>(ch))->CaInt[1];
		l += 2;
	}
	tempI[0]=0;

	QString s = "";
	for (int i = 0; i < n; ++i)
		s += QString::number(tempI[i]) + " ";
	s.remove(s.size() - 1, 1);
	qDebug() << endl << "I: " << s;
	
	//Вместо транспонирования и умножения сразу умножаем в виде строка*строка
	CStar = new int[k];
	for (int i = 0; i < k; ++i)
	{
		CStar[i] = 0;
		for (int j = 0; j < n; ++j)
			CStar[i] = f->pl[CStar[i]][f->um[tempI[j]][Hx[i][j]]];
	}

	s = "";
	for (int i = 0; i < k; ++i)
		s += QString::number(CStar[i]) + " ";
	s.remove(s.size() - 1, 1);
	
	qDebug() << endl << "Encode arr: " << s;

	delete[] tempI;
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
	//            temp = f->pl[temp][f->um[CStar[j]][Hdots[i][j]]];
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
	//                E[j][ij] = Hdots[j][i];
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
