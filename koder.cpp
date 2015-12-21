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
	inverseP = nullptr;
	D = nullptr;
	inverseD = nullptr;
	Hx = nullptr;
	Sx = nullptr;

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
	if (P)
	{
		delete[] P;
		delete[] inverseP;
	}
	if (D)
	{
		delete[] D;
		delete[] inverseD;
	}
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
	if (Sx) delete[] Sx;
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

	nBC = new nBinEqvCod(2, 2, (n + 1));


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
	inverseP = new int[n];

	int** tempInverseP = nullptr;

	for (int i = 0; i < n; ++i)
		P[i] = f->randIntInField();

	int** tempP = new int*[n];
	for (int i = 0; i < n; ++i)
	{
		tempP[i] = new int[n];
		tempP[i][i] = P[i];
	}

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			if (i != j) tempP[i][j] = 0;
		}
	}

	tempInverseP = f->inverseMatrix(tempP, n);

	for (int i = 0; i < n; ++i)
	{
		inverseP[i] = tempInverseP[i][i];
	}

	qDebug() << endl << "P";
	QString s = "";
	QString s2 = "";
	for (int i = 0; i < n; ++i)
	{
		s += QString::number(P[i]) + " ";
		s2 += QString::number(inverseP[i]) + " ";
	}
	s.remove(s.size() - 1, 1);
	s2.remove(s2.size() - 1, 1);
	qDebug() << s << endl << "inverseP" << endl << s2;

	for (int i = 0; i < n; ++i)
	{
		delete[] tempP[i];
		delete[] tempInverseP[i];
	}
	delete[] tempP; delete[] tempInverseP;

	return true;
}

bool Koder::MakeD()
{
	D = new int[n];
	inverseD = new int[n];

	for (int i = 0; i < n; ++i)
		D[i] = i + 1;

	std::random_shuffle(D, D + n);

	for (int i = 0; i < n; ++i)
	{
		inverseD[D[i] - 1] = i + 1;
	}

	qDebug() << endl << "D";
	QString s = "";
	QString s2 = "";
	for (int i = 0; i < n; ++i)
	{
		s += QString::number(D[i]) + " ";
		s2 += QString::number(inverseD[i]) + " ";
	}
	s.remove(s.size() - 1, 1);
	s2.remove(s2.size() - 1, 1);
	qDebug() << s;
	qDebug() << "inverseD" << endl << s2;

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

	qDebug() << endl << "H" << endl;
	for (int i = 0; i < k; ++i)
	{
		QString s = "";
		for (int j = 0; j < n; j++)
		{
			if (H[i][j] >= 0 && H[i][j] <= 9) s += " ";
			s += QString::number(H[i][j]) + " ";
		}
		qDebug() << s;
	}
	qDebug() << endl;

	return true;
}

bool Koder::MakeX()
{
	X = f->RandMatrixSxS(k);
	inverseX = f->inverseMatrix(X, k);

	qDebug() << endl << "X                                        inverseX" << endl;
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
		qDebug() << s << "|" << s2;
	}
	qDebug() << endl;

	int** tempX = new int*[k];
	for (int i = 0; i < k; ++i)
	{
		tempX[i] = new int[k];
	}
	tempX = f->MatrixMult(X, k, k, inverseX, k, k);

	//qDebug() << endl << "X cheking" << endl;
	//for (int i = 0; i < k; ++i)
	//{
	//	QString s = "";
	//	for (int j = 0; j < k; j++)
	//	{
	//		if (tempX[i][j] >= 0 && tempX[i][j] <= 9) s += " ";
	//		s += QString::number(tempX[i][j]) + " ";
	//	}
	//	qDebug() << s;
	//}
	//qDebug() << endl;

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

void Koder::encode(QString I, int** Hx)
{
	//QStringList splited = I.split(" ");
	//if (splited.size() != k) qDebug() << splited.size() << " wrong count!";

	int* tempI = new int[n];
	//for (int i = 0; i < k; ++i)
	//{
	//	tempI[i] = splited[i].toInt();
	//}

	if (I.length() != (n / 2))
	{
		qDebug() << "I wrong length!";
	}
	else
	{
		for (int i = 0, l = 1; i < I.length(); ++i)
		{
			ushort ch = I.at(i).toLatin1();
			tempI[l] = nBC->getEqvVecByNum(static_cast<int>(ch))->CaInt[0];
			tempI[l + 1] = nBC->getEqvVecByNum(static_cast<int>(ch))->CaInt[1];
			l += 2;
		}
		tempI[0] = 0;

		QString s = "";
		for (int i = 0; i < n; ++i)
			s += QString::number(tempI[i]) + " ";
		s.remove(s.size() - 1, 1);
		qDebug() << endl << "I: " << s;

		//Вместо транспонирования и умножения сразу умножаем в виде строка*строка
		Sx = new int[k];
		for (int i = 0; i < k; ++i)
		{
			Sx[i] = 0;
			for (int j = 0; j < n; ++j)
				Sx[i] = f->pl[Sx[i]][f->um[tempI[j]][Hx[i][j]]];
		}

		s = "";
		for (int i = 0; i < k; ++i)
			s += QString::number(Sx[i]) + " ";
		s.remove(s.size() - 1, 1);

		qDebug() << endl << "Sx: " << s;

	}

	delete[] tempI;
}

QString Koder::decode(int* Sx)
{
	int* Cx = new int[n];
	for (int i = 0; i < n; ++i)
	{
		if (i >= k) Cx[i] = 1;
		else Cx[i] = 0;
	}

	int** tempArr = new int*[k];
	for (int i = 0; i < k; ++i)
	{
		tempArr[i] = new int[k + 1];
	}

	for (int i = 0; i < k; ++i)
	{
		for (int j = 0; j < k; j++)
		{
			tempArr[i][j] = Hx[i][j];
		}
	}
	for (int i = 0; i < k; ++i)
	{
		tempArr[i][k] = Sx[i];
		for (int j = k; j < n; j++)
		{
			tempArr[i][k] = f->pl[tempArr[i][k]][Hx[i][j]];
		}
	}

	/*qDebug() << endl << "tempArr" << endl;
	for (int i = 0; i < k; ++i)
	{
	QString s = "";
	for (int j = 0; j < (k+1); j++)
	{
	if (tempArr[i][j] >= 0 && tempArr[i][j] <= 9) s += " ";
	s += QString::number(tempArr[i][j]) + " ";
	}
	qDebug() << s;
	}
	qDebug() << endl;*/

	int* tempSx = gauss(tempArr, k, k + 1);

	for (int i = 0; i < k; ++i)
	{
		Cx[i] = tempSx[i];
	}
	delete[] tempSx;


	QString s = "";
	for (int i = 0; i < n; ++i)
		s += QString::number(Cx[i]) + " ";
	s.remove(s.size() - 1, 1);
	qDebug() << endl << "Cx: " << s;


	int* tempCx = new int[n];
	
	for (int j = 0; j < n; ++j)
	{
		tempCx[j] = Cx[inverseD[j] - 1];
	}

	for (int j = 0; j < n; ++j)
		tempCx[j] = f->um[tempCx[j]][inverseP[j]];


	s = "";
	for (int i = 0; i < n; ++i)
		s += QString::number(tempCx[i]) + " ";
	s.remove(s.size() - 1, 1);
	qDebug() << endl << "Cx`: " << s;


	tempSx = new int[k];
	for (int i = 0; i < k; ++i)
	{
		tempSx[i] = 0;
		for (int j = 0; j < n; ++j)
			tempSx[i] = f->pl[tempSx[i]][f->um[tempCx[j]][H[i][j]]];
	}

	s = "";
	for (int i = 0; i < k; ++i)
		s += QString::number(tempSx[i]) + " ";
	s.remove(s.size() - 1, 1);
	qDebug() << endl << "Sx`: " << s;

	for (int i = 0; i < k; ++i)
	{
		delete[] tempArr[i];
	}
	delete[] tempArr;
	delete[] tempSx;
	delete[] tempCx;
	delete[] Cx;

	return "dec";
}

int* Koder::gauss(int** Array, int row, int col)
{
	int** resultPH = new int*[row]; //создание массива для прямого хода
	for (int i = 0; i < row; i++)
	{
		resultPH[i] = new int[col];
	}
	for (int i = 0; i < row; i++) //копирование данных, что бы не испортить входные данные
		for (int j = 0; j < col; j++)
			resultPH[i][j] = Array[i][j];
	int* resultOH = new int[row]; //массив для результата обратного хода
	for (int i = 0; i < row; i++)
		resultOH[i] = 0;

	//Прямой ход
	for (int i = 0; i < row; i++)
	{
		int k = i; //выбор наибольшего элемента в столбце
		for (int j = i; j < row; j++)
		{
			if (resultPH[j][i] > resultPH[k][i]) k = j;
		}
		if (k != i) //замена строк если это необходимо 
		{
			for (int j = i; j < col; j++)
			{
				int t = resultPH[i][j];
				resultPH[i][j] = resultPH[k][j];
				resultPH[k][j] = t;
			}
		}

		for (int j = i + 1; j < row; j++) //вычитание текущей строки из всех нижних, умножив на соответствующий коэффициент
		{
			int c = f->del(resultPH[j][i], resultPH[i][i]);
			for (k = i; k < col; k++)
			{
				resultPH[j][k] = f->pl[resultPH[j][k]][f->um[resultPH[i][k]][c]];
			}
		}
	}

	//Обратный ход
	for (int i = row - 1; i >= 0; i--)
	{
		resultOH[i] = resultPH[i][col - 1];
		for (int j = 0; j < col - 1; j++)
		{
			if (j != i)
			{
				resultOH[i] = f->pl[resultOH[i]][f->um[resultOH[j]][resultPH[i][j]]];
			}
		}
		resultOH[i] = f->del(resultOH[i], resultPH[i][i]);
	}

	//qDebug() << endl << "result" << endl;
	//for (int i = 0; i < row; ++i)
	//{
	//	QString s = "";
	//	for (int j = 0; j < (row + 1); j++)
	//	{
	//		if (resultPH[i][j] >= 0 && resultPH[i][j] <= 9) s += " ";
	//		s += QString::number(resultPH[i][j]) + " ";
	//	}
	//	qDebug() << s;
	//}
	//qDebug() << endl;

	for (int i = 0; i < row; ++i)
	{
		delete[] resultPH[i];
	}
	delete[] resultPH;

	return resultOH;
}