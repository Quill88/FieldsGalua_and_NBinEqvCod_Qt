#include "koder.h"

Koder::Koder(int t, int k) : k(k), t(t)
{
    d = k + (2 * t + 1);
    t2 = 2 * t;
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
    f = nullptr;
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
        for (int i = 0; i < t2; ++i)
        {
            delete[] Hpoly[i];
        }
        delete[] Hpoly;
    }
    if (H)
    {
        for (int i = 0; i < t2; ++i)
        {
            delete[] H[i];
        }
        delete[] H;
    }
    if (X)
    {
        for (int i = 0; i < t2; ++i)
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
        for (int i = 0; i < t2; ++i)
        {
            delete[] Hx[i];
        }
        delete[] Hx;
    }
    if (nBC) delete nBC;
    if (Sx)
        delete[] Sx;
    if (f) delete f;
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

    this->MakeField();

    fixedMch = QBitArray(10);
    fixedMch[0] = true; fixedMch[1] = true; fixedMch[2] = true;
    fixedMch[3] = true;

    degF = 3;

    nBC = new nBinEqvCod(n, t, (n + 1));

    //    if (!MakeHdots())qDebug() << "error in building Hdots";
    //    if (!MakeHpoly())qDebug() << "error in building Hpoly";
    //    if (!MakeH())qDebug() << "error in building H";
    //    if (!MakeX())qDebug() << "error in building X";
    //    if (!MakeP())qDebug() << "error in building P";
    //    if (!MakeD())qDebug() << "error in building D";
    //    if (!MakeHx())qDebug() << "error in building Hx";
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

    qDebug() << endl << "allPoints";
    for (int i = 0; i < allPoints.count(); ++i)
    {
        qDebug() << i << "\t {" << allPoints.at(i).x()
                 << ", " << allPoints.at(i).y()
                 << ", " << allPoints.at(i).z() << "}";
    }

    for (int i = 0; i < n; ++i)
    {
        int r = qrand() % (allPoints.count());
        Hdots[0][i] = allPoints.at(r).x();
        Hdots[1][i] = allPoints.at(r).y();
        Hdots[2][i] = allPoints.at(r).z();
    }

    qDebug() << endl << "Hdots";
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

    return true;
}

int degFField[7][2] = { { 3, 10 }, { 4, 15 }, { 5, 21 }, { 6, 28 }, { 7, 36 }, { 8, 45 }, { 9, 55 } };
bool Koder::MakeHpoly()
{
    if (degFField[0][1] < t2)
        for (int i = 1; i < 7; ++i)
        {
            degF = degFField[i][0];
            if (degFField[i][1] >= t2) i += 7;
        }

    if (degF < 3 || degF>9)
    {
        return false;
    }

    Hpoly = new int*[t2];
    for (int i = 0; i < t2; ++i)
    {
        Hpoly[i] = new int[3];
    }

    QVector<QVector3D> allFunctions;

    qDebug() << endl << "Functions";
    for (int i = degF; i >= 0; --i)
        for (int j = degF - i; j >= 0; --j)
        {
            allFunctions.push_back(QVector3D(i, j, degF - i - j));
            qDebug() << i << " " << j << " " << degF - i - j;
        }

    for (int i = 0; i < t2; ++i)
    {
        int r = qrand() % allFunctions.count();
        Hpoly[i][0] = allFunctions.at(r).x();
        Hpoly[i][1] = allFunctions.at(r).y();
        Hpoly[i][2] = allFunctions.at(r).z();
    }

    qDebug() << endl << "Hpoly";
    for (int i = 0; i < t2; ++i)
    {
        qDebug() << Hpoly[i][0] << " " << Hpoly[i][1] << " " << Hpoly[i][2];
    }

    return true;
}

//TODO
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

    for (int i = 0; i < n; ++i)
        P[i] = i + 1;

    std::random_shuffle(P, P + n);

    //-----------------------------------------------------------------
    //P[0] = 2; P[1] = 4; P[2] = 7; P[3] = 5; P[4] = 1; P[5] = 6; P[6] = 3;

    for (int i = 0; i < n; ++i)
    {
        inverseP[P[i] - 1] = i + 1;
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
    qDebug() << s;
    qDebug() << "inverseP" << endl << s2;

    return true;
}

bool Koder::MakeD()
{
    D = new int[n];
    inverseD = new int[n];

    int** tempInverseD = nullptr;

    for (int i = 0; i < n; ++i)
        D[i] = f->randIntInField();

    //-----------------------------------------------------------------
    //D[0] = 7; D[1] = 3; D[2] = 4; D[3] = 4; D[4] = 7; D[5] = 5; D[6] = 7;

    int** tempD = new int*[n];
    for (int i = 0; i < n; ++i)
    {
        tempD[i] = new int[n];
        tempD[i][i] = D[i];
    }

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (i != j) tempD[i][j] = 0;
        }
    }

    tempInverseD = f->inverseMatrix(tempD, n);

    for (int i = 0; i < n; ++i)
    {
        inverseD[i] = tempInverseD[i][i];
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
    qDebug() << s << endl << "inverseD" << endl << s2;

    for (int i = 0; i < n; ++i)
    {
        delete[] tempD[i];
        delete[] tempInverseD[i];
    }
    delete[] tempD; delete[] tempInverseD;

    return true;
}

bool Koder::MakeH()
{
    H = new int*[t2];
    for (int i = 0; i < t2; ++i)
    {
        H[i] = new int[n];
    }

    for (int i = 0; i < t2; ++i)
        for (int j = 0; j < n; ++j)
        {
            H[i][j] = f->powerNum(Hdots[0][j], Hpoly[i][0]);
            H[i][j] = f->um[H[i][j]][f->powerNum(Hdots[1][j], Hpoly[i][1])];
            H[i][j] = f->um[H[i][j]][f->powerNum(Hdots[2][j], Hpoly[i][2])];
        }

    //-----------------------------------------------------------------
    //	H[0][0] = 1; H[0][1] = 1; H[0][2] = 1; H[0][3] = 1; H[0][4] = 1; H[0][5] = 1; H[0][6] = 1;
    //	H[1][0] = 1; H[1][1] = 2; H[1][2] = 3; H[1][3] = 4; H[1][4] = 5; H[1][5] = 6; H[1][6] = 7;
    //	H[2][0] = 1; H[2][1] = 3; H[2][2] = 5; H[2][3] = 7; H[2][4] = 2; H[2][5] = 4; H[2][6] = 6;
    //	H[3][0] = 1; H[3][1] = 4; H[3][2] = 7; H[3][3] = 3; H[3][4] = 6; H[3][5] = 2; H[3][6] = 5;

    //  H[0][0] = 1; H[0][1] = 2; H[0][2] = 5; H[0][3] = 3; H[0][4] = 2; H[0][5] = 3; H[0][6] = 5;
    //  H[1][0] = 1; H[1][1] = 5; H[1][2] = 7; H[1][3] = 2; H[1][4] = 6; H[1][5] = 4; H[1][6] = 3;
    //  H[2][0] = 1; H[2][1] = 1; H[2][2] = 2; H[2][3] = 1; H[2][4] = 3; H[2][5] = 5; H[2][6] = 1;
    //  H[3][0] = 1; H[3][1] = 4; H[3][2] = 4; H[3][3] = 7; H[3][4] = 7; H[3][5] = 6; H[3][6] = 6;

    qDebug() << endl << "H";
    for (int i = 0; i < t2; ++i)
    {
        QString s = "";
        for (int j = 0; j < n; j++)
        {
            if (H[i][j] >= 0 && H[i][j] <= 9) s += " ";
            s += QString::number(H[i][j]) + " ";
        }
        qDebug() << s;
    }

    return true;
}

bool Koder::MakeX()
{
    do
    {
        X = f->RandMatrixSxS(t2);
        //-----------------------------------------------------------------
        //		X[0][0] = 1; X[0][1] = 7; X[0][2] = 4; X[0][3] = 6;
        //		X[1][0] = 4; X[1][1] = 6; X[1][2] = 3; X[1][3] = 5;
        //		X[2][0] = 0; X[2][1] = 4; X[2][2] = 0; X[2][3] = 5;
        //		X[3][0] = 7; X[3][1] = 4; X[3][2] = 5; X[3][3] = 4;

        inverseX = f->inverseMatrix(X, t2);
        if (inverseX == nullptr)
        {
            for (int i = 0; i < t2; ++i)
            {
                delete[] X[i];
            }
            delete[] X;
        }
    } while (inverseX == nullptr);

    qDebug() << endl << "X                               inverseX";
    for (int i = 0; i < t2; ++i)
    {
        QString s = "";
        QString s2 = "";
        for (int j = 0; j < t2; j++)
        {
            if (X[i][j] >= 0 && X[i][j] <= 9) s += " ";
            //if (X[i][j] >= 10 && X[i][j] <= 99) s += " ";
            s += QString::number(X[i][j]) + " ";
            if (inverseX[i][j] >= 0 && inverseX[i][j] <= 9) s2 += " ";
            //if (inverseX[i][j] >= 10 && inverseX[i][j] <= 99) s2 += " ";
            s2 += QString::number(inverseX[i][j]) + " ";
        }
        qDebug() << s << "|" << s2;
    }

    int** tempX = new int*[t2];
    for (int i = 0; i < t2; ++i)
    {
        tempX[i] = new int[t2];
    }
    tempX = f->MatrixMult(X, t2, t2, inverseX, t2, t2);

    qDebug() << endl << "X cheking";
    for (int i = 0; i < t2; ++i)
    {
        QString s = "";
        for (int j = 0; j < t2; j++)
        {
            if (tempX[i][j] >= 0 && tempX[i][j] <= 9) s += " ";
            s += QString::number(tempX[i][j]) + " ";
        }
        qDebug() << s;
    }

    if (tempX)
    {
        for (int i = 0; i < t2; ++i)
        {
            delete[] tempX[i];
        }
        delete[] tempX;
    }
    return true;
}

bool Koder::MakeHx()
{
    Hx = new int *[t2];
    for (int i = 0; i < t2; ++i)
    {
        Hx[i] = new int[n];
    }

    for (int i = 0; i < t2; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            Hx[i][j] = H[i][j];
        }
    }

    Hx = f->MatrixMult(X, t2, t2, H, t2, n);

    //--------------------------------------------------------
    qDebug() << endl << "XH";
    for (int i = 0; i < t2; ++i)
    {
        QString s = "";
        for (int j = 0; j < n; j++)
        {
            if (Hx[i][j] >= 0 && Hx[i][j] <= 9) s += " ";
            s += QString::number(Hx[i][j]) + " ";
        }
        qDebug() << s;
    }

    //--------------------------------------------------------

    int** tempHx = new int*[t2];
    for (int i = 0; i < t2; ++i)
    {
        tempHx[i] = new int[n];
    }

    for (int i = 0; i < t2; ++i)
        for (int j = 0; j < n; ++j)
        {
            tempHx[i][j] = Hx[i][P[j] - 1];
        }

    for (int i = 0; i < t2; ++i)
    {
        delete[] Hx[i];
    }
    delete[] Hx;

    Hx = tempHx;
    tempHx = nullptr;

    //--------------------------------------------------------
    qDebug() << endl << "XHP";
    for (int i = 0; i < t2; ++i)
    {
        QString s = "";
        for (int j = 0; j < n; j++)
        {
            if (Hx[i][j] >= 0 && Hx[i][j] <= 9) s += " ";
            s += QString::number(Hx[i][j]) + " ";
        }
        qDebug() << s;
    }

    for (int i = 0; i < t2; ++i)
        for (int j = 0; j < n; ++j)
            Hx[i][j] = f->um[Hx[i][j]][D[j]];

    qDebug() << endl << "Hx";
    for (int i = 0; i < t2; ++i)
    {
        QString s = "";
        for (int j = 0; j < n; j++)
        {
            if (Hx[i][j] >= 0 && Hx[i][j] <= 9) s += " ";
            s += QString::number(Hx[i][j]) + " ";
        }
        qDebug() << s;
    }

    return true;
}

void Koder::encode(QString I, int** Hx)
{
    //QStringList splited = I.split(" ");
    //if (splited.size() != t2) qDebug() << splited.size() << " wrong count!";

    int* tempI = new int[n];
    //for (int i = 0; i < t2; ++i)
    //{
    //	tempI[i] = splited[i].toInt();
    //}

    if (I.length() != 3)
    {
        qDebug() << "I wrong length!";
    }
    else
    {
        //QString temps = "";
        //for (int i = 0; i < I.length(); ++i)
        //{
        //	ushort c = I.at(i).toLatin1();
        //	int ch = c;
        //	if (ch >= 64 && ch <= 127) temps += "0";
        //	else if (ch >= 32 && ch <= 63) temps += "00";
        //	else if (ch >= 16 && ch <= 31) temps += "000";
        //	else if (ch >= 8 && ch <= 15) temps += "0000";
        //	else if (ch >= 4 && ch <= 7) temps += "00000";
        //	else if (ch >= 2 && ch <= 3) temps += "000000";
        //	else if (ch >= 0 && ch <= 1) temps += "0000000";
        //	temps += QString::number(ch, 2);
        //}
        //int A = temps.toInt(0, 2);

        nBinEqvVec* v = nBC->getEqvVecByNum(qrand() % 1029);
        for (int i = 0; i < n; ++i) tempI[i] = v->CaInt[i];
        delete v;
        tempI[0] = 3; tempI[1] = 1; tempI[2] = 0; tempI[3] = 0;	tempI[4] = 0; tempI[5] = 0; tempI[6] = 0;

        QString s = "";
        for (int i = 0; i < n; ++i)
            s += QString::number(tempI[i]) + " ";
        s.remove(s.size() - 1, 1);
        qDebug() << endl << "I: " << s;

        //Вместо транспонирования и умножения сразу умножаем в виде строка*строка
        Sx = new int[t2];
        for (int i = 0; i < t2; ++i)
        {
            Sx[i] = 0;
            for (int j = 0; j < n; ++j)
                Sx[i] = f->pl[Sx[i]][f->um[tempI[j]][Hx[i][j]]];
        }

        s = "";
        for (int i = 0; i < t2; ++i)
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
        //if (i >= t2) Cx[i] = 1;
        if (i < (n - t2)) Cx[i] = 1;
        else Cx[i] = 0;
    }
    //Cx[0] = 0; Cx[1] = 1; Cx[2] = 0; Cx[3] = 0; Cx[4] = 0; Cx[5] = 0; Cx[6] = 0;

    int** tempArr = new int*[t2];
    for (int i = 0; i < t2; ++i)
    {
        tempArr[i] = new int[t2 + 1];
    }

    for (int i = 0; i < t2; ++i)
    {
        //for (int j = 0; j < t2; ++j)
        for (int j = n - t2; j < n; ++j)
        {
            //tempArr[i][j] = Hx[i][j];
            tempArr[i][j - (n - t2)] = Hx[i][j];
        }
    }

    int* tempSx = new int[t2];
    //for (int i = 0; i < t2; ++i)
    //{
    //	tempSx[i] = 0;
    //	for (int j = 0; j < t2; ++j)
    //		tempSx[i] = f->pl[tempSx[i]][f->um[Sx[j]][inverseX[j][i]]];
    //}

    //for (int i = 0; i < t2; ++i)
    //{
    //	Sx[i] = tempSx[i];
    //}
    //delete[] tempSx;
    //tempSx = nullptr;

    QString s = "";
    //for (int i = 0; i < t2; ++i)
    //	s += QString::number(Sx[i]) + " ";
    //s.remove(s.size() - 1, 1);

    //qDebug() << endl << "Sx: " << s;

    for (int i = 0; i < t2; ++i)
    {
        tempArr[i][t2] = Sx[i];
        //for (int j = t2; j < n; ++j)
        for (int j = 0; j < (n - t2); ++j)
        {
            tempArr[i][t2] = f->pl[tempArr[i][t2]][f->um[Hx[i][j]][Cx[j]]];
        }
    }

    qDebug() << endl << "tempArr";
    for (int i = 0; i < t2; ++i)
    {
        QString s = "";
        for (int j = 0; j < t2 + 1; j++)
        {
            if (tempArr[i][j] >= 0 && tempArr[i][j] <= 9) s += " ";
            s += QString::number(tempArr[i][j]) + " ";
        }
        qDebug() << s;
    }

    tempSx = gauss(tempArr, t2, t2 + 1);

    //for (int i = 0; i < t2; ++i)
    for (int i = t2 - 1; i < n; ++i)
    {
        //Cx[i] = tempSx[i];
        Cx[i] = tempSx[i - (n - t2)];
    }
    delete[] tempSx;


    s = "";
    for (int i = 0; i < n; ++i)
        s += QString::number(Cx[i]) + " ";
    s.remove(s.size() - 1, 1);
    qDebug() << endl << "Cx: " << s;


    int* tempCx = new int[n];

    for (int j = 0; j < n; ++j)
        tempCx[j] = Cx[j];


    for (int j = 0; j < n; ++j)
        tempCx[j] = f->um[Cx[j]][inverseD[j]];

    s = "";
    for (int i = 0; i < n; ++i)
        s += QString::number(tempCx[i]) + " ";
    s.remove(s.size() - 1, 1);
    qDebug() << endl << "Cx * D^-1: " << s;

    for (int j = 0; j < n; ++j)
        Cx[j] = tempCx[inverseP[j] - 1];

    s = "";
    for (int i = 0; i < n; ++i)
        s += QString::number(Cx[i]) + " ";
    s.remove(s.size() - 1, 1);
    qDebug() << endl << "Cx`: " << s;

    tempSx = new int[t2];
    for (int i = 0; i < t2; ++i)
    {
        tempSx[i] = 0;
        for (int j = 0; j < n; ++j)
            tempSx[i] = f->pl[tempSx[i]][f->um[Cx[j]][H[i][j]]];
    }

    s = "";
    for (int i = 0; i < t2; ++i)
        s += QString::number(tempSx[i]) + " ";
    s.remove(s.size() - 1, 1);
    qDebug() << endl << "Sx`: " << s;

    for (int i = 0; i < t2; ++i)
    {
        delete[] tempArr[i];
    }
    tempArr = new int*[t];
    for (int i = 0; i < t; ++i)
    {
        tempArr[i] = new int[t + 1];
    }

    for (int i = 0; i < t; ++i)
        for (int j = 0; j < t + 1; ++j)
            tempArr[i][j] = tempSx[i + j];

    qDebug() << endl << "tempArr";
    for (int i = 0; i < t; ++i)
    {
        QString s = "";
        for (int j = 0; j < t + 1; j++)
        {
            if (tempArr[i][j] >= 0 && tempArr[i][j] <= 9) s += " ";
            s += QString::number(tempArr[i][j]) + " ";
        }
        qDebug() << s;
    }

    int* alpha = gauss(tempArr, t, t + 1);
    s = "";
    for (int i = 0; i < t; ++i)
        s += QString::number(alpha[i]) + " ";
    s.remove(s.size() - 1, 1);
    qDebug() << endl << "alpha: " << s;

    int* tempE = new int[n];
    for (int i = 0; i < n; ++i)
    {
        //tempE[i] = 0;
        //tempE[i] = f->pl[tempE[i]][f->powerNum(i + 1, t)];
        //for (int j = (t - 1); j > 0; --j)
        //{
        //	tempE[i] = f->pl[tempE[i]][f->um[alpha[j]][f->powerNum(i + 1, j)]];
        //}
        //tempE[i] = f->pl[tempE[i]][alpha[0]];
        tempE[i] = f->powerNum(i + 1, t);
        for (int j = 0; j < (t - 1); ++j)
        {
            tempE[i] = f->pl[tempE[i]][f->um[alpha[j]][f->powerNum(i + 1, t - j - 1)]];
        }
        tempE[i] = f->pl[tempE[i]][alpha[t - 1]];
    }

    s = "";
    for (int i = 0; i < n; ++i)
        s += QString::number(tempE[i]) + " ";
    s.remove(s.size() - 1, 1);
    qDebug() << endl << "e`: " << s;

    for (int i = 0; i < n; ++i)
    {
        if (!tempE[i]) tempE[i] = 1;
        else tempE[i] = 0;
    }
    s = "";
    for (int i = 0; i < n; ++i)
        s += QString::number(tempE[i]) + " ";
    s.remove(s.size() - 1, 1);
    qDebug() << endl << "e: " << s;


    int** foldArr = new int*[t2];
    for (int i = 0; i < t2; ++i)
    {
        foldArr[i] = new int[t + 1];
    }
    for (int i = 0; i < t2; ++i)
    {
        int fl = 0;
        for (int j = 0; j < n; ++j)
        {
            if (tempE[j])
            {
                foldArr[i][fl] = H[i][j];
                ++fl;
            }
        }
        foldArr[i][fl] = tempSx[i];
    }

    qDebug() << endl << "fold";
    for (int i = 0; i < t2; ++i)
    {
        QString s = "";
        for (int j = 0; j < t + 1; j++)
        {
            if (foldArr[i][j] >= 0 && foldArr[i][j] <= 9) s += " ";
            s += QString::number(foldArr[i][j]) + " ";
        }
        qDebug() << s;
    }


    for (int ii = 0; ii < 3; ++ii)
    {
        for (int i = 0; i < t; ++i)
            for (int j = 0; j < t + 1; ++j)
                tempArr[i][j] = foldArr[ii + i][j];
        qDebug() << endl << "tempArr";
        for (int i = 0; i < t; ++i)
        {
            QString s = "";
            for (int j = 0; j < t + 1; j++)
            {
                if (tempArr[i][j] >= 0 && tempArr[i][j] <= 9) s += " ";
                s += QString::number(tempArr[i][j]) + " ";
            }
            qDebug() << s;
        }
        int* alpha2 = gauss(tempArr, t, t + 1);
        s = "";
        for (int i = 0; i < t; ++i)
            s += QString::number(alpha2[i]) + " ";
        s.remove(s.size() - 1, 1);
        qDebug() << "alpha: " << s;

        int fl = 0;
        for (int j = 0; j < n; ++j)
        {
            if (tempE[j])
            {
                tempCx[j] = alpha2[fl];
                ++fl;
            }
            else tempCx[j] = 0;
        }
        s = "";
        for (int i = 0; i < n; ++i)
            s += QString::number(tempCx[i]) + " ";
        s.remove(s.size() - 1, 1);
        qDebug() << endl << "e`: " << s;

        delete[] alpha2;
    }

    //for (int j = 0; j < n; ++j)
    //	Cx[j] = tempE[P[j] - 1];
    //s = "";
    //for (int i = 0; i < n; ++i)
    //	s += QString::number(Cx[i]) + " ";
    //s.remove(s.size() - 1, 1);
    //qDebug() << endl << "e * P: " << s;

    //for (int j = 0; j < n; ++j)
    //	tempE[j] = f->um[Cx[j]][D[j]];
    //s = "";
    //for (int i = 0; i < n; ++i)
    //	s += QString::number(tempE[i]) + " ";
    //s.remove(s.size() - 1, 1);
    //qDebug() << endl << "e * P * D: " << s;

    //for (int i = 0; i < t2; ++i)
    //{
    //	tempSx[i] = 0;
    //	for (int j = 0; j < n; ++j)
    //		tempSx[i] = f->pl[tempSx[i]][f->um[tempE[j]][H[i][j]]];
    //}
    //s = "";
    //for (int i = 0; i < t2; ++i)
    //	s += QString::number(tempSx[i]) + " ";
    //s.remove(s.size() - 1, 1);
    //qDebug() << endl << "Sx: " << s;

    //for (int i = 0; i < n; ++i)
    //{
    //	tempE[i] = 0;
    //	tempE[i] = f->powerNum(i + 1, t);
    //	for (int j = 0; j < (t - 1); ++j)
    //	{
    //		tempE[i] = f->pl[tempE[i]][f->um[alpha2[j]][f->powerNum(i + 1, t - j - 1)]];
    //	}
    //	tempE[i] = f->pl[tempE[i]][alpha2[t - 1]];
    //}
    //for (int i = 0; i < n; ++i)
    //	s += QString::number(tempE[i]) + " ";
    //s.remove(s.size() - 1, 1);
    //qDebug() << "?: " << s;

    for (int i = 0; i < t; ++i)
    {
        delete[] tempArr[i];
    }
    delete[] tempArr;

    for (int i = 0; i < t2; ++i)
    {
        delete[] foldArr[i];
    }
    delete[] foldArr;

    delete[] tempSx;
    delete[] tempCx;
    delete[] alpha;
    delete[] tempE;
    delete[] Cx;

    return " ";
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
        int r = i; //выбор наибольшего элемента в столбце
        for (int j = i; j < row; j++)
        {
            if (resultPH[j][i] > resultPH[r][i]) r = j;
        }
        if (r != i) //замена строк если это необходимо
        {
            for (int j = i; j < col; j++)
            {
                int t = resultPH[i][j];
                resultPH[i][j] = resultPH[r][j];
                resultPH[r][j] = t;
            }
        }

        for (int j = i + 1; j < row; j++) //вычитание текущей строки из всех нижних, умножив на соответствующий коэффициент
        {
            int c = f->del(resultPH[j][i], resultPH[i][i]);
            //for (r = i; r < col; r++)
            for (r = i; r < col; r++)
            {
                resultPH[j][r] = f->pl[resultPH[j][r]][f->um[resultPH[i][r]][c]];
            }
        }
    }

    qDebug() << endl << "gauss...";
    for (int i = 0; i < row; ++i)
    {
        QString s = "";
        for (int j = 0; j < (row + 1); j++)
        {
            //if (resultPH[i][j] >= 0 && resultPH[i][j] <= 9) s += " ";
            s += QString::number(resultPH[i][j]) + " ";
        }
        qDebug() << s;
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

void Koder::test_encode_nonBEC()
{
    qDebug() <<endl<< "--------------------------Testing encode NBEC----------------------------------";

    QTime time;
    //
    nBC = new nBinEqvCod(n, t, (n + 1));

    int n_bits = nBC->getM().get_str(2).capacity();
    gmp_randclass s(gmp_randinit_default);
    s.seed(time.msecsSinceStartOfDay());

    for (int i = 0; i < 10; ++i)
    {
        mpz_class temp = s.get_z_bits(n_bits);
        if(temp>nBC->getM()) temp = s.get_z_bits(n_bits-1);

        time = QTime::currentTime();
        nBinEqvVec *v1 = nBC->getEqvVecByNum(temp);
        qDebug() << endl << v1->ToStr();

        nBinEqvVec *v2 = nBC->getEqvVecByCa(v1->Ca);
        qDebug() << v2->ToStr();

        qDebug() << ((v1->A == v2->A) ? "OK" : "---------------\n>>>>>ERROR<<<<<\n---------------");
        qDebug() << "Time spend: " << time.msecsTo(QTime::currentTime()) << " msecs";
        delete v1;
        delete v2;
    }
    if (nBC) delete nBC; nBC = nullptr;
}

#include <chrono>
void Koder::mult_test(int size, bool matrixPrint)
{
    qDebug() << endl << "---------------Matrix_multiple_testing--------------";
    int** C;

    srand(time(NULL));

    int **A = new int *[size];
    int **B = new int *[size];
    for(int i=0; i<size; ++i)
    {
        A[i] = new int[size];
        B[i] = new int[size];
        for(int j=0; j<size; ++j)
        {
            A[i][j] = 0+rand()%15;
            B[i][j] = 0+rand()%15;
        }
    }

    if(matrixPrint)
    {
        qDebug() << endl << "A";
        for (int i = 0; i < size; ++i)
        {
            QString s = "";
            for (int j = 0; j < size; j++)
            {
                s += QString::number(A[i][j]) + " ";
            }
            qDebug() << s;
        }
        qDebug() << endl << "B";
        for (int i = 0; i < size; ++i)
        {
            QString s = "";
            for (int j = 0; j < size; j++)
            {
                s += QString::number(B[i][j]) + " ";
            }
            qDebug() << s;
        }
    }

    qDebug() << endl << "-----------------Local_mult-------------------------";
    auto start = std::chrono::system_clock::now();
    C = f->MatrixMult(A, size, size, B, size, size);
    auto end = std::chrono::system_clock::now();
    if(matrixPrint)
    {
        qDebug() << "C";
        for (int i = 0; i < size; ++i)
        {
            QString s = "";
            for (int j = 0; j < size; j++)
            {
                s += QString::number(C[i][j]) + " ";
            }
            qDebug() << s;
        }
    }
    std::chrono::duration<double> diff = end - start;
    qDebug() << "Time spend: " << diff.count() << endl;
    for (int i = 0; i < size; ++i)
    {
        delete[] C[i];
    }
    delete[] C;

    qDebug() << endl << "-----------------OpenMP_mult------------------------";
    start = std::chrono::system_clock::now();
    C = f->MatrixMultOpenMP(A, size, size, B, size, size);
    end = std::chrono::system_clock::now();
    if(matrixPrint)
    {
        qDebug() << "C";
        for (int i = 0; i < size; ++i)
        {
            QString s = "";
            for (int j = 0; j < size; j++)
            {
                s += QString::number(C[i][j]) + " ";
            }
            qDebug() << s;
        }
    }
    diff = end - start;
    qDebug() << "Time spend: " << diff.count() << endl;
    for (int i = 0; i < size; ++i)
    {
        delete[] C[i];
    }
    delete[] C;

    qDebug() << endl << "------------------Fast_mult-------------------------";
    start = std::chrono::system_clock::now();
    C = f->MatrixMultFast(A, size, size, B, size, size);
    end = std::chrono::system_clock::now();
    if(matrixPrint)
    {
        qDebug() << "C";
        for (int i = 0; i < size; ++i)
        {
            QString s = "";
            for (int j = 0; j < size; j++)
            {
                s += QString::number(C[i][j]) + " ";
            }
            qDebug() << s;
        }
    }
    diff = end - start;
    qDebug() << "Time spend: " << diff.count() << endl;
    for (int i = 0; i < size; ++i)
    {
        delete[] C[i];
    }
    delete[] C;
}

void Koder::test(int n, int w, int q)
{
    qDebug() << endl << "-------------------Testing--------------------------";

    nBC = new nBinEqvCod(n, w, q);
    nBC->test();
    for (int i = 0; i < nBC->getM(); ++i)
        qDebug() << nBC->getStringEqvVec(i);
    if (nBC) delete nBC; nBC = nullptr;
}

void Koder::Fact_test()
{
    qDebug() << endl << "-------------------------MPZ_Testing---------------------------------";
    mpz_class b("15511187532873822802000000000301646930321106000000000000020016986112000000000000"),
            c("60167034436371883082406479598617053710361955327150763910740921313000000000375665858281472000000000000");
    qDebug() << "b: " << b.get_str().c_str() << endl << "c: " << c.get_str().c_str();
    c *= b;
    qDebug() << "c *= b res: " << c.get_str().c_str();
    qDebug() << endl << "-----------------------FactTreeTesting-------------------------------";
    int n = 1000;
    mpz_class a = FactTree(n);// FactTree(100);
    qDebug() << "FactTree("<< n <<"): " << a.get_str().c_str();
}
