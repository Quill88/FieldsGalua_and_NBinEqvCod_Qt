#pragma once
#ifndef NBINEQVVEC_H
#define NBINEQVVEC_H

#include <QString>

/*non�binary equivalent vector
���������� ����������� ������*/
class nBinEqvVec
{
public:
	nBinEqvVec();
	~nBinEqvVec();

    QString ToStr();
	
	int A;
	int n;
	int w;

	int* ab;
	int* a;

    QString Ca;
};

#endif /* NBINEQVVEC_H */
