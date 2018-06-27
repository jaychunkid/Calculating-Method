/*�����Դ����������㷨*/
#ifndef SYSTEM_OF_LINEAR_ALGEBRAIC_EQUATIONS_H
#define SYSTEM_OF_LINEAR_ALGEBRAIC_EQUATIONS_H

//˳���˹��ȥ��
bool Gauss1(double** A, double* b, double e, unsigned n, double* x);

//����Ԫ��˹��ȥ��
bool Gauss2(double** A, double* B, double e, unsigned n, double* x);

//ȫ��Ԫ��˹��ȥ��
bool Gauss3(double** A, double* b, double e, unsigned n, double* x);

//ֱ��LU�ֽⷨ
bool LU(double** A, double* b, double e, unsigned n, double* x);

//Jacobi����
int Jacobi(double** A, double* b, double* y, double e, unsigned n, 
unsigned M, double* x);

//Seidel����
int Seidel(double** A, double* b, double* y, double e, unsigned n, 
unsigned M, double* x);

//�ɳڵ���
int SOR(double** A, double* b, double* y, double w, double e, unsigned n, 
unsigned M, double* x);

#endif