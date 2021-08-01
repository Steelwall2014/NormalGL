#pragma once
#include "Utils.h"
class UserCamera;
class Matrix
{
public:
	Matrix() 
	{ 
		shape[0] = 0; shape[1] = 0;
		data = nullptr;
	}
	Matrix(const Matrix &);
	Matrix(int row, int col, double value);
	Matrix(int row, int col);
	Matrix(int col, ...);
	~Matrix();
	int shape[2];
	double **data; 
	void operator=(const Matrix &other);
	Matrix operator+(const Matrix &other);
	Matrix operator-();
	Matrix operator-(const Matrix &other);
	void operator+=(const Matrix &other);
	void operator-=(const Matrix &other);
	Matrix operator*(const Matrix &other);
	Matrix operator*(double other);
	void operator*=(const Matrix &other);
	Matrix GetTransposedInvMatrix();			// �������ת��
	Matrix GetInvMatrix();						// �������
	void Transpose();							// ת��
	Matrix CrossProduct(const Matrix &other);	// ���
	double DotProduct(const Matrix &other);		// ���
	void Normalize();							// ��һ��
	static Matrix CreateTranslationMatrix(double offset_x, double offset_y);	// ƽ�ƾ���
	static Matrix CreateRotateMatrix(double theta);								// ��ת����
	static Matrix CreateZoomMatrix(double);										// ����
	static Matrix CreateZoomWithPointMatrix(double, double, double);			// �Ƶ�����
	static Matrix CreateIdentityMatrix(int);									// ��λ��
	static Matrix CreateRotateWithCenterMatrix(double theta, int x, int y);		// �Ƶ���ת
	static Matrix CreateViewMatrix(UserCamera &userCam);						// �۲�任����
	static Matrix CreatePerspectiveMatrix(UserCamera &uc);						// ͶӰ�任����
	static Matrix CreateRotateWithAxisMatrix_3D(Matrix &axis, Point3D point, double theta);	// ��ά������ת
	static Matrix CreateTranslationMatrix_3D(double offset_x, double offset_y, double offset_z);// ��άƽ��
};