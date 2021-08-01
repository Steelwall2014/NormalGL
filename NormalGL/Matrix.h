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
	Matrix GetTransposedInvMatrix();			// 求逆矩阵并转置
	Matrix GetInvMatrix();						// 求逆矩阵
	void Transpose();							// 转置
	Matrix CrossProduct(const Matrix &other);	// 叉乘
	double DotProduct(const Matrix &other);		// 点乘
	void Normalize();							// 归一化
	static Matrix CreateTranslationMatrix(double offset_x, double offset_y);	// 平移矩阵
	static Matrix CreateRotateMatrix(double theta);								// 旋转矩阵
	static Matrix CreateZoomMatrix(double);										// 缩放
	static Matrix CreateZoomWithPointMatrix(double, double, double);			// 绕点缩放
	static Matrix CreateIdentityMatrix(int);									// 单位阵
	static Matrix CreateRotateWithCenterMatrix(double theta, int x, int y);		// 绕点旋转
	static Matrix CreateViewMatrix(UserCamera &userCam);						// 观察变换矩阵
	static Matrix CreatePerspectiveMatrix(UserCamera &uc);						// 投影变换矩阵
	static Matrix CreateRotateWithAxisMatrix_3D(Matrix &axis, Point3D point, double theta);	// 三维绕轴旋转
	static Matrix CreateTranslationMatrix_3D(double offset_x, double offset_y, double offset_z);// 三维平移
};