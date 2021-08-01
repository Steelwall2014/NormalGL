#include "Matrix.h"
#include "Utils.h"
#include "UserCamera.h"

#include <cmath>

Matrix::Matrix(int row, int col)
{
	shape[0] = row; shape[1] = col;
	data = new double *[row];
	for (int i = 0; i < row; i++)
	{
		data[i] = new double[col];
		for (int j = 0; j < col; j++)
		{
			data[i][j] = 0;
		}
	}
}
Matrix::Matrix(int row, int col, double value)
{
	shape[0] = row; shape[1] = col;
	data = new double *[row];
	for (int i = 0; i < row; i++)
	{
		data[i] = new double[col];
		for (int j = 0; j < col; j++)
		{
			data[i][j] = value;
		}
	}
}
Matrix::Matrix(const Matrix &other)
{
	shape[0] = other.shape[0]; shape[1] = other.shape[1];
	data = new double *[shape[0]];
	for (int i = 0; i < shape[0]; i++)
	{
		data[i] = new double[shape[1]];
		for (int j = 0; j < shape[1]; j++)
		{
			data[i][j] = other.data[i][j];
		}
	}
}
Matrix::Matrix(int col, ...)
{
	shape[0] = 1; shape[1] = col;
	data = new double *[1];
	data[0] = new double[col];
	va_list arg_ptr;
	va_start(arg_ptr, col);
	for (int i = 0; i < col; ++i)
	{
		
		data[0][i] = va_arg(arg_ptr, double);
	}
	va_end(arg_ptr);
}
Matrix::~Matrix()
{
	for (int i = 0; i < shape[0]; i++)
		delete[] data[i];
	delete[] data;
}
void Matrix::operator=(const Matrix &other)
{
	//Matrix *new_matrix = new Matrix(other.shape[0], other.shape[1]);
	if (shape[0] != other.shape[0] || shape[1] != other.shape[1])
	{
		for (int i = 0; i < shape[0]; i++)
			delete[] data[i];
		delete[] data;
		data = new double *[other.shape[0]];
		for (int i = 0; i < other.shape[0]; i++)
			data[i] = new double[other.shape[1]];

		shape[0] = other.shape[0];
		shape[1] = other.shape[1];
	}
	for (int i = 0; i < other.shape[0]; i++)
	{
		for (int j = 0; j < other.shape[1]; j++)
		{
			this->data[i][j] = other.data[i][j];
		}
	}
	//return *this;
}
Matrix Matrix::operator+(const Matrix &other)
{
	if (this->shape[0] != other.shape[0] || this->shape[1] != other.shape[1])
		throw("invalid shapes");
	Matrix new_matrix = Matrix(other.shape[0], other.shape[1]);
	for (int i = 0; i < other.shape[0]; i++)
	{
		for (int j = 0; j < other.shape[1]; j++)
		{
			new_matrix.data[i][j] = this->data[i][j] + other.data[i][j];
		}
	}
	return new_matrix;
}
Matrix Matrix::operator-()
{
	Matrix new_matrix = Matrix(this->shape[0], this->shape[1]);
	for (int i = 0; i < this->shape[0]; i++)
	{
		for (int j = 0; j < this->shape[1]; j++)
		{
			new_matrix.data[i][j] = -this->data[i][j];
		}
	}
	return new_matrix;
}
void Matrix::operator+=(const Matrix &other)
{
	if (this->shape[0] != other.shape[0] || this->shape[1] != other.shape[1])
		throw("invalid shapes");
	for (int i = 0; i < other.shape[0]; i++)
	{
		for (int j = 0; j < other.shape[1]; j++)
		{
			this->data[i][j] += other.data[i][j];
		}
	}
	//return *this;
}
Matrix Matrix::operator-(const Matrix &other)
{
	if (this->shape[0] != other.shape[0] || this->shape[1] != other.shape[1])
		throw("invalid shapes");
	Matrix new_matrix = Matrix(other.shape[0], other.shape[1]);
	for (int i = 0; i < other.shape[0]; i++)
	{
		for (int j = 0; j < other.shape[1]; j++)
		{
			new_matrix.data[i][j] = this->data[i][j] - other.data[i][j];
		}
	}
	return new_matrix;
}
void Matrix::operator-=(const Matrix &other)
{
	if (this->shape[0] != other.shape[0] || this->shape[1] != other.shape[1])
		throw("invalid shapes");
	for (int i = 0; i < other.shape[0]; i++)
	{
		for (int j = 0; j < other.shape[1]; j++)
		{
			this->data[i][j] -= other.data[i][j];
		}
	}
	//return *this;
}
Matrix Matrix::operator*(const Matrix &other)
{
	if (this->shape[1] != other.shape[0])
		throw("invalid shapes");
	Matrix new_matrix = Matrix(this->shape[0], other.shape[1], 0);
	for (int i = 0; i < this->shape[0]; i++)
		for (int j = 0; j < other.shape[1]; j++)
			for (int k = 0; k < this->shape[1]; k++)
				new_matrix.data[i][j] += this->data[i][k] * other.data[k][j];
	return new_matrix;
}
void Matrix::operator*=(const Matrix &other)
{
	if (this->shape[1] != other.shape[0])
		throw("invalid shapes");
	Matrix new_matrix(this->shape[0], other.shape[1], 0);
	for (int i = 0; i < this->shape[0]; i++)
		for (int j = 0; j < other.shape[1]; j++)
			for (int k = 0; k < this->shape[1]; k++)
				new_matrix.data[i][j] += this->data[i][k] * other.data[k][j];
	for (int i = 0; i < this->shape[0]; i++)
		for (int j = 0; j < other.shape[1]; j++)
			this->data[i][j] = new_matrix.data[i][j];
	//return *this;
}
Matrix Matrix::operator*(double other)
{
	Matrix new_matrix = Matrix(this->shape[0], this->shape[1], 0);
	for (int i = 0; i < this->shape[0]; i++)
		for (int j = 0; j < this->shape[1]; j++)
			new_matrix.data[i][j] = this->data[i][j] * other;
	return new_matrix;
}
Matrix Matrix::CrossProduct(const Matrix &other)
{
	Matrix m = Matrix(1, 3);
	m.data[0][0] = this->data[0][1] * other.data[0][2] - other.data[0][1] * this->data[0][2];
	m.data[0][1] = other.data[0][0] * this->data[0][2] - this->data[0][0] * other.data[0][2];
	m.data[0][2] = this->data[0][0] * other.data[0][1] - other.data[0][0] * this->data[0][1];
	return m;
}
double Matrix::DotProduct(const Matrix &other)
{
	double dpResult = 0;
	dpResult = this->data[0][0] * other.data[0][0] + this->data[0][1] * other.data[0][1] + this->data[0][2] * other.data[0][2];
	return dpResult;
}
void Matrix::Normalize()
{
	for (int i = 0; i < this->shape[0]; i++)
	{
		double len = 0;
		for (int j = 0; j < this->shape[1]; j++)
		{
			len += this->data[i][j] * this->data[i][j];
		}
		len = sqrt(len);
		for (int j = 0; j < this->shape[1]; j++)
		{
			this->data[i][j] /= len;
		}
	}
}
Matrix Matrix::CreateTranslationMatrix(double offset_x, double offset_y)
{
	Matrix matrix = CreateIdentityMatrix(3);
	matrix.data[2][0] = offset_x;
	matrix.data[2][1] = offset_y;
	return matrix;
}
Matrix Matrix::CreateRotateMatrix(double theta)
{
	Matrix matrix = CreateIdentityMatrix(3);
	matrix.data[0][0] = cos(radians(theta));
	matrix.data[0][1] = sin(radians(theta));
	matrix.data[1][0] = -sin(radians(theta));
	matrix.data[1][1] = cos(radians(theta));
	return matrix;
}
Matrix Matrix::CreateRotateWithCenterMatrix(double theta, int center_x, int center_y)
{
	Matrix translation1 = CreateTranslationMatrix(-center_x, -center_y);
	Matrix rotate = CreateRotateMatrix(theta);
	Matrix translation2 = CreateTranslationMatrix(center_x, center_y);
	return translation1 * rotate * translation2;
}
Matrix Matrix::CreateZoomMatrix(double zoom_scale)
{
	Matrix matrix = CreateIdentityMatrix(3);
	matrix.data[0][0] *= zoom_scale;
	matrix.data[1][1] *= zoom_scale;
	return matrix;
}
Matrix Matrix::CreateZoomWithPointMatrix(double zoom_scale, double center_x, double center_y)
{
	Matrix translation1 = CreateTranslationMatrix(-center_x, -center_y);
	Matrix zoom = CreateZoomMatrix(zoom_scale);
	Matrix translation2 = CreateTranslationMatrix(center_x, center_y);
	return translation1 * zoom * translation2;
}
Matrix Matrix::CreateIdentityMatrix(int side_length)
{
	Matrix matrix = Matrix(side_length, side_length);
	for (int i = 0; i < side_length; i++)
	{
		for (int j = 0; j < side_length; j++)
		{
			if (i == j)
				matrix.data[i][j] = 1;
			else
				matrix.data[i][j] = 0;
		}
	}
	return matrix;
}

Matrix Matrix::CreateViewMatrix(UserCamera &userCam)
{
	Matrix viewMatrix = Matrix(4, 4, 0);
	Matrix p = Matrix(3, -userCam.eye_x, -userCam.eye_y, -userCam.eye_z);//p为负视线向量

	//第一列
	viewMatrix.data[0][0] = userCam.u.data[0][0];
	viewMatrix.data[1][0] = userCam.u.data[0][1];
	viewMatrix.data[2][0] = userCam.u.data[0][2];
	//第二列
	viewMatrix.data[0][1] = userCam.v.data[0][0];
	viewMatrix.data[1][1] = userCam.v.data[0][1];
	viewMatrix.data[2][1] = userCam.v.data[0][2];
	//第三列
	viewMatrix.data[0][2] = userCam.n.data[0][0];
	viewMatrix.data[1][2] = userCam.n.data[0][1];
	viewMatrix.data[2][2] = userCam.n.data[0][2];
	//第四行
	viewMatrix.data[3][0] = p.DotProduct(userCam.u);
	viewMatrix.data[3][1] = p.DotProduct(userCam.v);
	viewMatrix.data[3][2] = p.DotProduct(userCam.n);

	viewMatrix.data[3][3] = 1;

	return viewMatrix;
}

Matrix Matrix::CreatePerspectiveMatrix(UserCamera &uc)
{
	Matrix perspectiveMatrix = Matrix(4, 4, 0);

	double dz_el = uc.z_near - uc.eye_z;//视点与视平面的z差值
	double reciprocal_tan_halffov = 1 / (double)tan(uc.fovy / 2);//tan的倒数

	perspectiveMatrix.data[0][0] = (double)reciprocal_tan_halffov / uc.aspect;
	perspectiveMatrix.data[1][1] = reciprocal_tan_halffov;

	double dz_nf = uc.z_near - uc.z_far;//近平面和远平面的z差,ppt中矩阵的n-f

	perspectiveMatrix.data[2][2] = uc.z_far / dz_nf;
	perspectiveMatrix.data[3][2] = uc.z_far * uc.z_near / dz_nf;

	perspectiveMatrix.data[2][3] = -1;

	return perspectiveMatrix;
}

Matrix Matrix::CreateRotateWithAxisMatrix_3D(Matrix &axis, Point3D point, double theta)
{
	Matrix temp_axis = axis;
	temp_axis.Normalize();
	Matrix rotation = CreateIdentityMatrix(4);
	double u = temp_axis.data[0][0], v = temp_axis.data[0][1], w = temp_axis.data[0][2];
	double Cos = cos(radians(theta)), Sin = sin(radians(theta));
	double a = point.x, b = point.y, c = point.z;

	rotation.data[0][0] = u * u + (v * v + w * w) * Cos;
	rotation.data[1][0] = u * v * (1 - Cos) - w * Sin;
	rotation.data[2][0] = u * w * (1 - Cos) + v * Sin;
	rotation.data[3][0] = (a * (v * v + w * w) - u * (b * v + c * w)) * (1 - Cos) + (b * w - c * v) * Sin;
	rotation.data[0][1] = u * v * (1 - Cos) + w * Sin;
	rotation.data[1][1] = v * v + (u * u + w * w) * Cos;
	rotation.data[2][1] = v * w * (1 - Cos) - u * Sin;
	rotation.data[3][1] = (b * (u * u + w * w) - v * (a * u + c * w)) * (1 - Cos) + (c * u - a * w) * Sin;
	rotation.data[0][2] = u * w * (1 - Cos) - v * Sin;
	rotation.data[1][2] = v * w * (1 - Cos) + u * Sin;
	rotation.data[2][2] = w * w + (u * u + v * v) * Cos;
	rotation.data[3][2] = (c * (u * u + v * v) - w * (a * u + b * v)) * (1 - Cos) + (a * v - b * u) * Sin;
	return rotation;
}

Matrix Matrix::CreateTranslationMatrix_3D(double offset_x, double offset_y, double offset_z)
{
	Matrix matrix = CreateIdentityMatrix(4);
	matrix.data[3][0] = offset_x;
	matrix.data[3][1] = offset_y;
	matrix.data[3][2] = offset_z;
	return matrix;
}

void Matrix::Transpose()
{
	int row = this->shape[0];
	int col = this->shape[1];
	int i = 0;

	for (i; i < row - 1; ++i)
	{
		for (int j = i + 1; j < col; ++j)
		{
			swap(this->data[i][j], this->data[j][i]);
		}
	}
}

Matrix Matrix::GetTransposedInvMatrix()
{
	return LUP_solve_inverse(*this);
}

Matrix Matrix::GetInvMatrix()
{
	Matrix trsd = this->GetTransposedInvMatrix();
	trsd.Transpose();
	return trsd;
}