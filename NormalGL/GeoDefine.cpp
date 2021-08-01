#include "GeoDefine.h"
#include "Matrix.h"
#include "Tin.h"
Point3D Point3D::operator*(Matrix &trans)
{
	Point3D result;
	double vector[4], result_vector[4];
	vector[0] = x; vector[1] = y; vector[2] = z; vector[3] = 1;
	for (int i = 0; i < 4; i++)
	{
		double value = 0;
		for (int j = 0; j < 4; j++)
			value += vector[j] * trans.data[j][i];
		result_vector[i] = value;
	}
	result.id = this->id;
	result.x = result_vector[0]; // / result_vector[3];
	result.y = result_vector[1]; // / result_vector[3];
	result.z = result_vector[2]; // / result_vector[3];
	result.depth = result_vector[3];
	return result;
}

void Point3D::Normalize()
{
	double length = sqrt(x * x + y * y + z * z);
	x /= length;
	y /= length;
	z /= length;
}

double Point3D::DotProduct(Matrix &other)
{
	double result = 0;
	result += this->x * other.data[0][0];
	result += this->y * other.data[0][1];
	result += this->z * other.data[0][2];
	return result;
}

//Point3D Face::operator[](int index)
//{
//	return (*this->Tin_grid)[point_index[index]];
//}
//Face::Face(Tin *Tin_grid, int pt1_index, int pt2_index, int pt3_index)
//{
//	this->point_index[0] = pt1_index;
//	this->point_index[1] = pt2_index;
//	this->point_index[2] = pt3_index;
//	this->Tin_grid = Tin_grid;
//	//DetermineTriShape();
//}
//
//Face::Face(Tin *Tin_grid, int pts_index[3])
//{
//	this->point_index[0] = pts_index[0];
//	this->point_index[1] = pts_index[1];
//	this->point_index[2] = pts_index[2];
//	this->Tin_grid = Tin_grid;
//	//DetermineTriShape();
//}
//
//Face::Face(const Face &other)
//{
//	this->point_index[0] = other.point_index[0];
//	this->point_index[1] = other.point_index[1];
//	this->point_index[2] = other.point_index[2];
//	this->Tin_grid = other.Tin_grid;
//	//DetermineTriShape();
//}

/// <summary>
/// 判断光栅化三角形形状以及顶点状态
/// </summary>
//void Face::DetermineTriShape()
//{
//	if ((*Tin_grid)[this->point_index[0]].y == (*Tin_grid)[this->point_index[1]].y)
//	{
//		this->spvertex = this->point_index[2];
//		this->flatline[0] = this->point_index[0];
//		this->flatline[1] = this->point_index[1];
//		if ((*Tin_grid)[this->point_index[2]].y > (*Tin_grid)[this->point_index[1]].y)
//		{
//			this->trishape = 1;//平底
//		}
//		else this->trishape = 2;//平顶
//	}
//	else if ((*Tin_grid)[this->point_index[1]].y == (*Tin_grid)[this->point_index[2]].y)
//	{
//		this->spvertex = this->point_index[0];
//		this->flatline[0] = this->point_index[1];
//		this->flatline[1] = this->point_index[2];
//		if ((*Tin_grid)[this->point_index[0]].y > (*Tin_grid)[this->point_index[1]].y)
//		{
//			this->trishape = 1;//平底
//		}
//		else this->trishape = 2;//平顶
//	}
//	else if ((*Tin_grid)[this->point_index[0]].y == (*Tin_grid)[this->point_index[2]].y)
//	{
//		this->spvertex = this->point_index[1];
//		this->flatline[0] = this->point_index[2];
//		this->flatline[1] = this->point_index[0];
//		if ((*Tin_grid)[this->point_index[1]].y > (*Tin_grid)[this->point_index[1]].y)
//		{
//			this->trishape = 1;//平底
//		}
//		else this->trishape = 2;//平顶
//	}
//	else//不平底或平顶
//	{
//		this->trishape = 0;
//		if ((*Tin_grid)[this->point_index[0]].y < (*Tin_grid)[this->point_index[1]].y)
//		{
//			if ((*Tin_grid)[this->point_index[0]].y > (*Tin_grid)[this->point_index[2]].y)//2在最下面
//			{
//				this->spvertex = this->point_index[0];
//			}
//			else if ((*Tin_grid)[this->point_index[2]].y > (*Tin_grid)[this->point_index[1]].y)//2在最上面
//			{
//				this->spvertex = this->point_index[1];
//			}
//			else this->spvertex = this->point_index[2];//2在中间
//		}
//		else
//		{
//			if ((*Tin_grid)[this->point_index[1]].y > (*Tin_grid)[this->point_index[2]].y)//2在最下面
//			{
//				this->spvertex = this->point_index[1];
//			}
//			else if ((*Tin_grid)[this->point_index[2]].y > (*Tin_grid)[this->point_index[0]].y)//2在最上面
//			{
//				this->spvertex = this->point_index[0];
//			}
//			else this->spvertex = this->point_index[2];//2在中间
//		}
//	}
//}