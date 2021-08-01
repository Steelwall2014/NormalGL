#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
//#include <cstdlib>
#include <list>
#include <windows.h>

#include "Graphic.h"


using namespace std;

class Painter;
class Matrix;
class Face;
class AETNode
{
public:
	double x;//�˵��н�Сx����
	int ymin;//�˵��Сy����
	int yM;//�˵����y����
	double rk;//�߶�б�ʵ���
	AETNode(Point2D p1, Point2D p2);
	//~AETNode();
};
// aet�㷨��Ͱ
class Bucket
{
public:
	int ymin, ymax, scanline_count;//��ǰɨ���ߵ�yֵ
	vector<int> min_ys;
	list<AETNode> *edge_tables;
	//void Init(int ymin, int ymax);
	Bucket(int ymin, int ymax);
	Bucket(Point2D *points, int point_count);
	~Bucket();
	void Sort();
	void Correct();
	void Draw(Color fillcolor);
	//AET* operator[](int index)
	//{
	//	return aets[index - ymin];
	//}
};

// դ���ͷ�ļ���������
struct FileHead
{
	int nCol;
	int nRow;
	int nBand;
	double CellSize_x = 1.0;
	double CellSize_y = 1.0;
	float UpperLeft_x = 0.0;
	float UpperLeft_y = 0.0;
	string units = "Meters";
	string imgformat;
	vector<string> *hdr_lines;
};
// դ��������
class Image
{
public:
	FileHead *filehead;
	int *array;
	int GetPixelCount() { return this->filehead->nBand * this->filehead->nCol * this->filehead->nRow; };
};

// ���ӵ㣬���ʱ�õ�
struct SeedPoint
{
	Point2D point;
	Color seed_color;
};


bool RectCross(double xmin1, double xmax1, double ymin1, double ymax1, 
			   double xmin2, double xmax2, double ymin2, double ymax2);
/**	���������ɫ
@return �������ɫ
*/
Color get_random_color();
/**	���������
@param  low, highֵ�򣬰�������
@return �����
*/
int randint(int low, int high);
/**	�Ƕ�ת����
@return ����, double
*/
double radians(double degree);
/**	����ת�Ƕ�
@return �Ƕ�, double
*/
double degrees(double radian);
/**	��ȡ�������� ���Ӻ�
@param  x1 ��һ����x����
@param  y1 ��һ����y����
@param  x2 �ڶ�����x����
@param  y2 �ڶ�����y����
@return ����
*/
float GetLength(int x1, int y1, int x2, int y2);
double GetLength(double x1, double y1, double x2, double y2);
/**	���������� ���Ӻ�
@return void
*/
void swap(int &x1, int &x2);
/**	�����괰��λ��
*/
void get_mouse_pos(int &x, int &y, LPARAM lParam);
/**	��ȡ����(x2-x1, y2-y1)�ķ���
@return ����double��ֵ��Ϊ0-360��0Ϊx��������
*/
double get_direction(double x1, double y1, double x2, double y2);
/**	����������ɫ
@param  x ��������x
@param  y ��������y
*/
void set_unit_pixel(int x, int y, Color color);
/**	���ø�����ɫ
@param  x ��������x
@param  y ��������y
*/
void set_unit_grid(int x, int y, Color color);
/**	��ȡ������ɫ
@param  x ��������x
@param  y ��������y
@return ��ɫ��Color
*/
Color get_unit_pixel(int x, int y);
/**	��ȡ������ɫ
@param  x ��������x
@param  y ��������y
@return ��ɫ��Color
*/
Color get_unit_grid(int grid_x, int grid_y);
/**	�Ƚ��������Ա߱���Ĵ�С
@param  node1 һ�����Ա߱���
@param  node2 ��һ�����Ա߱���
@return bool
*/
bool cmp_aetnode(AETNode &node1, AETNode &node2);
/**	�����ַ��ָ��ַ���
@param  s ���ָ��ַ���
@param  delim �ָ��ַ�
@return �ָ����ַ�����ÿһ������
*/
vector<string> split(const string &s, char delim);
FileHead *ReadHead(string hdr_filename);
Image *ReadRaster(string filename, FileHead *file_head);
Image *ReadData(string filename);
void WriteHead(string filename, FileHead *filehead);
void WriteRaster(string filename, Image *img);
void WriteData(string filename, Image *img);

/**	��������ת��������
@param  x ��������
@param  y ��������
@param  new_x ��������
@param  new_y ��������
*/
void PixelToGrid(double x, double y, double &new_x, double &new_y);
/**	��������ת��������
@param  points ��������
@param  new_points ��������
@param  point_count �����Ŀ
*/
void PixelToGrid(Point2D *points, Point2D *new_points, int point_count);
// ������ͱ任��������任��ľ���
Point2D *WindowCoordinateTransform(Point2D *points, int point_count, Matrix &trans);
// ��OGRGeometryת��Ϊpainter
vector<Painter *> FromOGRGeometriesToPainters(vector<OGRGeometry *> &geom_set, Matrix &geom_to_map_trans);
// ��ȡƽ�淨����
Matrix *GetPlaneNormalVector(Point3D p1, Point3D p2, Point3D p3);
Point3D GetPlaneNormalVector_Point3D(Point3D p1, Point3D p2, Point3D p3);
// û����
vector<Face> PresentationData_3D();

// �������
void LUP_Descomposition(Matrix A, Matrix &L, Matrix &U, Matrix &P);

Matrix LUP_Solve(Matrix L, Matrix U, Matrix P, double b[]);

Matrix LUP_solve_inverse(Matrix A);
//{x = 246.00000000000000 y = 463.00000000000000 }
//{x = 104.00000000000000 y = 277.00000000000000 }
//{x = 482.00000000000000 y = 147.00000000000000 }
//{x = 725.00000000000000 y = 332.00000000000000 }
//{x = 593.00000000000000 y = 410.00000000000000 }
//{x = 360.00000000000000 y = 252.00000000000000 }
//{x = 242.00000000000000 y = 289.00000000000000 }
//{x = 254.00000000000000 y = 339.00000000000000 }
//{x = 449.00000000000000 y = 418.00000000000000 }
//{x = 478.00000000000000 y = 450.00000000000000 }
//{x = 356.00000000000000 y = 451.00000000000000 }
//{x = 301.00000000000000 y = 399.00000000000000 }
//{x = 246.00000000000000 y = 463.00000000000000 }

