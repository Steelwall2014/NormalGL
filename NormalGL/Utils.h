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
	double x;//端点中较小x坐标
	int ymin;//端点较小y坐标
	int yM;//端点最大y坐标
	double rk;//线段斜率倒数
	AETNode(Point2D p1, Point2D p2);
	//~AETNode();
};
// aet算法的桶
class Bucket
{
public:
	int ymin, ymax, scanline_count;//当前扫描线的y值
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

// 栅格的头文件，已弃用
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
// 栅格，已弃用
class Image
{
public:
	FileHead *filehead;
	int *array;
	int GetPixelCount() { return this->filehead->nBand * this->filehead->nCol * this->filehead->nRow; };
};

// 种子点，填充时用的
struct SeedPoint
{
	Point2D point;
	Color seed_color;
};


bool RectCross(double xmin1, double xmax1, double ymin1, double ymax1, 
			   double xmin2, double xmax2, double ymin2, double ymax2);
/**	返回随机颜色
@return 随机的颜色
*/
Color get_random_color();
/**	返回随机数
@param  low, high值域，包括两端
@return 随机数
*/
int randint(int low, int high);
/**	角度转弧度
@return 弧度, double
*/
double radians(double degree);
/**	弧度转角度
@return 角度, double
*/
double degrees(double radian);
/**	获取两点间距离 汤子豪
@param  x1 第一个点x坐标
@param  y1 第一个点y坐标
@param  x2 第二个点x坐标
@param  y2 第二个点y坐标
@return 距离
*/
float GetLength(int x1, int y1, int x2, int y2);
double GetLength(double x1, double y1, double x2, double y2);
/**	交换两个数 汤子豪
@return void
*/
void swap(int &x1, int &x2);
/**	获得鼠标窗口位置
*/
void get_mouse_pos(int &x, int &y, LPARAM lParam);
/**	获取向量(x2-x1, y2-y1)的方向
@return 方向，double，值域为0-360，0为x轴正方向
*/
double get_direction(double x1, double y1, double x2, double y2);
/**	设置像素颜色
@param  x 窗口坐标x
@param  y 窗口坐标y
*/
void set_unit_pixel(int x, int y, Color color);
/**	设置格子颜色
@param  x 格子坐标x
@param  y 格子坐标y
*/
void set_unit_grid(int x, int y, Color color);
/**	获取像素颜色
@param  x 窗口坐标x
@param  y 窗口坐标y
@return 颜色，Color
*/
Color get_unit_pixel(int x, int y);
/**	获取格子颜色
@param  x 格子坐标x
@param  y 格子坐标y
@return 颜色，Color
*/
Color get_unit_grid(int grid_x, int grid_y);
/**	比较两个活性边表结点的大小
@param  node1 一个活性边表结点
@param  node2 另一个活性边表结点
@return bool
*/
bool cmp_aetnode(AETNode &node1, AETNode &node2);
/**	按照字符分割字符串
@param  s 待分割字符串
@param  delim 分割字符
@return 分割后的字符串的每一个部分
*/
vector<string> split(const string &s, char delim);
FileHead *ReadHead(string hdr_filename);
Image *ReadRaster(string filename, FileHead *file_head);
Image *ReadData(string filename);
void WriteHead(string filename, FileHead *filehead);
void WriteRaster(string filename, Image *img);
void WriteData(string filename, Image *img);

/**	窗口坐标转格子坐标
@param  x 像素坐标
@param  y 像素坐标
@param  new_x 格网坐标
@param  new_y 格网坐标
*/
void PixelToGrid(double x, double y, double &new_x, double &new_y);
/**	窗口坐标转格子坐标
@param  points 像素坐标
@param  new_points 格网坐标
@param  point_count 点的数目
*/
void PixelToGrid(Point2D *points, Point2D *new_points, int point_count);
// 给定点和变换矩阵，输出变换后的矩阵
Point2D *WindowCoordinateTransform(Point2D *points, int point_count, Matrix &trans);
// 将OGRGeometry转换为painter
vector<Painter *> FromOGRGeometriesToPainters(vector<OGRGeometry *> &geom_set, Matrix &geom_to_map_trans);
// 获取平面法向量
Matrix *GetPlaneNormalVector(Point3D p1, Point3D p2, Point3D p3);
Point3D GetPlaneNormalVector_Point3D(Point3D p1, Point3D p2, Point3D p3);
// 没用上
vector<Face> PresentationData_3D();

// 算逆矩阵
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

