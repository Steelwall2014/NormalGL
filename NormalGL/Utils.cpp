#define NOMINMAX
#include <windows.h>
#include <math.h>
#include <vector>
#include <assert.h>
#include <locale.h>
#include <iostream>
#include <fstream>
#include <string>
#include <float.h>
#include <cmath>
#include <limits>
#include <cstdlib>


#include "utils.h"
#include "Graphic.h"
#include "Status.h"
#include "Matrix.h"
#include "Painters.h"
#include "3DModel.h"


using namespace std;

bool RectCross(double xmin1, double xmax1, double ymin1, double ymax1, double xmin2, double xmax2, double ymin2, double ymax2)
{
	int zx = fabs(xmin1 + xmax1 - xmin2 - xmax2);
	int x = fabs(xmin1 - xmax1) + fabs(xmin2 - xmax2);
	int zy = fabs(ymin1 + ymax1 - ymin2 - ymax2);
	int y = fabs(ymin1 - ymax1) + fabs(ymin2 - ymax2);
	if (zx <= x && zy <= y)
		return true;
	else
		return false;
}

Color get_random_color()
{
	int r = randint(0, 255), g = randint(0, 255), b = randint(0, 255);
	return _RGB(r, g, b);
}

int randint(int low, int high)
{
	return (rand() % (high - low + 1)) + low;
}
double radians(double degree)
{
	return degree / 180.0 * PI;
}
double degrees(double radian)
{
	return radian / PI * 180.0;
}
double get_direction(double x1, double y1, double x2, double y2)
{
	double dy = y2 - y1, dx = x2 - x1, degree;
	degree = atan2(dy, dx) * 180 / PI;
	if (degree < 0.0)
		return degree + 360.0;
	else
		return degree;
}
void get_mouse_pos(int &x, int &y, LPARAM lParam)
{
	x = LOWORD(lParam);
	y = HIWORD(lParam);
	DPtToLPt(x, y, x, y);//�����豸����תΪ�����߼�����
}
void set_unit_pixel(int x, int y, Color color)
{
	/*
	* set_unit: ����ĳһ����С��Ԫ����С��Ԫ������һ�����أ�����ģʽ�£���Ҳ������һ��С���ӣ�����ģʽ�£�
	* ע�������������½�����x��y���Ǵ����߼����꣬��λΪһ������
	* params:
	*	x: �����߼�����
	*	y: �����߼�����
	*	color: ��ɫ������ΪGraphic.cpp�е�_RGB
	* ���ߣ��ž���
	*/
	//Matrix m(1, 3);
	//m.data[0][0] = x; m.data[0][1] = y; m.data[0][2] = 1;
	//Matrix n = m * trans;
	//setPixel(n.data[0][0], n.data[0][1], color);
	setPixel(x, y, color);
}
void set_unit_grid(int grid_x, int grid_y, Color color)
{
	/*
	* set_unit: ����ĳһ����С��Ԫ����С��Ԫ������һ�����أ�����ģʽ�£���Ҳ������һ��С���ӣ�����ģʽ�£�
	* ע�������������½�����x��y���Ǵ����߼����꣬��λΪһ������
	* params:
	*	grid_x: �����к�
	*	grid_y: �����к�
	*	color: ��ɫ������ΪGraphic.cpp�е�_RGB
	* ���ߣ��ž���
	*/
	for (int i = grid_x * status.GetCellSize(); i < grid_x * status.GetCellSize() + status.GetCellSize(); i++)
		for (int j = grid_y * status.GetCellSize(); j < grid_y * status.GetCellSize() + status.GetCellSize(); j++)
			setPixel(i, j, color);
}
Color get_unit_pixel(int x, int y)
{
	/*
	* get_unit: ��ȡĳһ����С��Ԫ����ɫ����С��Ԫ������һ�����أ�����ģʽ�£���Ҳ������һ��С���ӣ�����ģʽ�£�
	* ע�������������½�����x��y���Ǵ������꣬��λΪһ������
	* Ŀǰ�����Ǽ򵥵ص���getPixel��������ϣ���Ժ���Ըĳ�x��y������ģʽ�±�ʾ���꣬�ڸ���ģʽ�±�ʾ�������к�
	* params:
	*	x: �����߼�����
	*	y: �����߼�����
	* ���ߣ��ž���
	*/
	return getPixel(x, y);
}
Color get_unit_grid(int grid_x, int grid_y)
{
	/*
	* get_unit: ��ȡĳһ����С��Ԫ����ɫ����С��Ԫ������һ�����أ�����ģʽ�£���Ҳ������һ��С���ӣ�����ģʽ�£�
	* ע�������������½�����x��y���Ǵ������꣬��λΪһ������
	* Ŀǰ�����Ǽ򵥵ص���getPixel��������ϣ���Ժ���Ըĳ�x��y������ģʽ�±�ʾ���꣬�ڸ���ģʽ�±�ʾ�������к�
	* params:
	*	x: ��������
	*	y: ��������
	* ���ߣ��ž���
	*/

	return getPixel(grid_x * status.GetCellSize() + 1, grid_y * status.GetCellSize() + 1);
}

float GetLength(int x1, int y1, int x2, int y2)
{
	int dx = x2 - x1, dy = y2 - y1;
	return sqrt(dx * dx + dy * dy);
}
double GetLength(double x1, double y1, double x2, double y2)
{
	double dx = x2 - x1, dy = y2 - y1;
	return sqrt(dx * dx + dy * dy);
}

void swap(int &x1, int &x2)
{
	x1 += x2;
	x2 = x1 - x2;
	x1 = x1 - x2;
}


AETNode::AETNode(Point2D p1, Point2D p2)
{
	double x1 = p1.x, x2 = p2.x;
	int y1 = (int)p1.y, y2 = (int)p2.y;
	if (y2 == y1)//�߶�ƽ�У�0��������ĸ��������
	{
		this->rk = DBL_MAX;
		this->yM = y1;
		this->ymin = y1;
		if (x1 < x2) this->x = x1;
		else this->x = x2;
	}
	else if (x1 == x2)
	{
		this->rk = 0.0;
		this->yM = max(y1, y2);
		this->ymin = min(y1, y2);
		this->x = x1;
	}
	else
	{
		double xm = x2 - x1;//���ں����ж�x��С
		double k = xm / (y2 - y1);//б�ʵĵ���
		this->rk = k;
		//int sign = k * xm;
		if (k > 0)
		{
			if (xm > 0)
			{
				this->x = x1;
				this->yM = y2;
				this->ymin = y1;
			}
			else
			{
				this->x = x2;
				this->yM = y1;
				this->ymin = y2;
			}
		}
		else if (k < 0)
		{
			if (xm > 0)
			{
				this->x = x2;
				this->yM = y1;
				this->ymin = y2;
			}
			else
			{
				this->x = x1;
				this->yM = y2;
				this->ymin = y1;
			}
		}
	}
}
bool CheckAET(list<AETNode> &aet)
{
	if (aet.empty())
		return true;
	list<AETNode>::iterator iter, end_iter = aet.end();
	end_iter--;
	for (iter = aet.begin(); iter != end_iter; )
	{
		AETNode node1 = *iter++, node2 = *iter;
		if (!cmp_aetnode(node1, node2))
		{
			return false;
		}
	}
	return true;
}


Bucket::Bucket(int ymin, int ymax)
{
	this->ymin = ymin;
	this->ymax = ymax;
	this->scanline_count = ymax - ymin + 1;
	this->edge_tables = new list<AETNode>[scanline_count];
}
Bucket::Bucket(Point2D *points, int point_count)
{
	int y1 = numeric_limits<int>::max(), y2 = numeric_limits<int>::min();//y1��Сֵ��y2���ֵ
	for (int i = 0; i < point_count; i++)
	{
		y1 = min((int)points[i].y, y1);
		y2 = max((int)points[i].y, y2);
	}
	this->ymin = y1; this->ymax = y2;
	this->scanline_count = y2 - y1 + 1;
	this->edge_tables = new list<AETNode>[scanline_count];
	points[point_count] = points[0];//ʹ֮��β����
	for (int i = 0; i < point_count; i++)
	{
		if ((int)points[i].y == (int)points[i + 1].y)
			continue;
		AETNode line = AETNode(points[i], points[i + 1]);
		//if (line->rk == DBL_MAX)
		//{
		//	delete line;
		//	continue;
		//}

		if (this->edge_tables[line.ymin - this->ymin].empty())
			this->min_ys.push_back(line.ymin - this->ymin);
		this->edge_tables[line.ymin - this->ymin].push_back(line);
	}

	//return buck;
}
Bucket::~Bucket()
{
	for (int i = 0; i < scanline_count; i++)
	{
		this->edge_tables[i].clear();
	}
	delete[] this->edge_tables;
}

void Bucket::Sort()
{
	//list<AETNode*> &aet = edge_tables[0];
	for (int i = 0; i < scanline_count; i++)
	{
		//aet = edge_tables[i];
		edge_tables[i].sort(cmp_aetnode);
	}
}
void Bucket::Correct()
{
	for (int min_y : min_ys)
	{
		for (AETNode &aet_node : edge_tables[min_y])
		{
			if (!edge_tables[aet_node.yM - this->ymin].empty())
				aet_node.yM--;
		}
	}
}
void Bucket::Draw(Color color)
{
	list<AETNode> aet;
	list<AETNode>::iterator iter, end_iter;
	AETNode *curr_aetnode, *next_aetnode;
	int flag = 0;//flagΪ1ʱ���Ƹ�x����һ�ڵ�x
	int start, finish, j = 0;
	int scanline_index;
	for (int i = 0; i < this->scanline_count; i++)	// ����������1�ͻ������е����һ��ɨ���ߵ�ʱ�򱨴�ԭ�������ˣ����������һ�л��Ա߱���������
	{
		scanline_index = i + this->ymin;
		//if (scanline_index == 346)
		//	cout << 1;
		flag = 0;
		//if (this->edge_tables[i].empty() == false)
		if (CheckAET(aet))
			aet.merge(this->edge_tables[i], cmp_aetnode);
		// �ҵ���߱����ԭ���ˣ�����Ϊaet������������aet�������������ԭ�����ڣ�����һ��ɨ����aet������ģ�����x����rk�����������������һ�����
		// ԭ��xС���Ǹ�������һ���֮������Դ���ԭ��x����Ǹ���ͬʱ��rk����Ļ�Ҳ��������ģ����³�Ϊ�����������
		if (aet.empty())
			continue;
		end_iter = aet.end();
		end_iter--;
		j = 0;
		for (iter = aet.begin(); iter != end_iter; )
		{
			flag = !flag;
			curr_aetnode = &(*iter);
			iter++;
			next_aetnode = &(*iter);
			iter--;
			if (flag && i != 0)
			{
				start = min(curr_aetnode->x, next_aetnode->x);
				finish = max(curr_aetnode->x, next_aetnode->x);
				for (int x = start; x < finish; x++)
					set_unit(x, scanline_index, color);
			}
			curr_aetnode->x += curr_aetnode->rk;
			if (curr_aetnode->yM == scanline_index)
				iter = aet.erase(iter);
			else
				iter++;
		}
		// ���һ��Ԫ�ص�������һ��
		curr_aetnode = &(*end_iter);
		if (curr_aetnode->yM == scanline_index)
			aet.erase(end_iter);
		else
			curr_aetnode->x += curr_aetnode->rk;
	}
}




void PixelToGrid(double x, double y, double &new_x, double &new_y)
{
	new_x = int(x / status.GetCellSize());
	new_y = int(y / status.GetCellSize());
}

void PixelToGrid(Point2D *points, Point2D *new_points, int point_count)
{
	for (int i = 0; i < point_count; i++)
	{
		PixelToGrid(points[i].x, points[i].y, new_points[i].x, new_points[i].y);
	}
}

void PixelToGrid(vector<Point2D> points, vector<Point2D> &new_points)
{
	new_points.clear();
	double x, y;
	for (int i = 0; i < points.size(); i++)
	{
		PixelToGrid(points[i].x, points[i].y, x, y);
		new_points.push_back(Point2D{ x, y });
	}
}


//list<AETNode *> *list_merge(list<AETNode *> l1, list<AETNode *> l2)
//{
//
//}
bool cmp_aetnode(AETNode &node1, AETNode &node2)
{
	if (node1.x < node2.x) return true;
	if (node1.x > node2.x) return false;
	return node1.rk <= node2.rk;

	//if (fabs(node1->x - node2->x) < 0.00000001)
	//	return node1->rk <= node2->rk;
	//else
	//	return node1->x <= node2->x;
}
void quickSort(PixelPoint *arr, int left, int right)//�Ե����鰴��y��������
{
	if (left > right) {
		return;
	}
	int i = left;
	int j = right;
	PixelPoint base = arr[left];
	while (i != j) {
		while (arr[j].y >= base.y && i < j) {
			j--;
		}

		while (arr[i].y <= base.y && i < j) {
			i++;
		}

		if (i < j) {
			PixelPoint tmp = arr[i];
			arr[i] = arr[j];
			arr[j] = tmp;
		}
	}

	arr[left] = arr[i];
	arr[i] = base;

	quickSort(arr, left, i - 1);
	quickSort(arr, i + 1, right);
}
void sortPoints(PixelPoint *points, int point_count)//�Ե����鰴��y��������
{
	quickSort(points, 0, point_count - 1);
}


size_t find_first_not_delim(const string &s, char delim, size_t pos)
{
	for (size_t i = pos; i < s.size(); i++)
		if (s[i] != delim)
			return i;
	return string::npos;
}
vector<string> split(const string &s, char delim = ' ')
{
	vector<string> tokens;
	size_t lastPos = find_first_not_delim(s, delim, 0);
	size_t pos = s.find(delim, lastPos);
	while (lastPos != string::npos)
	{
		tokens.push_back(s.substr(lastPos, pos - lastPos));
		lastPos = find_first_not_delim(s, delim, pos);
		pos = s.find(delim, lastPos);
	}
	return tokens;
}
FileHead *ReadHead(string hdr_filename)
{
	fstream head_file;
	string line;
	vector<string> *lines = new vector<string>;
	FileHead *file_head = new FileHead;
	head_file.open(hdr_filename, ios::in);
	while (getline(head_file, line))
	{
		lines->push_back(line);
		vector<string> tokens = split(line, ' ');
		string value = tokens[tokens.size() - 1];
		if (tokens[0] == "samples")
			file_head->nCol = std::stoi(value);
		else if (tokens[0] == "lines")
			file_head->nRow = std::stoi(value);
		else if (tokens[0] == "bands")
			file_head->nBand = std::stoi(value);
		else if (tokens[0] == "interleave")
			file_head->imgformat = value;
		else if (tokens[0] == "pixel")
		{
			value = value.substr(0, value.size() - 1);
			file_head->CellSize_x = atof(tokens[3].substr(1, tokens[3].size() - 2).c_str());
			file_head->CellSize_y = atof(tokens[4].substr(0, tokens[4].size() - 1).c_str());
			file_head->units = split(value, '=')[1];
		}
	}
	head_file.close();
	file_head->hdr_lines = lines;
	return file_head;
}
Image *ReadRaster(string filename, FileHead *file_head)
{
	Image *img = new Image;
	fstream raster_file;
	int row = file_head->nRow, col = file_head->nCol, band = file_head->nBand;
	streamsize count = streamsize(row * col * band);
	raster_file.open(filename, ios::in | ios::binary);
	img->array = new int[count];
	char *temp = new char[count];
	raster_file.read(temp, count);
	for (int i = 0; i < count; i++)
		img->array[i] = (int)temp[i];
	img->filehead = file_head;
	raster_file.close();
	return img;
}
Image *ReadData(string filename)
{
	Image *img;
	FileHead *file_head;
	file_head = ReadHead(filename + string(".hdr"));
	img = ReadRaster(filename, file_head);
	return img;
}


void WriteHead(string filename, FileHead *filehead)
{
	fstream head_file;
	vector<string> *lines = filehead->hdr_lines;
	head_file.open(filename.c_str(), ios::out);
	for (size_t i = 0; i < lines->size(); i++)
	{
		string line = (*lines)[i];
		vector<string> tokens = split(line, ' ');
		string value = tokens[tokens.size() - 1];
		if (tokens[0] == "interleave")
			line = "interleave = " + filehead->imgformat;
		else if (tokens[0] == "description")
			(*lines)[i + 1] = "  File Exported from cpp}";
		else if (tokens[0] == "bands")
			line = "bands   = " + to_string(filehead->nBand);
		head_file << line << '\n';
	}
}
void WriteRaster(string filename, Image *img)
{
	fstream raster_file;
	int row = img->filehead->nRow, col = img->filehead->nCol, band = img->filehead->nBand;
	streamsize count = row * col * band;
	char *temp = new char[count];
	for (int i = 0; i < count; i++)
		temp[i] = img->array[i];
	raster_file.open(filename.c_str(), ios::out | ios::binary);
	raster_file.write(temp, count);
}
void WriteData(string filename, Image *img)
{
	WriteHead(filename + ".hdr", img->filehead);
	WriteRaster(filename, img);
}

Point2D *WindowCoordinateTransform(Point2D *points, int point_count, Matrix &trans)
{	// ����ĵ�����Ҫ��new
	Matrix point_matrix = Matrix(point_count, 3, 1);
	for (int i = 0; i < point_count; i++)
	{
		point_matrix.data[i][0] = points[i].x;
		point_matrix.data[i][1] = points[i].y;
	}
	Point2D *new_points = new Point2D[point_count];
	Matrix new_point_matrix = point_matrix * trans;
	for (int i = 0; i < point_count; i++)
	{
		new_points[i] = Point2D{ new_point_matrix.data[i][0], new_point_matrix.data[i][1] };
	}
	return new_points;
}

vector<Painter *> FromOGRGeometriesToPainters(vector<OGRGeometry *> &geom_set, Matrix &geom_to_map_trans)
{
	vector<Painter *> painter_set;
	for (auto geom : geom_set)
	{
		switch (geom->getGeometryType())
		{
		case wkbLineString:
		{
			OGRLineString *line_geom = (OGRLineString *)geom;
			int point_count = line_geom->getNumPoints();
			Point2D *points = new Point2D[point_count], *after_trans_points;
			for (int i = 0; i < point_count; i++)
			{
				points[i].x = line_geom->getX(i);
				points[i].y = line_geom->getY(i);
			}
			after_trans_points = WindowCoordinateTransform(points, point_count, geom_to_map_trans);
			Painter *painter = new PolylinePainter(after_trans_points, point_count);
			painter_set.push_back(painter);
			delete[] points;
			delete[] after_trans_points;
			break;
		}
		case wkbMultiLineString:
		{
			OGRMultiLineString *mline_geom = (OGRMultiLineString *)geom;
			int geom_count = mline_geom->getNumGeometries();
			for (int i = 0; i < geom_count; i++)
			{
				OGRLineString *line_geom = (OGRLineString *)mline_geom->getGeometryRef(i);
				int point_count = line_geom->getNumPoints();
				Point2D *points = new Point2D[point_count], *after_trans_points;
				for (int i = 0; i < point_count; i++)
				{
					points[i].x = line_geom->getX(i);
					points[i].y = line_geom->getY(i);
				}
				after_trans_points = WindowCoordinateTransform(points, point_count, geom_to_map_trans);
				painter_set.push_back(new PolylinePainter(after_trans_points, point_count));
				delete[] points;
				delete[] after_trans_points;
				//delete line_geom;
			}
			//delete mline_geom;
			break;
		}
		case wkbPolygon:
		{
			OGRPolygon *polygon_geom = (OGRPolygon *)geom;
			OGRLinearRing *ring = polygon_geom->getExteriorRing();
			int point_count = ring->getNumPoints();
			Point2D *points = new Point2D[point_count], *after_trans_points;
			for (int i = 0; i < point_count; i++)
			{
				points[i].x = ring->getX(i);
				points[i].y = ring->getY(i);
			}
			after_trans_points = WindowCoordinateTransform(points, point_count, geom_to_map_trans);
			painter_set.push_back(new PolygonPainter(after_trans_points, point_count - 1));
			delete[] points;
			delete[] after_trans_points;
			//delete ring;
			break;
		}
		}
	}
	return painter_set;
}

vector<Face> PresentationData_3D()
{
	vector<Face> solid = vector<Face>();
	return solid;
}
Matrix *GetPlaneNormalVector(Point3D p1, Point3D p2, Point3D p3)
{
	double na = (p2.y - p1.y) * (p3.z - p1.z) - (p2.z - p1.z) * (p3.y - p1.y);
	double nb = (p2.z - p1.z) * (p3.x - p1.x) - (p2.x - p1.x) * (p3.z - p1.z);
	double nc = (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
	Matrix *norm = new Matrix(1, 3);
	norm->data[0][0] = -na;
	norm->data[0][1] = -nb;
	norm->data[0][2] = -nc;
	norm->Normalize();
	return norm;
}
Point3D GetPlaneNormalVector_Point3D(Point3D p1, Point3D p2, Point3D p3)
{
	double na = (p2.y - p1.y) * (p3.z - p1.z) - (p2.z - p1.z) * (p3.y - p1.y);
	double nb = (p2.z - p1.z) * (p3.x - p1.x) - (p2.x - p1.x) * (p3.z - p1.z);
	double nc = (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
	Point3D norm;
	norm.x = -na;	// ��֪��Ϊʲô�Ӹ����ž�������ʵ����ɫ��Ч������������ͬ
	norm.y = -nb;
	norm.z = -nc;
	norm.Normalize();
	return norm;
}
//LUP�ֽ�
void LUP_Descomposition(Matrix A, Matrix &L, Matrix &U, Matrix &P)
{
	int num = A.shape[1];
	int row = 0;
	for (int i = 0; i < num; i++)
	{
		P.data[0][i] = i;
	}
	double p = 0.0;
	for (int i = 0; i < num - 1; i++)
	{
		for (int j = i; j < num; j++)
		{
			if (abs(A.data[j][i]) > p)
			{
				p = abs(A.data[j][i]);
				row = j;
			}
		}
		if (0 == p)
		{
			throw "�������죬�޷�������";
		}
		//����P[i]��P[row]
		int tmp = P.data[0][i];
		P.data[0][i] = P.data[0][row];
		P.data[0][row] = tmp;

		double tmp2 = 0.0;
		for (int j = 0; j < num; j++)
		{
			//����A.data[i][j]�� A.data[row][j]
			tmp2 = A.data[i][j];
			A.data[i][j] = A.data[row][j];
			A.data[row][j] = tmp2;
		}
		//����ͬLU�ֽ�
		double u = A.data[i][i], l = 0.0;
		for (int j = i + 1; j < num; j++)
		{
			l = A.data[j][i] / u;
			A.data[j][i] = l;
			for (int k = i + 1; k < num; k++)
			{
				A.data[j][k] = A.data[j][k] - A.data[i][k] * l;
			}
		}

	}
	//����L��U
	for (int i = 0; i < num; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			if (i != j)
			{
				L.data[i][j] = A.data[i][j];
			}
			else
			{
				L.data[i][j] = 1;
			}
		}
		for (int k = i; k < num; k++)
		{
			U.data[i][k] = A.data[i][k];
		}
	}
}

//LUP��ⷽ��
Matrix LUP_Solve(Matrix L, Matrix U, Matrix P, double b[])
{
	int num = L.shape[0];
	Matrix x(1, num);
	double *y = new double[num]();

	//�����滻
	for (int i = 0; i < num; i++)
	{
		y[i] = b[int(P.data[0][i])];
		for (int j = 0; j < i; j++)
		{
			y[i] = y[i] - L.data[i][j] * y[j];
		}
	}
	//�����滻
	for (int i = num - 1; i >= 0; i--)
	{
		x.data[0][i] = y[i];
		for (int j = num - 1; j > i; j--)
		{
			x.data[0][i] = x.data[0][i] - U.data[i][j] * x.data[0][j];
		}
		x.data[0][i] /= U.data[i][i];
	}
	return x;
}

//LUP����(��ÿ��b����ĸ���x������װ)
Matrix LUP_solve_inverse(Matrix A)
{
	//double* A_mirror = new double[N * N]();
	int num = A.shape[0];
	Matrix inv_A(num, num, 0);//���յ�����󣨻���Ҫת�ã�
	Matrix inv_A_each(1, num);//������ĸ���
	//double *B = new double[N*N]();
	double *b = new double[num]();//b��ΪB����о������

	for (int i = 0; i < num; i++)
	{
		Matrix L(num, num);
		Matrix U(num, num);
		Matrix P(1, num);

		//���쵥λ���ÿһ��
		for (int i = 0; i < num; i++)
		{
			b[i] = 0;
		}
		b[i] = 1;

		//ÿ�ζ���Ҫ���½�A����һ��
		/*for (int i = 0; i < N * N; i++)
		{
			A_mirror[i] = A[i];
		}*/

		LUP_Descomposition(A, L, U, P);

		inv_A_each = LUP_Solve(L, U, P, b);
		for (int j = 0; j < num; ++j)
		{
			inv_A.data[i][j] = inv_A_each.data[0][j];
		}
		//������ƴ������
	}
	//inv_A.Transpose();//�������ڸ���ÿ��b�����x���д洢�������ת��

	return inv_A;
}
