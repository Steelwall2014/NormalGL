#pragma once
#include "ogrsf_frmts.h"
#include "ogr_geometry.h"
//#include "Matrix.h"
//#include "..\\3rdParty\\gdal\\include\\ogrsf_frmts.h"
//#include "..\\3rdParty\\gdal\\include\\ogr_geometry.h"
#include <vector>
#include <algorithm>
#include <string>
using namespace std;
class Tin;
#define  PI  3.14159265358979323846

class Matrix;
//2D点
template<typename T>
struct _Point2D
{
	_Point2D(){ x = y = 0 ;}
	_Point2D( T x, T y ){ this->x = x, this->y = y; }

	T x,y;
};

typedef _Point2D<double> Point2D;
//
////2D包围盒
//template<typename T>
//struct _Box2D
//{
//	_Box2D(){ valid = false; }
//	_Box2D( T xmin,T ymin , T xmax, T ymax )
//	{
//		setBox(xmin, ymin, xmax, ymax);
//	}
//
//	T width(){ return _xmax - _xmin; }//包围盒宽度
//	T height(){ return _ymax - _ymin; }//包围盒高度
//	T centerX(){ return (_xmin + _xmax) /2 ; }//包围盒中心点x坐标
//	T centerY(){ return (_ymin + _ymax) /2 ; }//包围盒中心点y坐标
//	T xmin(){ return _xmin;}
//	T ymin(){ return _ymin;}
//	T xmax(){ return _xmax;}
//	T ymax(){ return _ymax;}
//
//	void setBox( T xmin,T ymin , T xmax, T ymax )//设置包围盒范围
//	{
//		this->_xmin = xmin, this->_xmax = xmax, this->_ymin = ymin, this->_ymax = ymax ;
//		valid = true;
//	}
//
//	void setBox( _Box2D& box )
//	{
//		setBox( box._xmin, box._ymin, box._xmax, box._ymax );
//	}
//
//	//根据添加的点扩展包围盒
//	void expand( T x, T y )
//	{
//		if( !valid ){
//			_xmin = _xmax = x;
//			_ymin = _ymax = y;
//			valid = true;
//		}
//		else
//		{
//			_xmin = min( _xmin, x );
//			_xmax = max( _xmax, x );
//			_ymin = min( _ymin, y );
//			_ymax = max( _ymax, y );
//		}
//	}
//
//	//根据添加的包围盒扩展包围盒
//	void expand( _Box2D<T>& box ){
//		if (box.valid)
//		{
//			if ( !valid )
//				setBox(box);
//			else
//			{
//				expand(box.xmin(), box.ymin());
//				expand(box.xmax(), box.ymax());
//			}
//		}
//	}
//
//protected:
//	T _xmin,_ymin;
//	T _xmax,_ymax;
//	bool valid;//是否有效
//};
//typedef _Box2D<double> Box2D;
//
////图形对象类型
//enum GeomType{ gtUnkown = 0, gtPoint = 1, gtPolyline = 2, gtPolygon = 3 };
//
////图形对象基类，供继承
//struct Geometry
//{
//	virtual ~Geometry(){}
//	virtual GeomType getGeomType() =  0; //获取图形对象类型
//	virtual Box2D getEnvelop() = 0 ;//获取图形对象包围盒
//	string label;
//};
//
////点图像对象
//struct PointGemetry:Geometry
//{
//	PointGemetry(  )
//	{
//		x = y = 0;
//	}
//	PointGemetry( double x, double y )
//	{
//		this->x= x;
//		this->y = y;
//	}
//
//	virtual GeomType getGeomType(){ return gtPoint; }
//	virtual Box2D getEnvelop(){ return Box2D( x,y,x,y ); }
//
//	double x,y;
//};
//
////线图形对象
//struct PolylineGeometry:Geometry
//{
//	virtual GeomType getGeomType(){ return gtPolyline; }
//
//	//运算符重载，获取第i个点
//	Point2D& operator[]( int i ){ return pts[i]; }
//
//	//添加点
//	void addPoint( double x, double y )
//	{
//		pts.push_back( Point2D( x, y ));
//		envelop.expand( x, y );
//	}
//
//	virtual Box2D getEnvelop(){ return envelop; }
//	const vector<Point2D>& getPts() const 
//	{
//		return pts;
//	}
//protected:
//	//所有点
//	vector<Point2D> pts;
//	Box2D envelop;
//};
//
////多边形图形对象
//struct PolygonGeometry:PolylineGeometry
//{
//	virtual GeomType getGeomType(){ return gtPolygon; }
//};

//图层
struct Layer
{
	Layer()
	{
		this->geomType = wkbNone;
	}

	Layer( OGRwkbGeometryType geomType )
	{
		this->geomType = geomType;
	}

	virtual ~Layer()
	{
		for( size_t i = 0 , size = geometrySet.size(); i < size ; ++i ) delete geometrySet[i];//析构函数删除所有图形对象
	}

	//运算符重载，返回第i个图形对象
	OGRGeometry* operator[]( int i ){ return geometrySet[i]; }

	//添加图像对象
	void addGeometry( OGRGeometry* pGeometry, bool updateEnvelop = false )
	{
		if (geometrySet.size() == 0)
		{
			layer_srs = pGeometry->getSpatialReference();
		}
		geometrySet.push_back( pGeometry );
		if (updateEnvelop)
		{
			OGREnvelope pEnvelope;
			pGeometry->getEnvelope(&pEnvelope);
			envelop.Merge(pEnvelope);
		}
			
	}

	void reset()
	{
		for (int i = 0; i < geometrySet.size(); i++)
			OGRGeometryFactory::destroyGeometry(geometrySet[i]);//delete geometrySet[i];//ogr的删除
		geometrySet.clear();
	}
	////设置图层范围
	//void setEnvelop( double xmin,double ymin , double xmax, double ymax  )
	//{
	//	envelop.setBox( xmin, ymin,xmax,ymax );
	//}

	//获取图层包含的图形对象数量
	int getGeometryCount(){ return geometrySet.size(); }
	//获取图层包围盒
	OGREnvelope getEnvelop() { return envelop; }

	vector<OGRGeometry*> geometrySet;//图形对象集合
	OGREnvelope envelop;//图层范围对应的包围盒
	OGRwkbGeometryType geomType;//图层类型
	OGRSpatialReference *layer_srs;
};

//数据集
struct Dataset
{
	virtual ~Dataset()
	{
		for( size_t i = 0, size = layerSet.size() ; i < size ; ++i ) delete layerSet[i];//析构函数删除图层
	}

	//运算符重载，返回第i个图层
	Layer* operator[]( int i ){ return layerSet[i]; }

	int indexOf( Layer* pLayer ){
		for ( int i = 0 ; i < layerSet.size(); ++i )
		{
			if( pLayer = layerSet[i]) return i;
		}
		return -1;
	}

	//获取图层数
	int getLayerCount(){ return layerSet.size(); }	

	//添加图层
	void addLayer( Layer* pLayer )
	{
		layerSet.push_back( pLayer );
		dataset_envelope.Merge(pLayer->envelop);
	}

	void reset()
	{
		for (int i = 0; i < layerSet.size(); i++)
			layerSet[i]->reset();
		layerSet.clear();
	}
	//图层集合
	vector<Layer*> layerSet;
	OGREnvelope dataset_envelope;
};

class Point3D : Point2D
{
public:
	int id;
	double x;
	double y;
	double z;
	double depth = 0;
	Point3D()
	{
		this->x = 0;
		this->y = 0;
		this->z = 0;
	}
	Point3D(const Point2D &point2d)
	{
		this->x = point2d.x;
		this->y = point2d.y;
		this->z = 0;
	}
	Point3D(double x, double y, double z)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}
	Point3D operator*(Matrix &trans);
	void Normalize();
	double DotProduct(Matrix &other);
};

class Line3D
{
public:
	Point3D pts[2];
	Line3D() {}
	Line3D(Point3D pt1, Point3D pt2)
	{
		this->pts[0] = pt1;
		this->pts[1] = pt2;
	}
};