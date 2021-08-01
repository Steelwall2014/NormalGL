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
//2D��
template<typename T>
struct _Point2D
{
	_Point2D(){ x = y = 0 ;}
	_Point2D( T x, T y ){ this->x = x, this->y = y; }

	T x,y;
};

typedef _Point2D<double> Point2D;
//
////2D��Χ��
//template<typename T>
//struct _Box2D
//{
//	_Box2D(){ valid = false; }
//	_Box2D( T xmin,T ymin , T xmax, T ymax )
//	{
//		setBox(xmin, ymin, xmax, ymax);
//	}
//
//	T width(){ return _xmax - _xmin; }//��Χ�п��
//	T height(){ return _ymax - _ymin; }//��Χ�и߶�
//	T centerX(){ return (_xmin + _xmax) /2 ; }//��Χ�����ĵ�x����
//	T centerY(){ return (_ymin + _ymax) /2 ; }//��Χ�����ĵ�y����
//	T xmin(){ return _xmin;}
//	T ymin(){ return _ymin;}
//	T xmax(){ return _xmax;}
//	T ymax(){ return _ymax;}
//
//	void setBox( T xmin,T ymin , T xmax, T ymax )//���ð�Χ�з�Χ
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
//	//������ӵĵ���չ��Χ��
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
//	//������ӵİ�Χ����չ��Χ��
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
//	bool valid;//�Ƿ���Ч
//};
//typedef _Box2D<double> Box2D;
//
////ͼ�ζ�������
//enum GeomType{ gtUnkown = 0, gtPoint = 1, gtPolyline = 2, gtPolygon = 3 };
//
////ͼ�ζ�����࣬���̳�
//struct Geometry
//{
//	virtual ~Geometry(){}
//	virtual GeomType getGeomType() =  0; //��ȡͼ�ζ�������
//	virtual Box2D getEnvelop() = 0 ;//��ȡͼ�ζ����Χ��
//	string label;
//};
//
////��ͼ�����
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
////��ͼ�ζ���
//struct PolylineGeometry:Geometry
//{
//	virtual GeomType getGeomType(){ return gtPolyline; }
//
//	//��������أ���ȡ��i����
//	Point2D& operator[]( int i ){ return pts[i]; }
//
//	//��ӵ�
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
//	//���е�
//	vector<Point2D> pts;
//	Box2D envelop;
//};
//
////�����ͼ�ζ���
//struct PolygonGeometry:PolylineGeometry
//{
//	virtual GeomType getGeomType(){ return gtPolygon; }
//};

//ͼ��
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
		for( size_t i = 0 , size = geometrySet.size(); i < size ; ++i ) delete geometrySet[i];//��������ɾ������ͼ�ζ���
	}

	//��������أ����ص�i��ͼ�ζ���
	OGRGeometry* operator[]( int i ){ return geometrySet[i]; }

	//���ͼ�����
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
			OGRGeometryFactory::destroyGeometry(geometrySet[i]);//delete geometrySet[i];//ogr��ɾ��
		geometrySet.clear();
	}
	////����ͼ�㷶Χ
	//void setEnvelop( double xmin,double ymin , double xmax, double ymax  )
	//{
	//	envelop.setBox( xmin, ymin,xmax,ymax );
	//}

	//��ȡͼ�������ͼ�ζ�������
	int getGeometryCount(){ return geometrySet.size(); }
	//��ȡͼ���Χ��
	OGREnvelope getEnvelop() { return envelop; }

	vector<OGRGeometry*> geometrySet;//ͼ�ζ��󼯺�
	OGREnvelope envelop;//ͼ�㷶Χ��Ӧ�İ�Χ��
	OGRwkbGeometryType geomType;//ͼ������
	OGRSpatialReference *layer_srs;
};

//���ݼ�
struct Dataset
{
	virtual ~Dataset()
	{
		for( size_t i = 0, size = layerSet.size() ; i < size ; ++i ) delete layerSet[i];//��������ɾ��ͼ��
	}

	//��������أ����ص�i��ͼ��
	Layer* operator[]( int i ){ return layerSet[i]; }

	int indexOf( Layer* pLayer ){
		for ( int i = 0 ; i < layerSet.size(); ++i )
		{
			if( pLayer = layerSet[i]) return i;
		}
		return -1;
	}

	//��ȡͼ����
	int getLayerCount(){ return layerSet.size(); }	

	//���ͼ��
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
	//ͼ�㼯��
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