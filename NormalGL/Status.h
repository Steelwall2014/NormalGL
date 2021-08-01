#pragma once
#include <windows.h>
#include <math.h>
#include <vector>
#include <locale.h>

#include "ogr_spatialref.h"
//#include "..\\3rdParty\\gdal\\include\\ogr_spatialref.h"
#include "Graphic.h"
#include "Painters.h"
#include "Utils.h"
#include "Filler.h"

using namespace std;

//���Ʋ�������
enum PaintType { ptNone = -999, ptDrawRaster = -1, ptDrawLine = 0, ptDrawPolyline, ptDrawPolygon, ptDrawCircle, ptDrawEllipse, ptDrawRect };

//��ʾģʽ
enum DisplayType { dtPixelMode = 1, dtGridMode_5 = 5, dtGridMode_10 = 10, dtGridMode_20 = 20 };

enum FillStatus { fsNone = 0, fsAutoAET, fsRecurse, fsStack, fsScanline };

enum OperationStatus { osNone, osTransform, osClip, os3D };// tsTranslation, tsRotate, tsZoom };

enum ClipType { ctCenterClip, ctSHClip, ctWAClip, ctClip };// �����㷨

enum Geom3DType { gtTeapot, gtPyramid,gtJunc }; // ��άģ������Щ
class Status
{
public:
	Status()
	{
		dttype = dtPixelMode;
		ptype = ptDrawLine;
		FillStatus = fsNone;
		GridColor = _RGB(220, 220, 220);
		OsStatus = osNone;
		//map_srs.importFromProj4("+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null +wktext +no_defs");
		map_srs.importFromEPSG(3857);
		geom_3d_type = gtPyramid;
	};

	//���û������ͱ����ߡ����ε�
	void SetPaintType(PaintType pType)
	{
		ptype = pType;
	}

	//������ʾģʽ�����ء����̸�
	void SetDisplayType(DisplayType dtType)
	{
		dttype = dtType;
	}

	//�����Ƿ����
	void SetFillStatus(FillStatus fill_status)
	{
		FillStatus = fill_status;
	}

	void SetGridColor(Color color)
	{
		GridColor = color;
	}
	PaintType GetPaintType()
	{
		return ptype;
	}
	DisplayType GetDisplayType()
	{
		return dttype;
	}
	int GetCellSize()
	{
		return (int)dttype;
	}
	FillStatus GetFillStatus()
	{
		return FillStatus;
	}
	Color GetGridColor()
	{
		return GridColor;
	}
	ClipType ctype;

	OperationStatus OsStatus;

	OGRSpatialReference map_srs;

	Geom3DType geom_3d_type;

	bool cull_face = false;
private:
	//���Ʋ�������
	PaintType ptype;

	//��ʾģʽ
	DisplayType dttype;

	//�Ƿ����
	FillStatus FillStatus;

	//����ɫ��
	Color GridColor;
};

extern Status status;
extern int winHeight, winWidth;
extern void (*set_unit)(int, int, Color);
extern Color(*get_unit)(int, int);
extern vector<SeedPoint> seed_points;
extern void(*FillFunction)(Point2D, Color);
extern double Kd, Ks, Ka;
extern int shininess;