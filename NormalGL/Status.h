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

//绘制操作类型
enum PaintType { ptNone = -999, ptDrawRaster = -1, ptDrawLine = 0, ptDrawPolyline, ptDrawPolygon, ptDrawCircle, ptDrawEllipse, ptDrawRect };

//显示模式
enum DisplayType { dtPixelMode = 1, dtGridMode_5 = 5, dtGridMode_10 = 10, dtGridMode_20 = 20 };

enum FillStatus { fsNone = 0, fsAutoAET, fsRecurse, fsStack, fsScanline };

enum OperationStatus { osNone, osTransform, osClip, os3D };// tsTranslation, tsRotate, tsZoom };

enum ClipType { ctCenterClip, ctSHClip, ctWAClip, ctClip };// 剪裁算法

enum Geom3DType { gtTeapot, gtPyramid,gtJunc }; // 三维模型有哪些
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

	//设置绘制类型比如线、矩形等
	void SetPaintType(PaintType pType)
	{
		ptype = pType;
	}

	//设置显示模式（像素、棋盘格）
	void SetDisplayType(DisplayType dtType)
	{
		dttype = dtType;
	}

	//设置是否填充
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
	//绘制操作类型
	PaintType ptype;

	//显示模式
	DisplayType dttype;

	//是否填充
	FillStatus FillStatus;

	//格网色彩
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