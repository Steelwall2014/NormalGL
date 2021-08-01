#pragma once
//#include "Graphic.h"
#include "Utils.h"
#include "GeoDefine.h"
#include "Painters.h"
#include "Clipper.h"
class DataManager;
class MapRenderer
{
public:
	MapRenderer();
	~MapRenderer();
	MapRenderer(DataManager *dm);
	//void AddLayer(Layer *);
	void AddPainterFromGeom(vector<OGRGeometry *> &);//, vector<Painter *> &, Matrix &);
	// 将图层缩放到适应地图画布
	void AdjustProjDatasetToMap();
	void AddPainter(Painter *);
	vector<Painter *> &GetDrawnPainters();
	void Reset();
	void Render();
	bool NothingToRender();		// 判断是否有需要渲染的painter
	void LinkDataManager(DataManager *dm);	// 和datamanager建立联系
	Clipper clipper;
	Matrix imported_painters_transformation = Matrix::CreateIdentityMatrix(3);	// 矢量文件的平移旋转缩放变换矩阵
	Matrix geom_to_map_transformation = Matrix::CreateIdentityMatrix(3);		// 手绘图形的平移旋转缩放变换矩阵
private:
	//OGREnvelope dataset_extent;
	//OGREnvelope map_extent;
	//Layer *initial_layer, *proj_layer;
	//Dataset *initial_dataset, *projected_dataset;
	DataManager *dm;
	vector<Painter *> imported_geoms, drawn_geoms;
};


//class Renderer3D
//{
//public:
//	Renderer3D();
//	~Renderer3D();
//	void Reset();
//	void Render();
//};



//{x = 479.00000000000000 y = 393.00000000000000 }
//{x = 394.00000000000000 y = 144.00000000000000 }
//{x = 691.00000000000000 y = 80.000000000000000 }
//{x = 708.00000000000000 y = 436.00000000000000 }
//{x = 393.00000000000000 y = 460.00000000000000 }
//{x = 479.00000000000000 y = 393.00000000000000 }