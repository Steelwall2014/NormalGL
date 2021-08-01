#pragma once
#include <vector>
#include "Tin.h"
#include "3DModel.h"
class MapRenderer;
class Model;
//class Tin;
class DataManager
{
private:
	
	//vector<Painter *> imported_geoms, drawn_geoms;
public:
	struct Pyramid
	{
		Tin *tin;
		Face pyramid_f1, pyramid_f2, pyramid_f3, pyramid_f4;
	};
	Line3D x_axis, y_axis, z_axis;
	Point3D origin = Point3D(0, 0, 0);
	Pyramid pyramid;
	Model teapot, junc;
	DataManager();
	MapRenderer *renderer;
	Dataset initial_dataset, projected_dataset;
	void AddLayer(Layer *pLayer);
	void SendToRender(MapRenderer *renderer);
	void Clear();
};

