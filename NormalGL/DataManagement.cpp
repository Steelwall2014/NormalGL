#include "DataManagement.h"
#include "Painters.h"
#include <windows.h>
#include "Graphic.h"
#include "Status.h"
#include "ogr_spatialref.h"
#include "gdal.h"
#include "ogrsf_frmts.h"
#include "ogr_geometry.h"
#include "Renderer.h"
//#include "Tin.h"
DataManager::DataManager()
{
	x_axis = Line3D(origin, Point3D(500, 0, 0));
	y_axis = Line3D(origin, Point3D(0, 500, 0));
	z_axis = Line3D(origin, Point3D(0, 0, 500));
	teapot = Model(".\\resource\\teapot.obj");
	junc = Model(".\\resource\\junc.obj");
	vector<Point3D> points;
	points.push_back(Point3D{ 250, -250, -250 });	// 0
	points.push_back(Point3D{ -250, -250, -250 });	// 1
	points.push_back(Point3D{ 0, 0, 250 });			// 2
	points.push_back(Point3D{ -250, 250, -250 });	// 3
	points.push_back(Point3D{ 250, 250, -250 });	// 4
	pyramid.tin = new Tin(points);
	
	pyramid.pyramid_f1 = Face(pyramid.tin, 0, 1, 2);
	pyramid.pyramid_f2 = Face(pyramid.tin, 1, 3, 2);
	pyramid.pyramid_f3 = Face(pyramid.tin, 3, 4, 2);
	pyramid.pyramid_f4 = Face(pyramid.tin, 4, 0, 2);

	//points.push_back(Point3D{ 250, -250, -250 });
	//points.push_back(Point3D{ -250, -250, -250 });
	//points.push_back(Point3D{ 0, 0, 250 });
	//pyramid_f1 = Face(points);
	//points[0] = Point3D{ -250, -250, -250 };
	//points[1] = Point3D{ -250, 250, -250 };
	//points[2] = Point3D{ 0, 0, 250 };
	//pyramid_f2 = Face(points);
	//points[0] = Point3D{ -250, 250, -250 };
	//points[1] = Point3D{ 250, 250, -250 };
	//points[2] = Point3D{ 0, 0, 250 };
	//pyramid_f3 = Face(points);
	//points[0] = Point3D{ 250, 250, -250 };
	//points[1] = Point3D{ 250, -250, -250 };
	//points[2] = Point3D{ 0, 0, 250 };
	//pyramid_f4 = Face(points);
}


void DataManager::AddLayer(Layer *pLayer)
{
	initial_dataset.addLayer(pLayer);
	OGRCoordinateTransformation *coord_trans = OGRCreateCoordinateTransformation(pLayer->layer_srs, &status.map_srs);
	Layer *projected_layer = new Layer(pLayer->geomType);
	for (auto geom : pLayer->geometrySet)
	{
		OGRGeometry *projected_geom = geom->clone();
		if (coord_trans != nullptr)
			projected_geom->transform(coord_trans);
		projected_layer->addGeometry(projected_geom, true);
	}
	projected_dataset.addLayer(projected_layer);
	if (projected_dataset.getLayerCount() == 1)
		this->renderer->AdjustProjDatasetToMap();
	//vector<Painter *> temp_painters = FromOGRGeometriesToPainters(projected_layer->geometrySet, this->renderer->geom_to_map_transformation);
	this->renderer->AddPainterFromGeom(projected_layer->geometrySet);//, this->imported_geoms, this->geom_to_map_transformation);
}
void DataManager::SendToRender(MapRenderer *renderer)
{
	renderer->LinkDataManager(this);
	this->renderer = renderer;
}

void DataManager::Clear()
{
	initial_dataset.reset();
	projected_dataset.reset();
}
