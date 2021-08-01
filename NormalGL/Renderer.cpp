#include "Renderer.h"
#include "Status.h"
#include "GeoDefine.h"
#include "Graphic.h"
#include "ogr_spatialref.h"
#include "gdal.h"
#include "ogrsf_frmts.h"
#include "ogr_geometry.h"
#include "DataManagement.h"

//#include "..\\3rdParty\\gdal\\include\\ogr_spatialref.h"
//#include "..\\3rdParty\\gdal\\include\\gdal.h"
//#include "..\\3rdParty\\gdal\\include\\ogrsf_frmts.h"
//#include "..\\3rdParty\\gdal\\include\\ogr_geometry.h"

// 注意：如果是在vs2015之前的进行编译，需要把文件cpl_config.h第20行的#define snprintf _snprintf取消注释
extern Status status;
extern int winHeight, winWidth;
extern vector<SeedPoint> seed_points;
MapRenderer::MapRenderer()
{
	geom_to_map_transformation = Matrix::CreateIdentityMatrix(3);
}

MapRenderer::MapRenderer(DataManager *dm)
{
	this->LinkDataManager(dm);
	geom_to_map_transformation = Matrix::CreateIdentityMatrix(3);
}

MapRenderer::~MapRenderer()
{
	Reset();
}

void MapRenderer::Render()
{
	//if (projected_dataset.getLayerCount() == 0)
	//	return;
	winHeight = getWindowHeight(); winWidth = getWindowWidth();
	if (status.GetDisplayType() == dtPixelMode)
	{
		for (auto painter : imported_geoms)
		{
			painter->Paint(this->imported_painters_transformation);
		}
	}

	int drawn_geoms_count = drawn_geoms.size();
	
	if (!clipper.GetBoundary().empty())
	{
		vector<Painter *> to_be_deleted;
		for (int i = 0; i < drawn_geoms_count; i++)
		{
			Painter *painter = drawn_geoms[i];
			string painter_type = painter->GetPaintShape();
			switch (status.ctype)
			{
			case ctCenterClip:
				if (painter_type == "POLYLINE" || painter_type == "LINE")
				{
					this->clipper.CenterLineClip(painter, drawn_geoms);
					to_be_deleted.push_back(painter);
				}
				break;
			case ctSHClip:
				if (painter_type == "POLYGON")
				{
					this->clipper.SutherlandHodgmanClip(painter, painter);
				}
				break;
			case ctWAClip:
				if (painter_type == "POLYGON")
				{
					int flag = this->clipper.WeilerAthertonClip(painter, drawn_geoms);
					if (flag)
						to_be_deleted.push_back(painter);
				}
				break;
			case ctClip:
				if (painter_type == "POLYLINE" || painter_type == "LINE")
				{
					this->clipper.CenterLineClip(painter, drawn_geoms);
					to_be_deleted.push_back(painter);
				}
				else if (painter_type == "POLYGON")
				{
					this->clipper.WeilerAthertonClip(painter, drawn_geoms);
					to_be_deleted.push_back(painter);
				}
				break;
			}
		}
		for (auto painter : to_be_deleted)
		{
			drawn_geoms.erase(remove(drawn_geoms.begin(), drawn_geoms.end(), painter), 
							  drawn_geoms.end());
		}
	}

	for (auto painter : drawn_geoms)
	{
		painter->Paint();
	}
	for (auto &seed_point : seed_points)
		FillFunction(seed_point.point, seed_point.seed_color);
	clipper.Reset();
}

void MapRenderer::Reset()
{
	this->geom_to_map_transformation = Matrix::CreateIdentityMatrix(3);
	//this->dm = nullptr;
	//this->projected_dataset = nullptr;
	for (int i = 0; i < this->imported_geoms.size(); i++)
		delete this->imported_geoms[i];
	this->imported_geoms.clear();
	for (int i = 0; i < this->drawn_geoms.size(); i++)
		delete this->drawn_geoms[i];
	this->drawn_geoms.clear();
}

//void MapRenderer::AddLayer(Layer *pLayer)
//{
//	dm->initial_dataset.addLayer(pLayer);
//	OGRCoordinateTransformation *coord_trans = OGRCreateCoordinateTransformation(pLayer->layer_srs, &status.map_srs);
//	Layer *projected_layer = new Layer(pLayer->geomType);
//	for (auto geom : pLayer->geometrySet)
//	{
//		OGRGeometry *projected_geom = geom->clone();
//		projected_geom->transform(coord_trans);
//		projected_layer->addGeometry(projected_geom, true);
//	}
//	dm->projected_dataset.addLayer(projected_layer);
//	if (dm->projected_dataset->getLayerCount() == 1)
//		AdjustProjDatasetToMap();
//	AddPainterFromGeom(projected_layer->geometrySet);//, this->imported_geoms, this->geom_to_map_transformation);
//	//map_extent.Merge(projected_layer->envelop);
//}

void MapRenderer::AddPainter(Painter *painter)
{
	this->drawn_geoms.push_back(painter);
}

void MapRenderer::AddPainterFromGeom(vector<OGRGeometry *> &geom_set)//, vector<Painter *> &painter_set, Matrix &geom_to_map_trans)
{
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
			after_trans_points = WindowCoordinateTransform(points, point_count, geom_to_map_transformation);
			Painter *painter = new PolylinePainter(after_trans_points, point_count);
			this->imported_geoms.push_back(painter);
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
				after_trans_points = WindowCoordinateTransform(points, point_count, geom_to_map_transformation);
				this->imported_geoms.push_back(new PolylinePainter(after_trans_points, point_count));
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
			after_trans_points = WindowCoordinateTransform(points, point_count, geom_to_map_transformation);
			this->imported_geoms.push_back(new PolygonPainter(after_trans_points, point_count-1));
			delete[] points;
			delete[] after_trans_points;
			//delete ring;
			break;
		}
		}
	}
}

bool MapRenderer::NothingToRender()
{
	return this->drawn_geoms.empty() && this->imported_geoms.empty();
}

vector<Painter *> &MapRenderer::GetDrawnPainters()
{
	return this->drawn_geoms;
}

void MapRenderer::AdjustProjDatasetToMap()
{
	winHeight = getWindowHeight(); winWidth = getWindowWidth();
	OGREnvelope env = dm->projected_dataset.dataset_envelope;
	double dataset_envcenter_x = (env.MaxX + env.MinX) / 2.0, dataset_envcenter_y = (env.MaxY + env.MinY) / 2.0;
	double mapUnitsPerPixelY = (env.MaxY - env.MinY) / (double)winHeight;
	double mapUnitsPerPixelX = (env.MaxX - env.MinX) / (double)winWidth * 0.95;
	double scale = mapUnitsPerPixelY > mapUnitsPerPixelX ? mapUnitsPerPixelY : mapUnitsPerPixelX;

	geom_to_map_transformation = 
		Matrix::CreateTranslationMatrix(-dataset_envcenter_x, -dataset_envcenter_y) * 
		Matrix::CreateZoomMatrix(1.0/scale*0.95) * 
		Matrix::CreateTranslationMatrix((double)winWidth/2.0, (double)winHeight/2.0);
}

void MapRenderer::LinkDataManager(DataManager *dm)
{
	this->dm = dm;
	dm->renderer = this;
}