#include "LayerImporter.h"
#include "GeoDefine.h"

#include "ogrsf_frmts.h"
#include "ogr_feature.h"
#include "gdal.h"
#include "ogr_api.h"
#include "ogr_spatialref.h"

//#include "..\\3rdParty\\gdal\\include\\ogrsf_frmts.h"
//#include "..\\3rdParty\\gdal\\include\\ogr_feature.h"
//#include "..\\3rdParty\\gdal\\include\\gdal.h"
//#include "..\\3rdParty\\gdal\\include\\ogr_api.h"
//#include "..\\3rdParty\\gdal\\include\\ogr_spatialref.h"
#include <assert.h>


#pragma comment(lib,".\\3rdParty\\gdal-3.1.4\\lib\\gdal_i.lib")

Layer* LayerImporter::importShpLayer(const char* path , const char* labelFieldName )
{
	GDALAllRegister();
	OGRRegisterAll();
	CPLSetConfigOption("GDAL_DATA", ".\\3rdParty\\gdal-3.1.4\\data");
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
	CPLSetConfigOption("SHAPE_ENCODING", "");

	//OGRSFDriver* poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
	//OGRDataSource* poDS = poDriver->Open( path );

	GDALDataset *poDS = (GDALDataset *)GDALOpenEx(path, GDAL_OF_READONLY, NULL, NULL, NULL);
	
	if( poDS == NULL ) return NULL;

	OGRLayer* poLayer = poDS->GetLayer(0);
	if( poLayer == NULL ) return NULL;

	OGREnvelope ev;
	poLayer->GetExtent(&ev);
	Layer* pLayer = new Layer( poLayer->GetGeomType());
	poLayer->GetExtent(&pLayer->envelop);
	//pLayer->setEnvelop( ev.MinX, ev.MinY, ev. MaxX, ev.MaxY );

	int labelFieldIndex = -1;
	if( strcmp( labelFieldName, "") != 0 )
		labelFieldIndex = poLayer->GetLayerDefn()->GetFieldIndex( labelFieldName );
	OGRFeature * pFeature;
	while ((pFeature = poLayer->GetNextFeature()) )
	{
		OGRGeometry* pOGRGeometry = pFeature->GetGeometryRef()->clone();
		pLayer->addGeometry(pOGRGeometry);
		//switch (pOGRGeometry->getGeometryType())
		//{
		//case wkbPoint:
		//	{
		//		OGRPoint* pOGRPoint = (OGRPoint*)pOGRGeometry;
		//		pLayer->addGeometry(pOGRPoint, false);
		//	}
		//	break;
		//case wkbLineString:
		//	{
		//		OGRLineString* pOGRLineString = (OGRLineString*)pOGRGeometry;
		//		PolylineGeometry* pGeometry = new PolylineGeometry( ) ;
		//		for (int i = 0, size = pOGRLineString->getNumPoints(); i < size ; ++i)
		//		{
		//			pGeometry->addPoint( pOGRLineString->getX(i),pOGRLineString->getY(i) );
		//		}
		//		if( labelFieldIndex >= 0 )
		//			pGeometry->label = pFeature->GetFieldAsString( labelFieldIndex );
		//		pLayer->addGeometry( pGeometry , false);
		//	}
		//	break;
		//case wkbPolygon:
		//	{
		//		OGRPolygon* pOGRPolygon = (OGRPolygon*)pOGRGeometry;
		//		PolygonGeometry* pGeometry = new PolygonGeometry( ) ;
		//		OGRLinearRing* ring = pOGRPolygon->getExteriorRing();
		//		for (int i = 0, size = ring->getNumPoints(); i < size ; ++i)
		//		{
		//			//assert(ring->getX(i) < 120 && ring->getY(i) < 32);
		//			pGeometry->addPoint( ring->getX(i),ring->getY(i) );
		//		}		
		//		if( labelFieldIndex >= 0 )
		//			pGeometry->label = pFeature->GetFieldAsString( labelFieldIndex );
		//		pLayer->addGeometry( pGeometry , false);
		//	}
		//	break;
		//}
		OGRFeature::DestroyFeature(pFeature);

		//static int ii = 0;
		//++ii;
		//if( ii > 0 ) break;
	}

	//OGRDataSource::DestroyDataSource(poDS);
	//OGRCleanupAll();
	GDALClose(poDS);
	return pLayer;
}
