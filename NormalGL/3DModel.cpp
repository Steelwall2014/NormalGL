#include "3DModel.h"
#include <fstream>
#include "Utils.h"
#include "Tin.h"
#include "GeoDefine.h"
#include "Matrix.h"

double SCALE_FACTOR = 150.0;

Model::Model(string filepath)
{
	this->vertice = new Tin();
	ifstream file = ifstream(filepath);
	//char line[50];
	string str_line;
	//string fff = "f";
	//vector<Point3D> vertice;
	while (!file.eof())
	{
		getline(file, str_line);
		//str_line = string(line);
		if (str_line.size() < 5)
			continue;
		vector<string> parts = split(str_line, ' ');
		if (parts[0] == "v")
		{
			double x = stod(parts[1]) * SCALE_FACTOR, z = stod(parts[2]) * SCALE_FACTOR, y = stod(parts[3]) * SCALE_FACTOR;
			Point3D point = Point3D{ x, y, z };
			vertice->push(point);
			vertex_norm_vecs.push_back(Point3D{ 0, 0, 0 });
			vertice->SetLastElemId();
			//(*vertice)[vertice->size() - 1].id = vertice->size() - 1;
		}
		else if (parts[0] == "f")
		{
			int index1 = stoi(parts[1])-1, index2 = stoi(parts[2])-1, index3 = stoi(parts[3])-1;
			Face tri = Face(vertice, index1, index2, index3);
			Point3D face_norm = GetPlaneNormalVector_Point3D((*vertice)[index1], (*vertice)[index2], (*vertice)[index3]);
			faces.push_back(tri);
			faces_norm_vecs.push_back(face_norm);
			faces[faces.size() - 1].face_id = faces.size() - 1;
			_add_vertex_norm(index1, face_norm);
			_add_vertex_norm(index2, face_norm);
			_add_vertex_norm(index3, face_norm);
		}
	}
	for (int i = 0; i < vertex_norm_vecs.size(); i++)
		vertex_norm_vecs[i].Normalize();
	for (int i = 0; i < faces.size(); i++)
		faces[i].face_id = i;
	//vertice = Tin(vertice);0x010ff9f4
	vertex_IR = vector<double>(vertice->size());
	vertex_IG = vector<double>(vertice->size());
	vertex_IB = vector<double>(vertice->size());
}

void Model::_add_vertex_norm(int index, Matrix &norm)
{
	vertex_norm_vecs[index].x += norm.data[0][0];
	vertex_norm_vecs[index].y += norm.data[0][1];
	vertex_norm_vecs[index].z += norm.data[0][2];
}
void Model::_add_vertex_norm(int index, Point3D norm)
{
	vertex_norm_vecs[index].x += norm.x;
	vertex_norm_vecs[index].y += norm.y;
	vertex_norm_vecs[index].z += norm.z;
}
