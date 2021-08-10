#ifndef MESHPROCESSING_MESHDEFINITION_H
#define MESHPROCESSING_MESHDEFINITION_H

#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
typedef OpenVolumeMesh::GeometryKernel<OpenVolumeMesh::Geometry::Vec3d> Src_VolumeMesh;

class VolumeMesh : public Src_VolumeMesh
{
public:
	VolumeMesh() :Src_VolumeMesh()
	{
	}

	OpenVolumeMesh::Geometry::Vec3d BBox_min;
	OpenVolumeMesh::Geometry::Vec3d BBox_max;
	OpenVolumeMesh::Geometry::Vec3d BBox_center;
	double avg_edge_length;
	std::vector<double> volume;
};

#endif