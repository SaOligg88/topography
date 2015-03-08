#include <iostream>
#include <fstream>
#include "geodesic_algorithm_exact.h"
#include <sstream>
using namespace std;

int main(int argc, char **argv) 
{
	std::vector<double> points;	
	std::vector<unsigned> faces;
	bool success = geodesic::read_mesh_from_file(argv[1],points,faces);
	geodesic::Mesh mesh;
	mesh.initialize_mesh_data(points, faces);
	geodesic::GeodesicAlgorithmExact algorithm(&mesh);	
	const char *path="/scr/litauen2/projects/distance/condor/";
  	std::ostringstream fileNameStream("");
    		fileNameStream << path << argv[2] << ".txt";
		std::string fileName = fileNameStream.str();
	ofstream outputFile;
	outputFile.open(fileName.c_str());
	
	unsigned source_vertex_index = atol(argv[2]);	
	geodesic::SurfacePoint source(&mesh.vertices()[source_vertex_index]);
	std::vector<geodesic::SurfacePoint> all_sources(1,source);
	algorithm.propagate(all_sources);	
	for(unsigned i=source_vertex_index; i<mesh.vertices().size(); ++i) 
	{
		geodesic::SurfacePoint p(&mesh.vertices()[i]);		
		double distance;
		unsigned best_source = algorithm.best_source(p,distance);
		outputFile << distance << std::endl;	
		std::cout << distance << std::endl;
	}
	outputFile.close();
	return 0;
}	
