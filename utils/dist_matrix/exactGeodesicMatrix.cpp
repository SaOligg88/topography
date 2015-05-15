#include <iostream>
#include <fstream>
#include "geodesic_algorithm_exact.h"
#include <sstream>
using namespace std;

// To compile: 
// g++ -Wall -I../geodesic exactGeodesicMatrix.cpp -o exactGeodesicMatrix

// Step 1: Convert surface to ascii
// mris_convert [surface] [surface_tmp.asc]
// Remove first line with:
// tail -n +2 [surface_tmp.asc] > [surface_tmp2.asc]
// sed -e 's/\( 0\)*$//g' [surface_tmp2.asc] > [surface.asc]
// rm -f [surface_tmp.asc] [surface_tmp2.asc]

// Step 2: 
// ./exactGeodesicMatrix [surface.asc] [/output directory/] [output filename prefix]

// Step 3:
// Reassemble full matrix and delete output vectors
// [script]

// BEGIN //
int main(int argc, char **argv) 
{
	if(argc < 6)
	{
		std::cout << "usage: ./exactGeodesicMatrix [surface.asc] [/output directory/] [output filename prefix] [index file] [node to run]" << std::endl; 
		return 0;
	}	
	std::vector<double> points;	
	std::vector<unsigned> faces;
	bool success = geodesic::read_mesh_from_file(argv[1],points,faces);
	geodesic::Mesh mesh;
	mesh.initialize_mesh_data(points, faces);
	geodesic::GeodesicAlgorithmExact algorithm(&mesh);
    	std::ostringstream fileNameStream("");
	fileNameStream << argv[2] << argv[3] << argv[5] << ".txt";
	std::string fileName = fileNameStream.str();	
	ofstream outputFile;
	outputFile.open(fileName.c_str());
	
	// Read index file:
	std::vector<unsigned> indices;
	ifstream indexfile(argv[4]);
	unsigned line;
	while ( indexfile >> line ) { //std::getline(indexfile,line)) {
		indices.push_back(line);
	}
	indexfile.close();

	unsigned source_vertex_index = atol(argv[5]);	
	geodesic::SurfacePoint source(&mesh.vertices()[indices[source_vertex_index]]);
	std::vector<geodesic::SurfacePoint> all_sources(1,source);
	algorithm.propagate(all_sources);

	// get distance:
	// change source_vertex_index to 0 to get all nodes for each run
	for(unsigned i=source_vertex_index; i<indices.size(); ++i) 
	{
		geodesic::SurfacePoint p(&mesh.vertices()[indices[i]]);		
		double distance;
		unsigned best_source = algorithm.best_source(p,distance);
		outputFile << distance << std::endl;			
	}
	outputFile.close();	
	return 0;
}

