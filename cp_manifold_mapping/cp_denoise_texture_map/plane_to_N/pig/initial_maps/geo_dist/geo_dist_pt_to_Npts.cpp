/* read mesh from file and 
 - if one vertex is specified, for all vertices of the mesh print their distances to this vertex
 - if two vertices are specified, print the shortest path between these vertices 

	Danil Kirsanov, 01/2008 
*/
#include <iostream>
#include <fstream>
#include "geodesic_algorithm_exact.h"


int main(int argc, char **argv) 
{
	if(argc < 2)
	{
		std::cout << "usage: mesh_file_name [first_vertex] [second_vertex]" << std::endl;
		return 0;
	}

	std::vector<double> points;	
	std::vector<unsigned> faces;

	bool success = geodesic::read_mesh_from_file(argv[1],points,faces);
	if(!success)
	{
		std::cout << "something is wrong with the input file" << std::endl;
		return 0;
	}

	geodesic::Mesh mesh;
	mesh.initialize_mesh_data(points, faces);		//create internal mesh data structure including edges

	geodesic::GeodesicAlgorithmExact algorithm(&mesh);	//create exact algorithm for the mesh

	unsigned source_vertex_index = (argc == 2) ? 0 : atol(argv[2]);

	geodesic::SurfacePoint source(&mesh.vertices()[source_vertex_index]);		//create source 
	std::vector<geodesic::SurfacePoint> all_sources(1,source);					//in general, there could be multiple sources, but now we have only one
		algorithm.propagate(all_sources);	//cover the whole mesh
		//algorithm.print_statistics();

		for(unsigned i=0; i<mesh.vertices().size(); ++i)
		{
			geodesic::SurfacePoint p(&mesh.vertices()[i]);		

			double distance;
			unsigned best_source = algorithm.best_source(p,distance);		//for a given surface point, find closets source and distance to this source

			std::cout << distance << " ";		//print geodesic distance for every vertex
		}
  		std::cout << std::endl;

	return 0;
}	
