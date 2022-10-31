#ifndef FLUX_PROJECT4_LOOP_H_
#define FLUX_PROJECT4_LOOP_H_

namespace flux {

/**
* \brief Gets the point opposite to an edge, given a triangle.
*/
int 
getOppositePoint( int p , int q , int t, std::vector<int> triangles );

/**
* \brief Creates a data structure of edge indeces and their associated triangles.
*/
int** 
edge_to_triangle( const Mesh<Triangle>& mesh );

/**
* \brief Creates a data structure of triangle indeces and their associated edges.
*/
int** 
triangle_to_edge( const Mesh<Triangle>& mesh );

/**
* \brief Applies Charles Loop's subdivision for triangle meshes in 2d or 3d.
*/
void
loop_subdivision( Mesh<Triangle>& mesh, int levels, int dim );

} // flux

#endif
