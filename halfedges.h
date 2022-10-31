/**********************************************************************
 * flux: framework for learning unstructured meshing
 * Copyright (c) 2022 Philip Caplan. All rights reserved.
 * Licensed under the MIT License (https://mit-license.org/)
 **********************************************************************/
#ifndef FLUX_HALFEDGES_H_
#define FLUX_HALFEDGES_H_

#include "mat.h"
#include "vec.h"

#include <memory>
#include <set>
#include <vector>

namespace flux {

// forward declaration of a mesh
template<typename type> class Mesh;

struct HalfEdge;

/**
 * \brief Represents a face in a half-edge based mesh.
 */
struct HalfFace {
  HalfEdge* edge; // pointer to one edge which has the face on the left
  int index = -1; // index of the face (i.e. element number)
  vec4d p;        // [nx, ny, nz, -dot(x,n)] where x is a point on the face
                  // n is the face normal (used in simplification algorithm)
};

/**
 * \brief Represents a vertex in a half-edge based mesh.
 */
struct HalfVertex {
  /**
   * \brief Constructs a half vertex from coordinates.
   *
   * \param[in] dim number of coordinates to copy (allows copying 2d or 3d coordiantes)
   * \param[in] x - pointer to coordinates
   */
  HalfVertex( int dim , const double* x ) {
    flux_assert( dim == 2 || dim == 3 );
    for (int d = 0; d < dim; d++)
      point[d] = x[d];
  }
  HalfEdge* edge; // pointer to one edge emanating from this vertex
  vec3d point;    // the coordinates of this vertex (always 3d!)
  int index = -1; // index of this vertex (i.e. vertex number)
  mat44d Q;       // error quadric (used in simplification algorithm)
};

/**
 * \brief Represents an edge in a half-edge based mesh.
 */
struct HalfEdge {
  HalfVertex* vertex; // pointer to the 'origin' vertex of this edge
  HalfEdge*   twin;   // pointer to the 'opposite' edge (i.e. parallel but in opposite direction)
  HalfEdge*   next;   // pointer to the next edge around the face in CCW order
  HalfEdge*   prev;   // pointer to the previous edge around the face in CCW order
  HalfFace*   face;   // pointer to the left face of this oriented edge

  // data used in simplification algorithm
  vec4d  vbar; // optimal point (in homogeneous coordinates) where the endpoints should be collapsed
  mat44d Qbar; // sum of endpoint Q's, i.e. Qbar = Q1 + Q2
  double cost;  // the cost (error) of this edge: vbar^T * Qbar * vbar
};

// the following operators are needed to find entities so they can be deleted from the sets
inline bool operator< ( const std::unique_ptr<HalfEdge>& lhs , const HalfEdge* rhs ) {
  return std::less<const HalfEdge*>()(lhs.get(),rhs);
}

inline bool operator< (const HalfEdge* lhs , const std::unique_ptr<HalfEdge>& rhs ) {
  return std::less<const HalfEdge*>()(lhs,rhs.get());
}

inline bool operator< ( const std::unique_ptr<HalfVertex>& lhs , const HalfVertex* rhs ) {
  return std::less<const HalfVertex*>()(lhs.get(),rhs);
}

inline bool operator< (const HalfVertex* lhs , const std::unique_ptr<HalfVertex>& rhs ) {
  return std::less<const HalfVertex*>()(lhs,rhs.get());
}

inline bool operator< ( const std::unique_ptr<HalfFace>& lhs , const HalfFace* rhs ) {
  return std::less<const HalfFace*>()(lhs.get(),rhs);
}

inline bool operator< (const HalfFace* lhs , const std::unique_ptr<HalfFace>& rhs ) {
  return std::less<const HalfFace*>()(lhs,rhs.get());
}

/**
 * \brief Container for half-edge entities (vertices, edges, faces)
 *        Can be used for Triangle or Quad meshes.
 */
template<typename type>
class HalfEdgeMesh {
public:

  // edge iterators
  typedef std::set< std::unique_ptr<HalfEdge> >::iterator edge_itr;
  typedef std::set< std::unique_ptr<HalfEdge> >::const_iterator const_edge_itr;

  // face iterators
  typedef std::set< std::unique_ptr<HalfFace> >::iterator face_itr;
  typedef std::set< std::unique_ptr<HalfFace> >::const_iterator const_face_itr;

  // vertex iterators
  typedef std::set< std::unique_ptr<HalfVertex> >::iterator vertex_itr;
  typedef std::set< std::unique_ptr<HalfVertex> >::const_iterator const_vertex_itr;

  /**
   * \brief Constructs all the half-edge entities from an incoming connectivity-based representation.
   *
   * \param[in] the mesh represented by vertices and element connectivity.
   */
  HalfEdgeMesh( const Mesh<type>& mesh );

  /**
   * \brief Converts this half-edge based mesh representation to vertices and element connectivity.
   *
   * \param[out] mesh - the mesh to be filled with element connectivity
   */
  void extract( Mesh<type>& mesh ) const;

  /**
   * \brief Computes the vertices surrounding a vertex (v) in CW order.
   */
  void get_onering( const HalfVertex* v , std::vector<HalfVertex*>& onering ) const;

  /**
   * \brief Computes the edges surrounding a vertex (v) in CW order.
   */
  void get_onering( const HalfVertex* v , std::vector<HalfEdge*>& onering ) const;

  /**
   * \brief Computes the faces surrounding a vertex (v) in CW order.
   */
  void get_onering( const HalfVertex* v , std::vector<HalfFace*>& onering ) const;

  /**
   * \brief Returns the set of HalfVertex pointers.
   */
  std::set<std::unique_ptr<HalfVertex>, std::less<> >& vertices() { return vertices_; }

  /**
   * \brief Returns the set of HalfEdge pointers.
   */
  std::set<std::unique_ptr<HalfEdge> , std::less<> >& edges() { return edges_; }

  /**
   * \brief Returns the set of HalfFace pointers.
   */
  std::set<std::unique_ptr<HalfFace>, std::less<> >& faces() { return faces_; }

  /**
   * \brief Creates a HalfEdge, storing it and returning a pointer to the created edge.
   */
  HalfEdge* create_edge() {
    auto re = edges_.insert( std::make_unique<HalfEdge>() );
    flux_assert( re.second );
    return re.first->get();
  }

  /**
   * \brief Creates a HalfVertex, storing it and returning a pointer to the created vertex.
   */
  HalfVertex* create_vertex( int dim , const double* x ) {
    auto rv = vertices_.insert( std::make_unique<HalfVertex>(dim,x) );
    flux_assert( rv.second );
    return rv.first->get();
  }

  /**
   * \brief Creates a HalfFace, storing it and returning a pointer to the created face.
   */
  HalfFace* create_face() {
    auto rf = faces_.insert( std::make_unique<HalfFace>() );
    flux_assert( rf.second );
    return rf.first->get();
  }

  /**
   * \brief Removes a HalfEdge from this container (useful when modifying meshes).
   */
  void remove( HalfEdge* edge ) {
    edge_itr it = edges_.find(edge);
    flux_assert( it != edges_.end() );
    edges_.erase(it);
  }

  /**
   * \brief Removes a HalfFace from this container (useful when modifying meshes).
   */
  void remove( HalfFace* face ) {
    face_itr it = faces_.find(face);
    flux_assert( it != faces_.end() );
    faces_.erase(it);
  }

  /**
   * \brief Removes a HalfVertex from this container (useful when modifying meshes).
   */
  void remove( HalfVertex* vertex ) {
    vertex_itr it = vertices_.find(vertex);
    flux_assert( it != vertices_.end() );
    vertices_.erase(it);
  }

  /**
   * \brief A helper function to retrieve a HalfVertex from an iterator.
   */
  HalfVertex* get( vertex_itr& x ) { return x->get(); }

  /**
   * \brief A helper function to retrieve a HalfEdge from an iterator.
   */
  HalfEdge* get( edge_itr& x ) { return x->get(); }

  /**
   * \brief A helper function to retrieve a HalfFace from an iterator.
   */
  HalfFace* get( face_itr& x ) { return x->get(); }

  /**
   * \brief Performs a few checks on the vertices and faces to make sure everything
   *        was constructed properly. No error is raised, but returns whether everything
   *        is constructed correctly or not.
   */
  bool check() const;

  /**
   * \brief Retrieves a single boundary edge.
   */
  HalfEdge* get_boundary_edge() const;

  /**
   * \brief Returns the current number of boundary edges in the mesh with an option to recompute them.
   *
   * \param[count] option to re-count (default = false)
   */
  int nb_boundary( bool count = false ) const {
    if (count) count_boundary_edges();
    return nb_boundary_;
  }

  /**
   * \brief Returns the number of connected boundaries in the mesh.
   *        If the mesh is not closed, there should be at least 1 outer boundary
   *        and any number of holes in the interior of the mesh.
   *        Returns 0 if the mesh is closed.
   */
  int nb_connected_boundaries() const;

  /**
   * \brief Flips a half-edge, updating all affected HalfEdge, HalfVertex and HalfFace's.
   *        Assumes a closed mesh (no boundary).
   *
   * \param[in] edge - pointer to the HalfEdge to be flipped
   */
 void flip( HalfEdge* edge );

private:
  /**
   * \brief Builds all the half-edge entities from the connectivity-based mesh
   */
  void build();

  const Mesh<type>& mesh_; // reference to the original mesh from which the entities were constructed

  // the half-edge vertices, edges and faces
  std::set< std::unique_ptr<HalfEdge> , std::less<> >   edges_;
  std::set< std::unique_ptr<HalfVertex> , std::less<> > vertices_;
  std::set< std::unique_ptr<HalfFace> , std::less<> >   faces_;

  /**
   * \brief Counts the number of boundary edges, updating nb_boundary_
   */
  void count_boundary_edges() const;
  mutable int nb_boundary_;
};

} // flux

#endif
