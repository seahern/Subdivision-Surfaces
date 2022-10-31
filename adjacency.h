/**********************************************************************
 * flux: framework for learning unstructured meshing
 * Copyright (c) 2022 Philip Caplan. All rights reserved.
 * Licensed under the MIT License (https://mit-license.org/)
 **********************************************************************/
#ifndef FLUX_ADJACENCY_H_
#define FLUX_ADJACENCY_H_

#include "array2d.h"

namespace flux {

// forward declaration of Mesh
template<typename type> class Mesh;

/**
 * \brief Represents the mesh adjacency information in a mesh,
 *        i.e. the element-to-element relations.
 *        Neighbours are stored according to the ordering of the faces.
 *        For example, adj(k,j) returns the neighbor of element k
 *        along the j-th face in the canonical ordering of the faces.
 *        See element.h and element.cpp for a description of this ordering.
 *        Important: this implementation only works for simplicial meshes (Line, Triangle, Tet)
 *                   and will not work for polygonal or polyhedral meshes.
 */
template<typename type>
class Adjacency : public array2d<int> {

public:
  /**
   * \brief Initializes and builds the neighbor relations.
   *        Calls the 'build' function (private) below.
   */
  Adjacency( const Mesh<type>& mesh );

  /**
   * \brief Returns the number of elements neighboring any given element.
   *        This is equal to the number of faces of the element type.
   */
  int nb_neighbors() const {
    return type::nb_faces;
  }

  /**
   * \brief Returns the neighbor index of element k1 in element k0.
   *        For example, if k0 and k1 are adjacent to each other,
   *        and k1 is in k0's j'th neighbor, this function returns j.
   *        In other words, adj(k0,j) = k1.
   *
   * \return the index of k1 as the neighbor of k0, but -1 if k0 and k1 are not neighbors.
   */
  int indexof( int k0 , int k1 ) const;

private:

  /**
   * \brief Builds the neighbor relations.
   *        This class inherits from array2d<int> so the data in array2d<int>
   *        is filled with the adjacency data.
   */
  void build();

  const Mesh<type>& mesh_; // keep a const reference to the mesh although it's only needed during construction
};

} // flux

#endif
