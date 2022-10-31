/**********************************************************************
 * flux: framework for learning unstructured meshing
 * Copyright (c) 2022 Philip Caplan. All rights reserved.
 * Licensed under the MIT License (https://mit-license.org/)
 **********************************************************************/
#ifndef FLUX_VERTICES_H_
#define FLUX_VERTICES_H_

#include "array2d.h"

namespace flux  {

/**
 * \brief Represents mesh vertices with optional metadata
 *        that may include vertex normals, vertex parameters (uv)
 *        and vertex-facet incidence relations (possibly for polygonal
 *        or polyhedral meshes)
 *        Note that normals are always in 3d and uv are always in 2d
 *        since these are mostly used for visualization.
 */
class Vertices : public array2d<double> {

public:
  /**
   * \brief Constructor from a specified dimension.
   *        Despite the vertices being an any ambient dimension,
   *        the normals are always 3d and uv (texture coordinates) are always 2d
   *
   * \param[in] dim - ambient dimension of the vertices
   */
  Vertices( int dim = 0 ) :
    array2d(dim),
    uv_(2),
    normals_(3),
    incidence_(-1)
  {}

  /**
   * \brief sets the ambient dimension of the vertices
   *        by setting the stride in the base array2d<double> class
   *
   * \param[in] dim - the ambient dimension.
   */
  void set_dim( int dim ) { set_stride(dim); }

  /**
   * \brief returns the dimension of the vertices.
   *
   * \return dimension
   */
  int dim() const { return stride(); }

  /**
   * \brief Copies the input vertices into this.
   *
   * \param[in] vertices - vertices to be copied
   */
  void copy( const Vertices& vertices );

  /**
   * \brief Retrieves a const reference to the vertex normals.
   */
  const array2d<double>& normals() const { return normals_; }

  /**
   * \brief Retrieves a non-const reference to the vertex normals.
   */
  array2d<double>& normals() { return normals_; }

  /**
   * \brief Retrieves a const reference to the vertex parameters.
   */
  const array2d<double>& uv() const { return uv_; }

  /**
   * \brief Retrieves a non-const reference to the vertex parameters.
   */
  array2d<double>& uv() { return uv_; }


  /**
   * \brief Clears any data stored at vertices (coordinates, uv, normals, etc.)
   */
  void clear();

private:

  // this class inherits from array2d<double> which is where the actual coordinate data is stored
  array2d<double> uv_;        // (u,v) or texture coordinate data useful for parametrizations
  array2d<double> normals_;   // 3d normals useful for visualization
  array2d<int>    incidence_; // vertex-facet incidence matrix useful for polygonal or polyhedral meshes
};

} // flux

#endif
