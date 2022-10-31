/**********************************************************************
 * flux: framework for learning unstructured meshing
 * Copyright (c) 2022 Philip Caplan. All rights reserved.
 * Licensed under the MIT License (https://mit-license.org/)
 **********************************************************************/
#ifndef FLUX_ELEMENT_H_
#define FLUX_ELEMENT_H_

#include <vector>

namespace flux {

/**
* \brief Vertex (Node)
*
*  0  (no edges, faces or triangles)
*/
struct Vertex {
  static const int dimension = 0;
  static const int nb_vertices = 1;
  static const int nb_edges = 0;
  static const int nb_faces = 0;
  static const int nb_triangles = 0;
  static int* edges;
  static int* triangles;
  static int* faces;
};

/**
* \brief Line element
*
*   face/edge 0
* 0 ---------- 1
*/
struct Line {
  static const int dimension = 1;
  static const int nb_vertices = 2;
  static const int nb_edges = 1;
  static const int nb_faces = 2;
  static const int nb_triangles = 0;
  static int edges[2];
  static int* triangles;
  static int faces[2];
  typedef Vertex face_type;
};

/**
* \brief Triangle element:
*
*        2
*        |
*        | \
* face/  |   \     face/edge 1
* edge   |     \
*  2     |  t0   \
*        0 ------- 1
*         face/edge 0
*/
struct Triangle {
  static const int dimension = 2;
  static const int nb_vertices = 3;
  static const int nb_edges = 3;
  static const int nb_faces = 3;
  static const int nb_triangles = 1;
  static int edges[6];
  static int triangles[3];
  static int faces[6];
  typedef Line face_type;
};

/**
* \brief Quadrilateral element:
*
*         face/edge 2
*        3 --------- 2
*        |   t1    / |
*        |       /   |
* face/  |     /     |   face/edge 1
* edge   |   /  t0   |
*  3     | /         |
*        0 --------- 1
*         face/edge 0
*/
struct Quad {
  static const int dimension = 2;
  static const int nb_vertices = 4;
  static const int nb_edges = 4;
  static const int nb_faces = 4;
  static const int nb_triangles = 2;
  static int edges[8];
  static int triangles[6];
  static int faces[8];
  typedef Line face_type;
};

/**
* \brief Tetrahedron element: (please see element.cpp for the triangles & edges)
*
*                    y
*                  .
*                ,/
*               /
*            2
*          ,/|`\
*        ,/  |  `\
*      ,/    '.   `\
*    ,/       |     `\
*  ,/         |       `\
* 0-----------'.--------1 --> x
*  `\.         |      ,/
*     `\.      |    ,/
*        `\.   '. ,/
*           `\. |/
*              `3
*                 `\.
*                    ` z
*/
struct Tet {
  static const int dimension = 3;
  static const int nb_vertices = 4;
  static const int nb_edges = 6;
  static const int nb_faces = 4;
  static const int nb_triangles = 4;
  static int edges[12];
  static int triangles[12];
  static int faces[12];
  typedef Triangle face_type;
};

/**
* \brief Polygon and Polyhedral elements
*
* By not defining certain attributes of a Polygon and Polyhedron
* we will get a compile-time error if we try to use functions that
* are not designed to be used with them
*/
struct Polygon {
  static const int dimension = 2;
  static const int nb_vertices = -1; // the number of vertices per polygon is not a known constant
  typedef Line facet_type;
};

struct Polyhedron {
  static const int dimension = 3;
  static const int nb_vertices = -1; // the number of vertices per polyhedron is not a known constant
  typedef Polygon face_type;
};

/**
* \brief returns a random point inside an element (specified by the template parameter 'type')
*
* \param[out] alpha the vector of barycentric coordinates (size = type::dimension)
*/
template<typename type> void random_barycentric( std::vector<double>& alpha );

} // flux

#endif
