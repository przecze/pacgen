#ifndef _MESH_H_
#define _MESH_H_

namespace PG
{

 /**
  * Creates a triangle mesh.
  * @param v   Vertex array: x_1, y_1, z_1, x_2, y_2, z_2, ..., x_nv, y_nv, z_nv.
  * @param nv  Number of vertices: sizeof(v)/3.
  * @param id  Index array: t1_a, t1_b, t1_c, t2_a, t2_b, t2_c, ..., tnid_a, tnid_b, tnid_c
  * @param nid Number of indices: sizeof(id)
  * @param The mesh.
  */
 void* Mesh3DCreate(const double* v, const int& nv, const int* id, const int& nid);
 void Mesh3DDestroy(void* mesh);
 /**
  * Gets the minimum corner of the AABB.
  * @param mesh Input: The mesh.
  * @param pmin Output: Sets the minimum corner.
  */
 void Mesh3DMin(void* mesh, double* pmin);
 /**
  * Gets the maximum corner of the AABB.
  * @param mesh Input: The mesh.
  * @param pmin Output: Sets the maximum corner.
  */
 void Mesh3DMax(void* mesh, double* pmax);
 /**
  * @param mesh Input: The mesh.
  * @return The number of vertices.
  */
 int Mesh3DNumVertices(void* mesh);
 /**
  * @param mesh Input: The mesh.
  * @return The pointer to the vertex vector (x1, y1, z1, ..., xnv, ynv, znv).
  */
 double* Mesh3DVertexPointer(void* mesh);
 /**
  * @param mesh Input: The mesh.
  * @return The pointer to the triangle normals vector
  *         (nx_1, ny_1, nz_1, ..., nx_(nid/3), ny_(nid/3), nz_(nid/3)).
  */
 double* Mesh3DTnormalPointer(void* mesh);
 /**
  * @param mesh Input: The mesh.
  * @return The pointer to the vertex normals vector.
  *         (nx_1, ny_1, nz_1, ..., nx_nv, ny_nv, nz_nv).
  */
 double* Mesh3DVnormalPointer(void* mesh);
 /**
  * @param mesh Input: The mesh.
  * @return The number of indices.
  */
 int Mesh3DNumIndices(void* mesh);
 /**
  * @param mesh Input: The mesh.
  * @return The pointer to the index vector.
  */
 unsigned int* Mesh3DIndexPointer(void* mesh);

}

#endif /* _MESH_H_ */