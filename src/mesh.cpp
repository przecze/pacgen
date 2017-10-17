#include "mesh.h"
#include "gen_utils.h"
#include <algorithm>
#include <cfloat>
#include <map>

namespace PG
{
#define COMPUTE_VERTEX_NORMALS

 struct Mesh3D
 {
  Mesh3D();
  void computeNormals();

  std::vector<vector3d> m_vertices;
  std::vector<vector3d> m_tnormals; /*!< triangle normals. */
 #ifdef COMPUTE_VERTEX_NORMALS
  std::vector<vector3d> m_normals;  /*!< vertex normals.   */
 #endif
  vector3d            m_min; /*!< minimum corner of the mesh AABB. */
  vector3d            m_max; /*!< maximum corner of the mesh AABB. */
  std::vector<uint32> m_indices;
 };

 Mesh3D::Mesh3D() : m_min({  DBL_MAX,  DBL_MAX,  DBL_MAX }),
                    m_max({ -DBL_MAX, -DBL_MAX, -DBL_MAX }) {}

 void Mesh3D::computeNormals()
 {
  int n_vertices((int) m_vertices.size());
  for (int i = 0; i < n_vertices; ++i)
   for (int ww = 0; ww < 3; ++ww)
   {
    m_max[ww] = PG_MAX(m_vertices[i][ww], m_max[ww]);
    m_min[ww] = PG_MIN(m_vertices[i][ww], m_min[ww]);
   }

#ifdef COMPUTE_VERTEX_NORMALS
  m_normals.resize(n_vertices, vector3d(0.0, 0.0, 0.0));
  std::map<std::pair<uint32, uint32>, vector3d> edgemap;
#endif
  uint32 n_indices = (uint32) m_indices.size();
  m_tnormals.resize(n_indices / 3);

  for (uint32 i = 0; i < n_indices; i += 3)
  {
   const uint32 ic = m_indices[i + 2];
   const uint32 ib = m_indices[i + 1];
   const uint32 ia = m_indices[i];
   const vector3d va = m_vertices[ia];
   const vector3d vb = m_vertices[ib];
   const vector3d vc = m_vertices[ic];

   const vector3d tnormal = vector3d::Normalize(vector3d::CrossProduct(vb - va, vc - va));

#ifdef COMPUTE_VERTEX_NORMALS
   vector3d w, k;
   double angle;

   w = vector3d::Normalize(vb - va);
   k = vector3d::Normalize(vc - va);
   angle = acos( vector3d::DotProduct(w, k) );
   m_normals[ia] = vector3d::Normalize(m_normals[ia] + tnormal * angle);

   w = vector3d::Normalize(va - vb);
   k = vector3d::Normalize(vc - vb);
   angle = acos( vector3d::DotProduct(w, k) );
   m_normals[ib] = vector3d::Normalize(m_normals[ib] + tnormal * angle);

   w = vector3d::Normalize(va - vc);
   k = vector3d::Normalize(vb - vc);
   angle = acos( vector3d::DotProduct(w, k) );
   m_normals[ic] = vector3d::Normalize(m_normals[ic] + tnormal * angle);

   const uint32 v[] = {ia, ib, ic};
   for (int j = 0; j < 3; j++)
   {
    uint32 _a = v[j];
    uint32 _b = v[(j == 2) ? 0 : (j + 1)];
    if(_b < _a) std::swap(_a, _b);
    vector3d en;
    auto it = edgemap.find(std::pair<uint32, uint32>(_a, _b));
    if(it == edgemap.end())
     en = tnormal * M_PI;
    else
     en = it->second + tnormal * M_PI;
    edgemap[std::pair<uint32, uint32>(_a, _b)] = vector3d::Normalize(en);
   }
#endif

   m_tnormals[i / 3] = tnormal;
  }
 }

 /***********************************************************************/

 void* Mesh3DCreate(const double* v, const int& nv, const int* id, const int& nid)
 {
  Mesh3D* mesh = new Mesh3D();

  mesh->m_vertices.reserve(nv);
  for (int i = 0; i < nv; ++i)
   mesh->m_vertices.push_back(vector3d(v[i*3], v[i*3+1], v[i*3+2]));

  mesh->m_indices.reserve(nid);
  for (int i = 0; i < nid; ++i)
   mesh->m_indices.push_back(id[i]);

  mesh->computeNormals();
  return static_cast<void*>(mesh);
 }

 void Mesh3DDestroy(void* mesh)
 {
  Mesh3D* m = static_cast<Mesh3D*>(mesh);
  if (m) delete m;
 }

 void Mesh3DMin(void* mesh, double* pmin)
 {
  Mesh3D* m = static_cast<Mesh3D*>(mesh);
  if (m)
   for (int i = 0; i < 3; ++i)
    pmin[i] = m->m_min[i];
 }

 void Mesh3DMax(void* mesh, double* pmax)
 {
  Mesh3D* m = static_cast<Mesh3D*>(mesh);
  if (m)
   for (int i = 0; i < 3; ++i)
    pmax[i] = m->m_max[i];
 }

 int Mesh3DNumVertices(void* mesh)
 {
  Mesh3D* m = static_cast<Mesh3D*>(mesh);
  if (m) return (int) m->m_vertices.size();
  return 0;
 }

 double* Mesh3DVertexPointer(void* mesh)
 {
  Mesh3D* m = static_cast<Mesh3D*>(mesh);
  if (m) return &m->m_vertices[0].x;
  return nullptr;
 }

 double* Mesh3DTnormalPointer(void* mesh)
 {
  Mesh3D* m = static_cast<Mesh3D*>(mesh);
  if (m) return &m->m_tnormals[0].x;
  return nullptr;
 }

 double* Mesh3DVnormalPointer(void* mesh)
 {
#ifdef COMPUTE_VERTEX_NORMALS
  Mesh3D* m = static_cast<Mesh3D*>(mesh);
  if (m) return &m->m_normals[0].x;
#endif
  return nullptr;
 }

 int Mesh3DNumIndices(void* mesh)
 {
  Mesh3D* m = static_cast<Mesh3D*>(mesh);
  if (m) return (int) m->m_indices.size();
  return 0;
 }

 unsigned int* Mesh3DIndexPointer(void* mesh)
 {
  Mesh3D* m = static_cast<Mesh3D*>(mesh);
  if (m) return &m->m_indices[0];
  return nullptr;
 }

}