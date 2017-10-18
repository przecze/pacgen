#include "gen_usdf.h"
#include "mesh.h"
#include "gen_kdtree.h"
#include <map>

namespace PG
{

 /**
 * Check if a number is positive or negative.
 * @param num The number.
 * @return 1 if the number if greater than zero, -1 otherwise.
 */
 static double GetSign(double num)
 {
  return (num < 0.0) ? -1.0 : 1.0;
 }

 /**
 * Finds the triangle's closest point from another point.
 * @param a First vertex.
 * @param b Second vertex.
 * @param c Third vertex.
 * @param p Point.
 * @param type 0->a. 1->b, 2->c, 3->ab, 4->bc, 5->ca, 6->abc
 * @return The triangle's closest point.
 */
 static vector3d ClosestPointToTriangle(const vector3d& a, const vector3d& b,
  const vector3d& c, const vector3d& p, int& type)
 {
  // Check if P in vertex region outside A
  vector3d ab = b - a;
  vector3d ac = c - a;
  vector3d ap = p - a;

  double d1 = vector3d::DotProduct(ab, ap);
  double d2 = vector3d::DotProduct(ac, ap);
  if (d1 <= 0.0f && d2 <= 0.0f)
  {
   type = 0;
   return a; // barycentric coordinates (1,0,0)
  }
  // Check if P in vertex region outside B
  vector3d bp = p - b;
  double d3 = vector3d::DotProduct(ab, bp);
  double d4 = vector3d::DotProduct(ac, bp);
  if (d3 >= 0.0f && d4 <= d3)
  {
   type = 1;
   return b; // barycentric coordinates (0,1,0)
  }
  // Check if P in edge region of AB, if so return projection of P onto AB
  double vc = d1 * d4 - d3 * d2;
  if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f)
  {
   double v = d1 / (d1 - d3);
   type = 3;
   return a + ab * v; // barycentric coordinates (1-v,v,0)
  }
  // Check if P in vertex region outside C
  vector3d cp = p - c;
  double d5 = vector3d::DotProduct(ab, cp);
  double d6 = vector3d::DotProduct(ac, cp);
  if (d6 >= 0.0f && d5 <= d6)
  {
   type = 2;
   return c; // barycentric coordinates (0,0,1)
  }
  // Check if P in edge region of AC, if so return projection of P onto AC
  double vb = d5 * d2 - d1 * d6;
  if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f)
  {
   double w = d2 / (d2 - d6);
   type = 5;
   return a + ac * w; // barycentric coordinates (1-w,0,w)
  }
  // Check if P in edge region of BC, if so return projection of P onto BC
  double va = d3 * d6 - d5 * d4;
  if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f)
  {
   double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
   type = 4;
   return b + (c - b) * w; // barycentric coordinates (0,1-w,w)
  }
  // P inside face region. Compute Q through its barycentric coordinates (u,v,w)
  double denom = 1.0f / (va + vb + vc);
  double v = vb * denom;
  double w = vc * denom;
  type = 6;
  return a + ab * v + ac * w; //=u*a+v*b+w*c,u=va*denom = 1.0f-v-w
 }

 /**
  * Point p intersects Segment v1-v2
  */
 static double GetXCoordSegment(const vector3d& v1, const vector3d& v2,
                                const vector3d& p)
 {
  double t;
  if (fabs(v2.y - v1.y) < 1.e-8)  t = (p.z - v1.z) / (v2.z - v1.z);
  else t = (p.y - v1.y) / (v2.y - v1.y);
  return v1.x + t * (v2.x - v1.x);
 }

 /**
  * Point p intersects Triangle plane v1-n
  */
 static double GetXCoordTriangle(const vector3d& v1, const vector3d& n,
                                 const vector3d& p)
 {
  return v1.x + (n.y *(v1.y - p.y) + n.z * (v1.z - p.z) ) / n.x;
 }

 static double LinearInterpolation(const double* const p,
  const double* const pos1, const double& dist1,
  const double* const pos2, const double& dist2, const int& dim)
 {
  const double t((p[dim] - pos1[dim]) / (pos2[dim] - pos1[dim]));
  return dist1 * (1.0 - t) + dist2 * t;
 }

 static double
  BilinearInterpolation(const double* const p,
  const double* const pos1, const double& dist1,
  const double* const pos2, const double& dist2,
  const double* const pos3, const double& dist3,
  const double* const pos4, const double& dist4, const int& dim1, const int& dim2)
 {
   double a_pos[3];
   for (int i = 0; i < 3; ++i) a_pos[i] = p[i];
   a_pos[dim2] = pos1[dim2];
   double a_dist(LinearInterpolation(a_pos, pos1, dist1, pos2, dist2, dim1));
   double b_pos[3];
   for (int i = 0; i < 3; ++i) b_pos[i] = p[i];
   b_pos[dim2] = pos3[dim2];
   double b_dist(LinearInterpolation(b_pos, pos3, dist3, pos4, dist4, dim1));
   return LinearInterpolation(p, a_pos, a_dist, b_pos, b_dist, dim2);
 }

 /**
 *        *3        *7                    011          111
 *
 *    *1        *5                     001          101
 *
 *                       z
 *        *2        *6   |  y             010          110
 *                       |/
 *    *0        *4       ----> x       000          100
 */
 static double Interpolation(const double* p, const double& x0, const double& y0,
  const double& z0, const double& x1, const double& y1, const double& z1,
  const double& v0, const double& v1, const double& v2, const double& v3,
  const double& v4, const double& v5, const double& v6, const double& v7)
 {
  const double xd = (p[0] - x0) / (x1 - x0);
  const double yd = (p[1] - y0) / (y1 - y0);
  const double zd = (p[2] - z0) / (z1 - z0);
  // along x
  const double c00 = v0 * (1 - xd) + v4 * xd;
  const double c01 = v1 * (1 - xd) + v5 * xd;
  const double c10 = v2 * (1 - xd) + v6 * xd;
  const double c11 = v3 * (1 - xd) + v7 * xd;
  // along y
  const double c0 = c00 * (1 - yd) + c10 * yd;
  const double c1 = c01 * (1 - yd) + c11 * yd;
  // along z
  return c0 * (1 - zd) + c1 * zd;
 }

/*************************************************************************/

USDFP::USDFP(const USDFP& s) : pos(s.pos), dist(s.dist) {}
USDFP::USDFP(double x, double y, double z, double dist) : pos(x, y, z), dist(dist) {}
bool USDFP::operator<(const USDFP& other) const
{
 return pos < other.pos;
}

/*************************************************************************/

typedef std::pair<vector3d, vector3d> KeyEdge;
struct KeyEdgeCmp
{
 bool operator()(const KeyEdge &l, const KeyEdge &r) const {
  return (l.first < r.first) || ( l.first == r.first && l.second < r.second);
 }
};

typedef std::map<KeyEdge, std::vector<double>, KeyEdgeCmp> EdgeMap;
typedef std::map<vector3d, std::vector<double>> VertexMap;

static inline double round(double number)
{
 return floor(number + 0.5);
}

/*************************************************************************/

USDF::USDF(void* mesh_container, const double& rmax, const int* const numpoints)
{
 for(int i = 0; i < 3; ++i) num_points[i] = numpoints[i];

 container = mesh_container;
 num_xz = num_points[0] * num_points[2];
 total_points = num_points[0] * num_points[1] * num_points[2];
 points.reserve(total_points);

 // Consider the surface displacement to create the grid of points
 Mesh3DMin(container, &points_min.x);
 Mesh3DMax(container, &points_max.x);
 points_min -= vector3d(rmax, rmax, rmax);
 points_max += vector3d(rmax, rmax, rmax);

 steps.x = (points_max.x - points_min.x) / static_cast<double>(num_points[0] - 1);
 steps.y = (points_max.y - points_min.y) / static_cast<double>(num_points[1] - 1);
 steps.z = (points_max.z - points_min.z) / static_cast<double>(num_points[2] - 1);

 inv_steps.x = 1.0 / steps.x;
 inv_steps.y = 1.0 / steps.y;
 inv_steps.z = 1.0 / steps.z;

 double diagonal = sqrt(steps.x * steps.x + steps.y * steps.y + steps.z * steps.z);
 neg_epsilon = (rmax + diagonal)*1.5;
 pos_epsilon = diagonal*1.5;

 initial_distances();
 compute_signed_distances();
}

USDF::~USDF()
{
 for (USDFP* p : points)
  delete p;
}

double USDF::GetDistance(const double* const cp)
{
 vector3d remain;
 for (int i = 0; i < 3; ++i) remain[i] = fmod(cp[i] - points_min[i], steps[i]);
 bool rx = remain.x < 1.e-8 || fabs(remain.x - steps.x) < 1.e-8;
 bool ry = remain.y < 1.e-8 || fabs(remain.y - steps.y) < 1.e-8;
 bool rz = remain.z < 1.e-8 || fabs(remain.z - steps.z) < 1.e-8;
 if ( rx && ry && rz ) // a point
 {
  return points[static_cast<int>(round( (cp[2] - points_min[2]) * inv_steps[2] ) +
                                 round( (cp[0] - points_min[0]) * inv_steps[0] ) * num_points[2] +
                                 round( (cp[1] - points_min[1]) * inv_steps[1] ) * num_xz)]->dist;
 }
 else if ( rx && ry ) // an edge - linear interpolation - on z axis
 {
  int p1 = static_cast<int>(floor( (cp[2] - points_min[2]) * inv_steps[2] ) +
                            round( (cp[0] - points_min[0]) * inv_steps[0] ) * num_points[2] +
                            round( (cp[1] - points_min[1]) * inv_steps[1] ) * num_xz);
  int p2 = p1 + 1;
  return LinearInterpolation(cp,
   &points[p1]->pos[0], points[p1]->dist,
   &points[p2]->pos[0], points[p2]->dist, 2);
 }
 else if ( ry && rz ) // an edge - linear interpolation - on x axis
 {
  int p1 = static_cast<int>(round( (cp[2] - points_min[2]) * inv_steps[2] ) +
                            floor( (cp[0] - points_min[0]) * inv_steps[0] ) * num_points[2] +
                            round( (cp[1] - points_min[1]) * inv_steps[1] ) * num_xz);
  int p2 = p1 + num_points[2];
  return LinearInterpolation(cp,
   &points[p1]->pos[0], points[p1]->dist,
   &points[p2]->pos[0], points[p2]->dist, 0);
 }
 else if ( rx && rz ) // an edge - linear interpolation - on y axis
 {
  int p1 = static_cast<int>(round( (cp[2] - points_min[2]) * inv_steps[2] ) +
                            round( (cp[0] - points_min[0]) * inv_steps[0] ) * num_points[2] +
                            floor( (cp[1] - points_min[1]) * inv_steps[1] ) * num_xz);
  int p2 = p1 + static_cast<int>(num_xz);
  return LinearInterpolation(cp,
   &points[p1]->pos[0], points[p1]->dist,
   &points[p2]->pos[0], points[p2]->dist, 1);
 }
 else if ( rx ) // a face - bilinear interpolation
 {
  int p1 = static_cast<int>(floor( (cp[2] - points_min[2]) * inv_steps[2] ) +
                            round( (cp[0] - points_min[0]) * inv_steps[0] ) * num_points[2] +
                            floor( (cp[1] - points_min[1]) * inv_steps[1] ) * num_xz);
  int p2 = p1 + 1;
  int p3 = p1 + static_cast<int>(num_xz);
  int p4 = p3 + 1;
  return BilinearInterpolation(cp,
   &points[p1]->pos[0], points[p1]->dist,
   &points[p2]->pos[0], points[p2]->dist,
   &points[p3]->pos[0], points[p3]->dist,
   &points[p4]->pos[0], points[p4]->dist, 2, 1);
 }
 else if ( ry ) // a face - bilinear interpolation
 {
  int p1 = static_cast<int>(floor( (cp[2] - points_min[2]) * inv_steps[2] ) +
                            floor( (cp[0] - points_min[0]) * inv_steps[0] ) * num_points[2] +
                            round( (cp[1] - points_min[1]) * inv_steps[1] ) * num_xz);
  int p2 = p1 + num_points[2];
  int p3 = p1 + 1;
  int p4 = p3 + num_points[2];
  return BilinearInterpolation(cp,
   &points[p1]->pos[0], points[p1]->dist,
   &points[p2]->pos[0], points[p2]->dist,
   &points[p3]->pos[0], points[p3]->dist,
   &points[p4]->pos[0], points[p4]->dist, 0, 2);
 }
 else if ( rz ) // a face - bilinear interpolation
 {
  int p1 = static_cast<int>(round( (cp[2] - points_min[2]) * inv_steps[2] ) +
                            floor( (cp[0] - points_min[0]) * inv_steps[0] ) * num_points[2] +
                            floor( (cp[1] - points_min[1]) * inv_steps[1] ) * num_xz);
  int p2 = p1 + num_points[2];
  int p3 = p1 + static_cast<int>(num_xz);
  int p4 = p3 + num_points[2];
  return BilinearInterpolation(cp,
   &points[p1]->pos[0], points[p1]->dist,
   &points[p2]->pos[0], points[p2]->dist,
   &points[p3]->pos[0], points[p3]->dist,
   &points[p4]->pos[0], points[p4]->dist, 0, 1);
 }
 else // trilinear interpolation
 {
  int p1 = static_cast<int>(floor( (cp[2] - points_min[2]) * inv_steps[2] ) +
                            floor( (cp[0] - points_min[0]) * inv_steps[0] ) * num_points[2] +
                            floor( (cp[1] - points_min[1]) * inv_steps[1] ) * num_xz);
  int p2 = p1 + static_cast<int>(num_xz);
  int p3 = p1 + 1;
  int p4 = p3 + static_cast<int>(num_xz);
  int p5 = p1 + num_points[2];
  int p6 = p5 + static_cast<int>(num_xz);
  int p7 = p5 + 1;
  int p8 = p7 + static_cast<int>(num_xz);
  return Interpolation(cp,
   points[p1]->pos.x, points[p1]->pos.y, points[p1]->pos.z,
   points[p8]->pos.x, points[p8]->pos.y, points[p8]->pos.z,
   points[p1]->dist, points[p3]->dist,
   points[p2]->dist, points[p4]->dist,
   points[p5]->dist, points[p7]->dist,
   points[p6]->dist, points[p8]->dist);
 }
}

void USDF::GetMin(double* m)
{
 Mesh3DMin(container, m);
}

void USDF::GetMax(double* m)
{
 Mesh3DMax(container, m);
}

void USDF::GetVolume(double& vol)
{
 int numIndices = Mesh3DNumIndices(container);
 double* meshVertP = Mesh3DVertexPointer(container);
 unsigned int* meshIdxP = Mesh3DIndexPointer(container);
 vol = 0.0;
 for (int i = 0; i < numIndices; i += 3)
 {
  const int ia = meshIdxP[i];
  const int ib = meshIdxP[i+1];
  const int ic = meshIdxP[i+2];
  const vector3d p1(meshVertP[3*ia], meshVertP[3*ia+1], meshVertP[3*ia+2]);
  const vector3d p2(meshVertP[3*ib], meshVertP[3*ib+1], meshVertP[3*ib+2]);
  const vector3d p3(meshVertP[3*ic], meshVertP[3*ic+1], meshVertP[3*ic+2]);
  double v321 = p3.x * p2.y * p1.z;
  double v231 = p2.x * p3.y * p1.z;
  double v312 = p3.x * p1.y * p2.z;
  double v132 = p1.x * p3.y * p2.z;
  double v213 = p2.x * p1.y * p3.z;
  double v123 = p1.x * p2.y * p3.z;
  vol += (1.0f / 6.0f) * (-v321 + v231 + v312 - v132 - v213 + v123);
 }
}

/**
 * Creates the grid of points and assign them an initial distance
 */
void USDF::initial_distances()
{
 struct Triangle {
  vector3d v[3];
  vector3d n;
 };

 std::vector<Triangle> m_triangles;
 std::vector<KDBox<3>*> boxes;

 /******************KD TREE************************/
 int numIndices = Mesh3DNumIndices(container);
 double* meshVertP = Mesh3DVertexPointer(container);
 double* meshTNormP = Mesh3DTnormalPointer(container);
 unsigned int* meshIdxP = Mesh3DIndexPointer(container);

 m_triangles.reserve(numIndices / 3);
 boxes.reserve(numIndices / 3);
 for (int i = 0; i < numIndices; i += 3)
 {
  //triangle vertex indices
  const int ia = meshIdxP[i];
  const int ib = meshIdxP[i + 1];
  const int ic = meshIdxP[i + 2];

  Triangle trian;
  trian.v[0] = vector3d(meshVertP[3*ia], meshVertP[3*ia+1], meshVertP[3*ia+2]);
  trian.v[1] = vector3d(meshVertP[3*ib], meshVertP[3*ib+1], meshVertP[3*ib+2]);
  trian.v[2] = vector3d(meshVertP[3*ic], meshVertP[3*ic+1], meshVertP[3*ic+2]);
  trian.n = vector3d(meshTNormP[(i/3)*3], meshTNormP[(i/3)*3+1], meshTNormP[(i/3)*3+2]);
  m_triangles.push_back(trian);

  KDBox<3>* box = new KDBox<3>();
  box->bmin[0] = box->bmax[0] = trian.v[0].x;
  box->bmin[1] = box->bmax[1] = trian.v[0].y;
  box->bmin[2] = box->bmax[2] = trian.v[0].z;
  for (int w = 1; w < 3; ++w)
  {
   for (int j = 0; j < 3; ++j) box->bmin[j] = PG_MIN(box->bmin[j], trian.v[w][j]);
   for (int j = 0; j < 3; ++j) box->bmax[j] = PG_MAX(box->bmax[j], trian.v[w][j]);
  }
  box->id = i / 3;
  boxes.push_back(box);
 }
 KDTNode<3>* kdtree = KDTNode<3>::Build(boxes);
 /******************KD TREE************************/

 // Distance points creation
 double y, z;
 y = points_min.y;
 for(int j = 0; j < num_points[1]; ++j)
 {
  double x = points_min.x;
  for(int i = 0; i < num_points[0]; ++i)
  {
   z = points_min.z;
   for(int k = 0; k < num_points[2]; ++k)
   {
    /***********************************************************/
    EdgeMap edges;
    VertexMap vertices;
    //For each point we create a ray
    vector3d ray_ini(x, y, z);
    vector3d ray_end(points_max.x, y, z);
    vector3d ray_dir(1,0,0);
    int count = 0;

    std::vector<int> ids;
    KDTNode<3>::Hit(kdtree, &ray_ini[0], &ray_end[0], ids);
    for (int it : ids)
    {
     Triangle f = m_triangles[it];

     vector3d newVa = f.v[0];
     vector3d newVb = f.v[1];
     vector3d newVc = f.v[2];

     //vector3d AA(0.0, 0.0, 0.0);
     vector3d BB(newVb.y - newVa.y, newVb.z - newVa.z, 0.0);
     vector3d CC(newVc.y - newVa.y, newVc.z - newVa.z, 0.0);
     vector3d PP(ray_ini.y - newVa.y, ray_ini.z - newVa.z, 0.0);
     double dd = BB.x * CC.y - CC.x * BB.y;
     double inv_dd = 1.0 / dd;
     double u = (PP.x * (BB.y - CC.y) + PP.y * (CC.x - BB.x) + dd) * inv_dd;
     double v = (PP.x * CC.y - PP.y * CC.x) * inv_dd;
     double w = (PP.y * BB.x - PP.x * BB.y) * inv_dd;

     if (u - 1.0 <= 1.e-8 && u >= -1.e-8 &&
         v - 1.0 <= 1.e-8 && v >= -1.e-8 &&
         w - 1.0 <= 1.e-8 && w >= -1.e-8)
     {
      bool is_edge = false;
      bool is_vertex = false;
      KeyEdge edge_key;
      vector3d vertex_key;

      if (fabs(u) <= 1.e-8 && fabs(v) <= 1.e-8)
      {
       if (newVc.x >= ray_ini.x)
       {
        is_vertex = true;
        vertex_key = newVc;
       }
      }
      else if (fabs(u) <= 1.e-8 && fabs(w) <= 1.e-8)
      {
       if (newVb.x >= ray_ini.x)
       {
        is_vertex = true;
        vertex_key = newVb;
       }
      }
      else if (fabs(w) <= 1.e-8 && fabs(v) <= 1.e-8)
      {
       if (newVa.x >= ray_ini.x)
       {
        is_vertex = true;
        vertex_key = newVa;
       }
      }
      else if ( fabs(u) <= 1.e-8 )
      {
       if (GetXCoordSegment(newVb, newVc, ray_ini) >= ray_ini.x)
       {
        is_edge = true;
        if (newVb < newVc) edge_key = std::make_pair(newVb, newVc);
        else               edge_key = std::make_pair(newVc, newVb);
       }
      }
      else if ( fabs(v) <= 1.e-8 )
      {
       if (GetXCoordSegment(newVa, newVc, ray_ini) >= ray_ini.x)
       {
        is_edge = true;
        if (newVa < newVc) edge_key = std::make_pair(newVa, newVc);
        else               edge_key = std::make_pair(newVc, newVa);
       }
      }
      else if ( fabs(w) <= 1.e-8 )
      {
       if (GetXCoordSegment(newVa, newVb, ray_ini) >= ray_ini.x)
       {
        is_edge = true;
        if (newVa < newVb) edge_key = std::make_pair(newVa, newVb);
        else               edge_key = std::make_pair(newVb, newVa);
       }
      }
      else
      {
       if (GetXCoordTriangle(newVa, f.n, ray_ini) >= ray_ini.x)
        count++;
      }

      if ( is_edge )
      {
       std::vector<double> values;
       if ( edges.count(edge_key) ) values = edges[edge_key];
       double dotRayTriNormal = vector3d::DotProduct( ray_dir, f.n );
       values.push_back( dotRayTriNormal );
       edges[edge_key] = values;
      }

      if ( is_vertex )
      {
       std::vector<double> values;
       if ( vertices.count(vertex_key) ) values = vertices[vertex_key];
       double dotRayTriNormal = vector3d::DotProduct( ray_dir, f.n );
       values.push_back( dotRayTriNormal );
       vertices[vertex_key] = values;
      }
     }
    }
    // Checking special cases when intersecting multiple vertices or edges
    for (auto& iterator : edges)
    {
     std::vector<double> tempt = iterator.second;
     int sis = (int) tempt.size();
     if ( sis == 1 ) count++;
     else
     {
      double sign = GetSign(tempt[0]);
      bool diff = false;
      for (int iw = 1; iw < sis && !diff; iw++)
       if ( GetSign(tempt[iw]) != sign )
        diff = true;
      if (!diff) count++;
     }
    }
    for (auto& iterator : vertices)
    {
     std::vector<double> tempt = iterator.second;
     int sis = (int) tempt.size();
     if ( sis == 1 ) count++;
     else
     {
      double sign = GetSign(tempt[0]);
      bool diff = false;
      for (int iw = 1; iw < sis && !diff; iw++)
       if ( GetSign(tempt[iw]) != sign )
        diff = true;
      if (!diff) count++;
     }
    }
    double ini_distance = ( count % 2 == 0 ) ? MAX_DISTANCE : MIN_DISTANCE;
    /***********************************************************/
    points.push_back(new USDFP(x, y, z, ini_distance));
    z += steps.z;
   }
   x += steps.x;
  }
  y += steps.y;
 }
 delete kdtree;
 for (int i = 0; i < numIndices/3; ++i) delete boxes[i];
}

std::vector<USDFP*> USDF::getIntersectingPoints(
 const vector3d& bmin, const vector3d& bmax)
{
 std::vector<USDFP*> found;

 int dminx = static_cast<int>(floor((bmin.x - points_min.x) * inv_steps[0]));
 int dminy = static_cast<int>(floor((bmin.y - points_min.y) * inv_steps[1]));
 int dminz = static_cast<int>(floor((bmin.z - points_min.z) * inv_steps[2]));

 int dmaxx = static_cast<int>(floor((bmax.x - points_min.x) * inv_steps[0]));
 int dmaxy = static_cast<int>(floor((bmax.y - points_min.y) * inv_steps[1]));
 int dmaxz = static_cast<int>(floor((bmax.z - points_min.z) * inv_steps[2]));

 for (int i = dminx; i <= dmaxx; ++i)
 {
  for (int j = dminy; j <= dmaxy; ++j)
  {
   for (int k = dminz; k <= dmaxz; ++k)
   {
    int c = k + i * num_points[2] + j * static_cast<int>(num_xz);
    if (c < total_points && c >= 0)
    {
     USDFP* point = points[c];
     if ((point->pos.x >= bmin.x && point->pos.x <= bmax.x) &&
         (point->pos.y >= bmin.y && point->pos.y <= bmax.y) &&
         (point->pos.z >= bmin.z && point->pos.z <= bmax.z))
      found.push_back(point);
    }
   }
  }
 }
 return found;
}

void USDF::compute_signed_distances()
{
 std::vector<vector3d> prism_min;
 std::vector<vector3d> prism_max;

 { // prisms creation
  double* meshVertP = Mesh3DVertexPointer(container);
  double* meshTNormP = Mesh3DTnormalPointer(container);
  unsigned int* meshIdxP = Mesh3DIndexPointer(container);
  size_t n_indices = Mesh3DNumIndices(container);

  prism_min.reserve(n_indices / 3);
  prism_max.reserve(n_indices / 3);
  for (size_t i = 0; i < n_indices; i += 3)
  {
   const int ia = meshIdxP[i];
   const int ib = meshIdxP[i + 1];
   const int ic = meshIdxP[i + 2];

   const vector3d va(meshVertP[3*ia], meshVertP[3*ia+1], meshVertP[3*ia+2]);
   const vector3d vb(meshVertP[3*ib], meshVertP[3*ib+1], meshVertP[3*ib+2]);
   const vector3d vc(meshVertP[3*ic], meshVertP[3*ic+1], meshVertP[3*ic+2]);

   const vector3d normal(meshTNormP[(i/3)*3], meshTNormP[(i/3)*3+1], meshTNormP[(i/3)*3+2]);

   vector3d prism_points[6];
   prism_points[0] = va + normal * pos_epsilon;
   prism_points[1] = vb + normal * pos_epsilon;
   prism_points[2] = vc + normal * pos_epsilon;
   prism_points[3] = va - normal * neg_epsilon;
   prism_points[4] = vb - normal * neg_epsilon;
   prism_points[5] = vc - normal * neg_epsilon;

   double xmin, ymin, zmin, xmax, ymax, zmax;
   xmin = xmax = prism_points[0].x;
   ymin = ymax = prism_points[0].y;
   zmin = zmax = prism_points[0].z;
   for (int i = 1; i < 6; i++)
   {
    if (prism_points[i].x < xmin) xmin = prism_points[i].x;
    if (prism_points[i].x > xmax) xmax = prism_points[i].x;
    if (prism_points[i].y < ymin) ymin = prism_points[i].y;
    if (prism_points[i].y > ymax) ymax = prism_points[i].y;
    if (prism_points[i].z < zmin) zmin = prism_points[i].z;
    if (prism_points[i].z > zmax) zmax = prism_points[i].z;
   }
   prism_min.push_back({xmin - steps.x * 0.5,
                        ymin - steps.y * 0.5, zmin - steps.z * 0.5});
   prism_max.push_back({xmax + steps.x * 0.5,
                        ymax + steps.y * 0.5, zmax + steps.z * 0.5});
  }
 }

 double* meshVertP = Mesh3DVertexPointer(container);
 double* meshTNormP = Mesh3DTnormalPointer(container);
 unsigned int* meshIdxP = Mesh3DIndexPointer(container);

 std::vector<int> triangles(points.size(), -1);
 /**
  * Calculating the signed distance from each point to the mesh. We read each
  * triangle and use the bounding box of its prism to get all the distance points inside.
  * For each point we calculate the distance to the current triangle.
  * If the absolute value of the distance is lower than the current point distance, we save it.
  */
 const int n = (int) prism_min.size();//one prism for each triangle
 for (int i = 0 ; i < n; i++)
 {
  vector3d normal(meshTNormP[i*3], meshTNormP[i*3+1], meshTNormP[i*3+2]);

  std::vector<USDFP*> closer = getIntersectingPoints(prism_min[i], prism_max[i]);
  int type;
  for (USDFP* point : closer)
  {
   double xx = (point->pos.x - points_min.x) / steps[0];
   double yy = (point->pos.y - points_min.y) / steps[1];
   double zz = (point->pos.z - points_min.z) / steps[2];
   int tindex = static_cast<int>(round(num_points[2] * xx + num_points[0] * num_points[2] * yy + zz));

   vector3d p = point->pos;
   const int ia(meshIdxP[3*i]);
   const int ib(meshIdxP[3*i + 1]);
   const int ic(meshIdxP[3*i + 2]);
   vector3d closest(ClosestPointToTriangle({meshVertP[3*ia], meshVertP[3*ia+1], meshVertP[3*ia+2]},
                                           {meshVertP[3*ib], meshVertP[3*ib+1], meshVertP[3*ib+2]},
                                           {meshVertP[3*ic], meshVertP[3*ic+1], meshVertP[3*ic+2]},
                                           p, type));
   const double newDistance = (p - closest).length() * GetSign(vector3d::DotProduct(normal, p - closest));
   if (fabs(newDistance) < fabs(point->dist))
   {
    if ((point->dist == MAX_DISTANCE || point->dist == MIN_DISTANCE) ||
         GetSign(newDistance) == GetSign(point->dist))
    {
     point->dist = newDistance;
     triangles[tindex] = i;
    }
    else
    {
     vector3d v(meshTNormP[triangles[tindex]*3],
                meshTNormP[triangles[tindex]*3+1], meshTNormP[triangles[tindex]*3+2]);

     double angle = PG_RAD2DEG(acos(vector3d::DotProduct(normal,v) / (normal.length() * v.length())));
     if (point->dist > 0.0 && newDistance < 0.0) //positive to negative case
     {
      if (angle >= -180.0 && angle <= 0.0) // valid between [0, -180]
      {
       point->dist = newDistance;
       triangles[tindex] = i;
      }
     }
     else if (point->dist < 0.0 && newDistance > 0.0) //negative to positive case
     {
      if (angle >= 0.0 && angle <= 180.0) // valid between [0, 180]
      {
       point->dist = newDistance;
       triangles[tindex] = i;
      }
     }
    }
   }
  }
 }

 for (USDFP* p : points)
 {
  if      (p->dist == MAX_DISTANCE) p->dist = MIN_DISTANCE;
  else if (p->dist == MIN_DISTANCE) p->dist = MAX_DISTANCE;
  else                              p->dist = -p->dist;
 }
}

}