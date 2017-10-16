#include "gen_utils.h"
#ifdef _WIN32
#include <sys\timeb.h>
#else
#include <sys/time.h>
#endif
#include <time.h>

namespace PG
{

/*
* When two spheres intersect they form a circle whose center is stored at "c3"
* and radius is stored at "r3". The circle also has a support which is defined
* by a plane called as "radical plane" and its equation is stored in the variables
* "A, B, C and D".
*
* @param err Limit value allowed to consider contacts.
* @param c1 First sphere center.
* @param r1 First sphere radius.
* @param c2 Second sphere center.
* @param r2 Second sphere radius.
* @param points List of intersection points.
* @parma p number of points in the halo intersection.
*/
bool FindIntersectionPoints(const double& err, const double* c1,
                            const double& r1, const double* c2,
                            const double& r2, std::vector<vector3d>& points,
                            const unsigned int& p)
{
 vector3d c3, u, v, n3;
 double r3;
 bool special;
 if (!FindIntersectingCircle(
  err, c1, r1, c2, r2, &c3[0], r3, &n3[0], &u[0], &v[0], special))
  return false;
 if (special)
 {
  points.push_back(c3);
 }
 else
 {
  double segments = 180.0;
  double step = segments / (double)p;
  double tmp = 2.0 * double(M_PI) / segments;
  for (double ii = 0; ii <= segments; ii = ii + step)
  {
   double theta = ii * tmp;
   vector3d vecc = c3 + (u * cos(theta) + v * sin(theta)) * r3;
   points.push_back(vecc);
  }
 }
 return true;
}


/*
* Determines if two spheres have points of intersections.
* @param err Limit value allowed to consider contacts.
* @param c1 First sphere center.
* @param r1 First sphere radius.
* @param c2 Second sphere center.
* @param r2 Second sphere radius.
* @param c3 Intersection circle center.
* @param r3 Intersection circle radius.
* @param n3 Intersection circle normal.
* @param u Orthogonal vector on the circle's plane.
* @param v Orthogonal vector on the circle's plane.
* @param special If the intersection is only a point.
* @return True if an intersection exits, false otherwise.
*/
bool FindIntersectingCircle(const double& err, const double* p1,
                            const double& r1, const double* p2,
                            const double& r2, double* c3, double& r3,
                            double* n3, double* u, double* v,
                            bool& special)
{
 double er(err);
 if (fabs(p1[0]-p2[0])<er && fabs(p1[1]-p2[1])<er && fabs(p1[2]-p2[2])<er)
  return false;
 special = false;

 double rr1(r1);
 double rr2(r2);
 double dist2(0.0);
 for (int i = 0; i < 3; ++i) dist2 += (p1[i] - p2[i]) * (p1[i] - p2[i]);

 if (((rr1 + rr2) * (rr1 + rr2) - dist2) > 0.0)
 {
  const double a0 = (1.0 - ((rr1 * rr1 - rr2 * rr2) / dist2)) * 0.5;
  double p3[] = {0, 0, 0};
  for (int i = 0; i < 3; ++i) p3[i] = p1[i] * a0 + p2[i] * (1.0 - a0);

  double _dist_(0.0);
  for (int i = 0; i < 3; ++i) _dist_ += (p3[i] - p1[i]) * (p3[i] - p1[i]);
  r3 = sqrt((rr1 * rr1) - _dist_);

  double sqrt_dist = sqrt(dist2);
  double vv[] = {0, 0, 0};
  double nn[] = {0, 0, 0};
  for (int i = 0; i < 3; ++i) nn[i] = (p2[i] - p1[i]) / sqrt_dist;

  /**
   * Calculating the orthogonal vectors U and V
   */
  if (nn[0] <= nn[1] && nn[0] <= nn[2])
  {
   vv[0] = 0;
   vv[1] = nn[2];
   vv[2] = -nn[1];
   if (fabs(vv[1]) <= 1.e-10 && fabs(vv[2]) <= 1.e-10)
   {
    vv[0] = -nn[2];
    vv[1] = 0;
    vv[2] = nn[0];
   }
  }
  else if (nn[1] <= nn[0] && nn[1] <= nn[2])
  {
   vv[0] = -nn[2];
   vv[1] = 0;
   vv[2] = nn[0];
   if (fabs(vv[0]) <= 1.e-10 && fabs(vv[2]) <= 1.e-10)
   {
    vv[0] = nn[1];
    vv[1] = -nn[0];
    vv[2] = 0;
   }
  }
  else if (nn[2] <= nn[0] && nn[2] <= nn[1])
  {
   vv[0] = nn[1];
   vv[1] = -nn[0];
   vv[2] = 0;
   if (fabs(vv[0]) <= 1.e-10 && fabs(vv[1]) <= 1.e-10)
   {
    vv[0] = 0;
    vv[1] = nn[2];
    vv[2] = -nn[1];
   }
  }

  _dist_ = 0.0;
  for (int i = 0; i < 3; ++i) _dist_ += vv[i] * vv[i];
  sqrt_dist = sqrt(_dist_);
  for (int i = 0; i < 3; ++i) vv[i] = vv[i] / sqrt_dist;

  double uu[] = {0, 0, 0};
  uu[0] = vv[1] * nn[2] - nn[1] * vv[2];
  uu[1] = vv[2] * nn[0] - nn[2] * vv[0];
  uu[2] = vv[0] * nn[1] - nn[0] * vv[1];
  _dist_ = 0.0;
  for (int i = 0; i < 3; ++i) _dist_ += uu[i] * uu[i];
  sqrt_dist = sqrt(_dist_);
  for (int i = 0; i < 3; ++i) uu[i] = uu[i] / sqrt_dist;

  for (int i = 0; i < 3; ++i) n3[i] = nn[i];
  for (int i = 0; i < 3; ++i) u[i] = uu[i];
  for (int i = 0; i < 3; ++i) v[i] = vv[i];
  for (int i = 0; i < 3; ++i) c3[i] = p3[i];

  return true;
 }
 return false;
}

/*
 * Determines if three spheres have points of intersections.
 * Firts intersect two spheres and then intersect the intersection circle with
 * the third sphere.
 * @param err Limit value allowed to consider contacts.
 * @param c1 First sphere center.
 * @param r1 First sphere radius.
 * @param c2 Second sphere center.
 * @param r2 Second sphere radius.
 * @param c3 Third sphere center.
 * @param r3 Third sphere radius.
 * @param points Intersection points.
 */
bool FindIntersectionPoints(const double& err, const double* c1,
                            const double& r1, const double* c2,
                            const double& r2, const vector3d& c3,
                            const double& r3, std::vector<vector3d>& points)
{
 /**
 * Data for the intersection circle between the first two spheres.
 */
 vector3d c4, u, v, n4;
 double r4;
 bool special;

 /**
 * Certifying that the sphere 1 and 2 intersect forming a circle o single point.
 */
 if (!FindIntersectingCircle(
  err, c1, r1, c2, r2, &c4[0], r4, &n4[0], &u[0], &v[0], special))
  return false;

 if (special) //sphere 1 and 2 intersection generate a single point c4.
 {
  //If the distance from c4 to c3 is equal to the r3 radius.
  if (fabs(vector3d::DotProduct(c4 - c3, c4 - c3) - r3 * r3) <= err)
  {
   points.push_back(c4);
   return true;
  }
 }
 else //sphere 1 and 2 intersection generate the circle (c4, r4)
 {
  /**
  * Must perform the intersection of the circle's(c4, r4) plane with
  * the third sphere(c3,r3)
  */
  double fact1 = n4.x * c3.x + n4.y * c3.y + n4.z * c3.z + -vector3d::DotProduct(n4, c4);
  double fact2 = n4.x * n4.x + n4.y * n4.y + n4.z * n4.z;
  double d = fabs(fact1) / sqrt(fact2);

  //the plane of the circle (c4, r4) is tangent to the sphere 3
  if (fabs(r3 - d) <= err)
  {
   /**
   * Computing the touching point c on sphere 3.
   */
   vector3d c = c3 - n4 * (fact1 / fact2);
   //radius r4 is equal to the distance between of circle (c4, r4) and point c
   if (fabs((r4 * r4) - vector3d::DotProduct(c - c4, c - c4)) <= err)
   {
    points.push_back(c);
    return true;
   }
  }
  else if (r3 - d > 0) //the plane of the circle (c4, r4) intersects the sphere 3
   //must intersect the circle(c4, r4) vs circle (c,r) (both in the same plane)
  {
   /**
   * Computing the circle center c
   */

   vector3d c;
   if (vector3d::DotProduct(c4-c3, n4) >= 0)
    c = c3 - n4 * (fact1 / fact2);
   else
    c = c3 + n4 * (fact1 / fact2);

   /**
   * If circles are in the same position then when do not consider the intersection.
   */
   if (fabs(c.x - c4.x) < err &&
       fabs(c.y - c4.y) < err &&
       fabs(c.z - c4.z) < err)
    return false;

   /**
   * Computing the radius center r
   */
   double r = sqrt(r3 * r3 - d * d);

   /**
   * Distance between the circle(c4, r4) and the circle(c, r)
   */
   double d2 = vector3d::DotProduct(c - c4, c - c4);
   double val = (r + r4) * (r + r4) - d2;

   if (fabs(val) < err) //when the circles are touching they form a single position.
   {
    points.push_back(c + vector3d::Normalize(c4 - c) * r);
    return true;
   }
   else if (val > 0.0) // when the circles are intersecting
   {
    double val2 = (r - r4) * (r - r4) - d2;

    if (fabs(val2) <= err) // one circle is inside another and both are touching
    {
     if (r < r4) //c inside c4 and touching it
     {
      points.push_back(c + vector3d::Normalize(c - c4) * r);
     }
     else //c4 inside c and touching it
     {
      points.push_back(c4 + vector3d::Normalize(c4 - c) * r4);
     }
     return true;
    }
    else
    {
     double r42 = r4 * r4;
     double r2 = r * r;
     double check_val = r42 - (r42 - r2 + d2) * (r42 - r2 + d2) / (4.0 * d2);

     //when the circles (c4, r4) and (c, r) are intersecting they form two positions.
     if (check_val > 0)
     {
      vector3d support = vector3d::Normalize(vector3d::CrossProduct(c - c4, n4));
      double h = sqrt(check_val);
      vector3d p2 = c4 + (c - c4) * (r42 - r2 + d2) / (2.0 * d2);
      points.push_back(p2 + support * h);
      points.push_back(p2 - support * h);
      return true;
     } //when the circles (c4, r4) and (c, r) are intersecting they form two positions.
    }
   } // when the circles are intersecting
  } //must intersect the circle(c4, r4) vs circle (c,r) (both in the same plane)
 } //sphere 1 and 2 intersection generate the circle (c4, r4)
 return false;
}

/*************************************************************/

void FindRotationMatrix(vector3d va, vector3d vb, double* tr)
{
 vector3d v1;
 double k, EPS = 0.000001;
 k = vector3d::DotProduct(va, vb) + 1;
 if (k < EPS)
 {
  k = 0;
  if (fabs(va.x) > fabs(va.z))
   v1 = vector3d(-va.y, va.x, 0.0);
  else
   v1 = vector3d(0.0, -va.z, va.y);
 }
 else v1 = vector3d::CrossProduct(va, vb);
 double q_x = v1.x;
 double q_y = v1.y;
 double q_z = v1.z;
 double q_w = k;
 double l = sqrt(q_x * q_x + q_y * q_y + q_z * q_z + q_w * q_w);
 if (l == 0)
 {
  q_x = 0;
  q_y = 0;
  q_z = 0;
  q_w = 1;
 }
 else
 {
  l = 1 / l;
  q_x = q_x * l;
  q_y = q_y * l;
  q_z = q_z * l;
  q_w = q_w * l;
 }

 double x = q_x, y = q_y, z = q_z, w = q_w;
 double x2 = x + x, y2 = y + y, z2 = z + z;
 double xx = x * x2, xy = x * y2, xz = x * z2;
 double yy = y * y2, yz = y * z2, zz = z * z2;
 double wx = w * x2, wy = w * y2, wz = w * z2;

 tr[0] = 1 - (yy + zz);
 tr[4] = xy - wz;
 tr[8] = xz + wy;
 tr[1] = xy + wz;
 tr[5] = 1 - (xx + zz);
 tr[9] = yz - wx;
 tr[2] = xz - wy;
 tr[6] = yz + wx;
 tr[10] = 1 - (xx + yy);
 // last column
 tr[3] = 0;
 tr[7] = 0;
 tr[11] = 0;
 // bottom row
 tr[12] = 0;
 tr[13] = 0;
 tr[14] = 0;
 tr[15] = 1;
}

/*************************************************************/

#ifdef _WIN32
double GetTime()
{
 struct _timeb t;
#ifdef _MSC_VER
 _ftime64_s(&t);
#else
 _ftime(&t);
#endif
 double uru = ((double)t.time) + 0.001 * t.millitm;
 return uru;
}
#else
double GetTime()
{
 struct timeval t;
 gettimeofday(&t, NULL);
 double uru = ((double)t.tv_sec) + 0.000001*t.tv_usec;
 return uru;
}
#endif

}
