#ifndef _GEN_CONTAINER_H_
#define _GEN_CONTAINER_H_

#include "gen_utils.h"

namespace PG
{

/*! Pack container. */
struct Container {
 virtual void GetMin(double* m) = 0; /*!< m: pointer to array of 3 doubles. */
 virtual void GetMax(double* m) = 0; /*!< m: pointer to array of 3 doubles. */
 virtual void GetVolume(double& vol) = 0;
 /**
  * Positive distance to the closest surface point.
  * For a point p outside the container return PA::MIN_DISTANCE
  */
 virtual double GetDistance(const double* const p) = 0;
 virtual ~Container() {};
};

struct Box : public Container {
 vector3d bmin, bmax; /*!< Defined by the user. */
 inline Box();
 inline Box(vector3d bmin, vector3d bmax);
 inline void GetMin(double* m);
 inline void GetMax(double* m);
 inline void GetVolume(double& vol);
 inline double GetDistance(const double* const p);
};

struct Cylinder : public Container {
 vector3d p0, p1;/*!< Defined by the user. */
 double r;       /*!< Defined by the user. */

 vector3d bmin, bmax;   /*!< Do not modify this! */
 bool     box_computed; /*!< Do not modify this! */

 inline Cylinder();
 inline Cylinder(vector3d p0, vector3d p1, double r);
 inline void GetMin(double* m);
 inline void GetMax(double* m);
 inline void ComputeBox(); /*!< To update bmin and bmax in case p0,p1,r change. */
 inline void GetVolume(double& vol);
 inline double GetDistance(const double* const p);
};

/***************************************************************************/

Box::Box() : bmin(0,0,0), bmax(1,1,1) {}

Box::Box(vector3d bmin, vector3d bmax) : bmin(bmin), bmax(bmax) {}

void Box::GetMin(double* m) { for (int i = 0; i < 3; ++i) m[i] = bmin[i]; }

void Box::GetMax(double* m) { for (int i = 0; i < 3; ++i) m[i] = bmax[i]; }

void Box::GetVolume(double& vol)
{
 vol = 1.0;
 for (int i = 0; i < 3; ++i) vol *= (bmax[i] - bmin[i]);
}

double Box::GetDistance(const double* const p)
{
 double sqDist = MAX_DISTANCE;
 for (int i = 0; i < 3; ++i)
 {
  double v = p[i];
  if (v > bmax[i]) return MIN_DISTANCE;
  if (v < bmin[i]) return MIN_DISTANCE;
  sqDist = PA_MIN(sqDist, bmax[i] - v);
  sqDist = PA_MIN(sqDist, v - bmin[i]);
 }
 return sqDist;
}

/***************************************************************************/

void Cylinder::ComputeBox()
{
 vector3d axis_from(0, 1, 0); // base axis
 vector3d axis_to = vector3d::Normalize(p1 - p0); // current axis
 double tr[16]; // rotation matrix
 FindRotationMatrix(axis_from, axis_to, &tr[0]);

 const double half_h = (p1 - p0).length() / 2.0;
 vector3d corners[] = { vector3d(-r, -half_h, -r),
                        vector3d( r, -half_h, -r),
                        vector3d( r, -half_h,  r),
                        vector3d(-r, -half_h,  r),
                        vector3d(-r,  half_h, -r),
                        vector3d( r,  half_h, -r),
                        vector3d( r,  half_h,  r),
                        vector3d(-r,  half_h,  r) };
 const int n_corners(8);
 vector3d new_corners[n_corners];
 for (int i = 0; i < n_corners; ++i)
 {
  double x = corners[i].x, y = corners[i].y, z = corners[i].z, w = 1.0;
  new_corners[i][0] = tr[0] * x + tr[4] * y + tr[8] * z + tr[12] * w;
  new_corners[i][1] = tr[1] * x + tr[5] * y + tr[9] * z + tr[13] * w;
  new_corners[i][2] = tr[2] * x + tr[6] * y + tr[10] * z + tr[14] * w;
  corners[i] = new_corners[i];
 }
 double tt[16]; // translation matrix
 tt[0] = 1; tt[4] = 0; tt[8] = 0; tt[12] = (p1.x + p0.x) * 0.5;
 tt[1] = 0; tt[5] = 1; tt[9] = 0; tt[13] = (p1.y + p0.y) * 0.5;
 tt[2] = 0; tt[6] = 0; tt[10] = 1; tt[14] = (p1.z + p0.z) * 0.5;
 tt[3] = 0; tt[7] = 0; tt[11] = 0; tt[15] = 1;
 for (int i = 0; i < n_corners; ++i)
 {
  double x = corners[i].x, y = corners[i].y, z = corners[i].z, w = 1.0;
  new_corners[i][0] = tt[0] * x + tt[4] * y + tt[8] * z + tt[12] * w;
  new_corners[i][1] = tt[1] * x + tt[5] * y + tt[9] * z + tt[13] * w;
  new_corners[i][2] = tt[2] * x + tt[6] * y + tt[10] * z + tt[14] * w;
 }
 bmin[0] = bmin[1] = bmin[2] = MAX_DISTANCE;
 bmax[0] = bmax[1] = bmax[2] = MIN_DISTANCE;
 for (int i = 0; i < n_corners; ++i)
 {
  for (int j = 0; j < 3; ++j) bmin[j] = PA_MIN(bmin[j], new_corners[i][j]);
  for (int j = 0; j < 3; ++j) bmax[j] = PA_MAX(bmax[j], new_corners[i][j]);
 }
 box_computed = true;
}

Cylinder::Cylinder() : p0(0.0, 0.0, 0.0), p1(0.0, 2.0, 0.0), r(1.0)
{
 box_computed = false;
}

Cylinder::Cylinder(vector3d p0, vector3d p1, double r) : p0(p0), p1(p1), r(r)
{
 box_computed = false;
}

void Cylinder::GetMin(double* m)
{
 if (!box_computed) ComputeBox();
 for (int i = 0; i < 3; ++i) m[i] = bmin[i];
}

void Cylinder::GetMax(double* m)
{
 if (!box_computed) ComputeBox();
 for (int i = 0; i < 3; ++i) m[i] = bmax[i];
}

void Cylinder::GetVolume(double& vol)
{
 vol = M_PI * r * r * (p1 - p0).length();
}

//Fast distance computation between a point and cylinders, cones, line swept
//spheres and cone-spheres (http://liris.cnrs.fr/Documents/Liris-1297.pdf)
double Cylinder::GetDistance(const double* const p)
{
 const vector3d c = (p0 + p1) * 0.5;
 const double x = vector3d::DotProduct(
  /*n:*/vector3d(p[0] - c[0], p[1] - c[1], p[2] - c[2]),
  /*u:*/vector3d::Normalize(p1 - p0));
 const double l = (p1 - p0).length();
 if (fabs(x) <= l / 2.0)
 {
  const vector3d c_p(c[0] - p[0], c[1] - p[1], c[2] - p[2]);
  const double y2 = /*n2:*/vector3d::DotProduct(c_p, c_p) - x * x;
  if (y2 <= r * r) return PA_MIN(l - fabs(x), r - sqrt(y2));
  else return MIN_DISTANCE;
 }
 else return MIN_DISTANCE;
}

}

#endif /* _GEN_CONTAINER_H_ */