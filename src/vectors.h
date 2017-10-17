#ifndef _VECTORS_H_
#define _VECTORS_H_

#include <assert.h>
#define _USE_MATH_DEFINES
#include <math.h>

namespace PG
{
 #define PG_MIN(x,y)   ((x) < (y) ? (x) : (y))
 #define PG_MAX(x,y)   ((x) > (y) ? (x) : (y))
 #define PG_DEG2RAD(x) ((x * M_PI) / 180.0)
 #define PG_RAD2DEG(x) ((x * 180.0) /M_PI);

 /*********************************************************************/

 template<class T>
 struct vec3
 {
  T x, y, z;

  vec3() : x(0), y(0), z(0) {}

  template<class FromT>
  vec3(const FromT& x, const FromT& y, const FromT& z)
   : x(x), y(y), z(z){}

  template<class FromT>
  vec3<T>& operator=(const vec3<FromT>& t)
  {
   if (this == &t)
    return *this;
   x = t.x;
   y = t.y;
   z = t.z;
   return *this;
  }
  vec3<T>& operator +=(const vec3<T>& t)
  {
   x += t.x;
   y += t.y;
   z += t.z;
   return *this;
  }
  vec3<T> operator +(const vec3<T>& t) const
  {
   return vec3<T>(x + t.x, y + t.y, z + t.z);
  }
  vec3<T>& operator -=(const vec3<T>& t)
  {
   x -= t.x;
   y -= t.y;
   z -= t.z;
   return *this;
  }
  vec3<T> operator -(const vec3<T>& t) const
  {
   return vec3<T>(x - t.x, y - t.y, z - t.z);
  }
  vec3<T> operator -() const
  {
   return vec3<T>(-x, -y, -z);
  }
  vec3<T> operator *(double t) const
  {
   return vec3<T>(x * t, y * t, z * t);
  }
  vec3<T> operator /(double t) const
  {
   return vec3<T>(x / t, y / t, z / t);
  }
  T length() const
  {
   return sqrt(x * x + y * y + z * z);
  }
  T & operator[](int n)
  {
   assert(n >= 0 && n <= 2);
   if      (n == 0) return x;
   else if (n == 1) return y;
   else             return z;
  }
  const T & operator[](int n) const
  {
   assert(n >= 0 && n <= 2);
   if      (n == 0) return x;
   else if (n == 1) return y;
   else             return z;
  }
  bool const operator ==(const vec3<T>& t) const
  {
   if (fabs(y - t.y) > 1.e-15) return false;
   if (fabs(z - t.z) > 1.e-15) return false;
   if (fabs(x - t.x) > 1.e-15) return false;
   return true;
  }
  bool const operator <(const vec3<T>& t) const
  {
   if (fabs(y - t.y) > 1.e-15) return (y < t.y);
   if (fabs(z - t.z) > 1.e-15) return (z < t.z);
   if (fabs(x - t.x) > 1.e-15) return (x < t.x);
   return false;
  }

  static vec3<T> Normalize(const vec3<T>& v)
  {
   T l = v.length();
   if (l == 0.0)
    return vec3<T>(0.0, 0.0, 0.0);
   return vec3<T>(v.x / l, v.y / l, v.z / l);
  }
  static void Normalize(T* v)
  {
   T l = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
   for (int i = 0; i < 3; ++i)
    v[i] = (l == 0.0) ? 0.0 : (v[i] / l);
  }
  static vec3<T> CrossProduct(const vec3<T>& v, const vec3<T>& rhs)
  {
   vec3<T> r;
   r[0] = v[1] * rhs[2] - rhs[1] * v[2];
   r[1] = v[2] * rhs[0] - rhs[2] * v[0];
   r[2] = v[0] * rhs[1] - rhs[0] * v[1];
   return r;
  }
  static void CrossProduct(T* v, T* rhs, T* res)
  {
   res[0] = v[1] * rhs[2] - rhs[1] * v[2];
   res[1] = v[2] * rhs[0] - rhs[2] * v[0];
   res[2] = v[0] * rhs[1] - rhs[0] * v[1];
  }
  static T DotProduct(const vec3<T>& v, const vec3<T>& t)
  {
   return v.x * t.x + v.y * t.y + v.z * t.z;
  }
  static T DotProduct(T* v, T* t)
  {
   return v[0] * t[0] + v[1] * t[1] + v[2] * t[2];
  }

 };

 /*********************************************************************/

 typedef vec3<double> vector3d;
}

#endif /* _VECTORS_H_ */