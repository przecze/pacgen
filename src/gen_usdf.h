#ifndef _GEN_USDF_H_
#define _GEN_USDF_H_

#include "gen_container.h"

namespace PG
{

struct USDFP /*! Uniform Signed Distance Field Point. */
{
 vector3d pos; /*!< Point position. */
 double   dist;/*!< Distance.       */
 USDFP(const USDFP& s);
 USDFP(double x, double y, double z, double dist);
 bool operator<(const USDFP& other) const;
};

class USDF : public Container /*! Uniform Signed Distance Field. */
{
public:
 std::vector<USDFP*> points;  /*!< List of distance points.                 */
 void*     container;     /*!< The mesh boundary.                           */
 vector3d  points_min;    /*!< AABB of distance pointa.                     */
 vector3d  points_max;    /*!<                                              */
 vector3d  steps;         /*!< Delta step size between distance points.     */
 vector3d  inv_steps;     /*!< Inverse delta distances.                     */
 double    pos_epsilon;   /*!< Displacement factor for the outer shell.     */
 double    neg_epsilon;   /*!< Displacement factor for the inner shell.     */
 double    num_xz;        /*!< Num point in x times the num points in z     */
 int       num_points[3]; /*!< Number of distance points per dimension.     */
 int       total_points;  /*!< Total number of points in the distance field.*/

public:
 /**
  * Creates a uniform grid of points.
  * Computes the exact distance for points inside an inner and outer shell.
  * @param mesh Triangle mesh.
  * @param rmax Maximum radius in the pack.
  * @param npoints Number of distance points in each dimension.
  */
 USDF(void* mesh, const double& rmax, const int* const npoints);
 ~USDF();
 /**
  * Interpolates the distance of the point p using the grid of points.
  * @param p A pointer to the query point.
  * @return The interpolated distance from p to the triangle surface.
  */
 double GetDistance(const double* const p);
 void GetMin(double* m);
 void GetMax(double* m);
 void GetVolume(double& vol);

private:
 std::vector<USDFP*>
  getIntersectingPoints(const vector3d& bmin, const vector3d& bmax);

 void initial_distances();
 void compute_signed_distances();
};

}
#endif /* _GEN_USDF_H_ */