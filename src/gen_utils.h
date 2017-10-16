#ifndef _GEN_UTILS_H_
#define _GEN_UTILS_H_

#include <vector>
#include "vectors.h"
#include <limits>
#include <stdint.h>

namespace PG
{

typedef unsigned short uint16;
typedef unsigned int   uint32;
typedef short          int16;
typedef int            int32;

static double MAX_DISTANCE = std::numeric_limits<double>::infinity();
static double MIN_DISTANCE = -std::numeric_limits<double>::infinity();

/**
 * For two spheres
 */
bool FindIntersectionPoints(const double& err, const double* c1,
                            const double& r1, const double* c2,
                            const double& r2, std::vector<vector3d>& points,
                            const unsigned int& p);

/**
 * Between two spheres
 */
bool FindIntersectingCircle(const double& err, const double* c1,
                            const double& r1, const double* c2,
                            const double& r2, double* c3, double& r3,
                            double* n3, double* u, double* v,
                            bool& special);

/**
 * For three spheres
 */
bool FindIntersectionPoints(const double& err, const double* c1,
                            const double& r1, const double* c2,
                            const double& r2, const vector3d& c3,
                            const double& r3, std::vector<vector3d>& points);


/****************************************************************************/
//4x4 Rotation matrix to rotate from va to vb
void FindRotationMatrix(vector3d va, vector3d vb, double* tr);

/****************************************************************************/

double GetTime();

}
#endif /* _GEN_UTILS_H_ */