#ifndef _GEN_PACK_H_
#define _GEN_PACK_H_

#include <set>
#include "gen_container.h"
#include "gen_num.h"

namespace PG
{

struct Sphere {
 double x,y,z,r; /*!< Position and radius. */
 Sphere(double x, double y, double z, double r);
};

struct SpherePack {
 std::vector<Sphere> s;
};

/*! Pack statistics */
struct SpherePackStat {
 double mcn;             /*!< mean coordination number.         */
 double vf;              /*!< volume fraction.                  */
 double porosity;        /*!< porosity.                         */
 double time;            /*!< generation time in seconds.       */
 uint32 nspheres;        /*!< number of generated spheres.      */
 std::vector<uint16> cn; /*!< per particle coordination number. */
};

/*! Uniform grid. Domain subdivision to accelerate collision detection */
struct Grid3d {
 double emin[3];         /*!< Left bottom corner.             */
 double emax[3];         /*!< Right upper corner.             */
 double invsvoxel[3];    /*!< Inverse size of each voxel.     */
 uint32 totalVoxels;     /*!< Total number of voxels.         */
 std::vector<std::set<uint32> > voxels;
 uint16 nvoxel[3];       /*!< Num of voxels on each dimension.*/
};

/***************************************************************************/

/**
 * Fills a SpherePack with particles inside a Container following a distribution.
 * @param df   Input: Distance field container.
 * @param ng   Input: Number generator.
 * @param dom  Input: Grid3d
 * @param pack Input/Output: Pack.
 * @return SpherePackStat: Number of particles and generation time in seconds.
 */
SpherePackStat GenerateSpherePack(
 Container* const df, NG* const ng, Grid3d* const dom, SpherePack* const pack);

/**
 * Complements SpherePackStatistical values of interest.
 * @param df   Input: Distance field container.
 * @param dom  Input: Grid3d
 * @param pack Input: Pack.
 * @return SpherePackStat: Porosity, volume fraction and
 *                         coordination number(mean and per particle).
 */
void SpherePackStatistics(
 Container* const df, Grid3d* const dom, SpherePack* const pack,
 SpherePackStat* const SpherePackStat);

}

#endif /* _GEN_PACK_H_ */