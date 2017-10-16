#include "gen_pack.h"
#include <fstream>
#include <queue>
#include <time.h>
#include <stdio.h>
#include <iomanip>

namespace PG
{
//both thresholds are applied to the sq of the radii sum
#define CNTC_THRES 0.002500 // ( 0.5 / 10 ) ^ 2
#define COL_THRES  0.000625 // ( 0.5 / 20 ) ^ 2
#define _HANDLES_REJECTION_

/***************************************************************************/

Sphere::Sphere(double p_x, double p_y, double p_z, double p_r) :
 x(p_x), y(p_y), z(p_z), r(p_r) {}

/***************************************************************************/

double GetNextRadius(NG* const ng, std::queue<double>& lprs)
{
#ifdef _HANDLES_REJECTION_
 if (!lprs.empty())
 {
  const double r(lprs.front());
  lprs.pop();
  return r;
 }
#endif
 return ng->GetValue();
};

void Insert(Grid3d* const dom, SpherePack* const pack,
 std::queue<uint32>& fidQ, const vector3d& p, const double& rad,
 SpherePackStat* const stat)
{
 pack->s.push_back(Sphere(p.x, p.y, p.z, rad));
 int32 dminx, dminy, dminz, dmaxx, dmaxy, dmaxz;
 dminx = static_cast<int32>(floor((p.x - rad - dom->emin[0]) * dom->invsvoxel[0]));
 dminy = static_cast<int32>(floor((p.y - rad - dom->emin[1]) * dom->invsvoxel[1]));
 dminz = static_cast<int32>(floor((p.z - rad - dom->emin[2]) * dom->invsvoxel[2]));
 dmaxx = static_cast<int32>(floor((p.x + rad - dom->emin[0]) * dom->invsvoxel[0]));
 dmaxy = static_cast<int32>(floor((p.y + rad - dom->emin[1]) * dom->invsvoxel[1]));
 dmaxz = static_cast<int32>(floor((p.z + rad - dom->emin[2]) * dom->invsvoxel[2]));
 for (int32 i = dminx; i <= dmaxx; ++i)
  for (int32 j = dminy; j <= dmaxy; ++j)
   for (int32 k = dminz; k <= dmaxz; ++k)
   {
    int32 c = (i + j * dom->nvoxel[0] + k * dom->nvoxel[0] * dom->nvoxel[1]);
    if (0 <= c && c < (int32) dom->totalVoxels)
     dom->voxels[c].insert(stat->nspheres);
   }
   fidQ.push(stat->nspheres++);
};

double my_rand(const double& min, const double& max)
{
 return min + ((double) rand() / (double) RAND_MAX) * (max - min);
};

SpherePackStat GenerateSpherePack(
 Container* const df, NG* const ng, Grid3d* const dom, SpherePack* const pack)
{
 if (!(df && ng && pack)) return SpherePackStat();

 SpherePackStat     stat;
 uint32             fid(100); // front id
 std::queue<uint32> fidQ;     // queue of front ids
 std::set<uint32>   closerFronts;
#ifdef _HANDLES_REJECTION_
 std::queue<double> lnrs; /*!< Queue of newly rejected sphere sizes.      */
 std::queue<double> lprs; /*!< Queue of previously rejected sphere sizes. */
#endif
 double             contactThreshold = ng->m_min * 1.e-4;

 /*************** INIT: UNIFORM GRID INITIALIZATION ****************/
 {
  const double voxel_size = 2.0 * ng->m_max;
  df->GetMin(&dom->emin[0]);
  df->GetMax(&dom->emax[0]);
  for (int i = 0; i < 3; ++i) dom->invsvoxel[i] = 1.0/voxel_size;
  for (int i = 0; i < 3; ++i)
   dom->nvoxel[i] = (int) (ceil(fabs(dom->emax[i] - dom->emin[i])/voxel_size));
  dom->totalVoxels = dom->nvoxel[0] * dom->nvoxel[1] * dom->nvoxel[2];
  dom->voxels.clear();
  dom->voxels.resize(dom->totalVoxels);
 }
 /*************** END: UNIFORM GRID INITIALIZATION ****************/

 stat.time = GetTime();
 stat.nspheres = 0;

 /*************** INIT: INSERT THE FIRST THREE SPHERES ****************/
 {
  time_t now;
  time(&now);
  srand((unsigned int)now);
  vector3d pos1((dom->emin[0] + dom->emax[0])*0.5,
                (dom->emin[1] + dom->emax[1])*0.5,
                (dom->emin[2] + dom->emax[2])*0.5);
  const double r1 = GetNextRadius(ng, lprs);
  uint32 exit_loop(0);
  while (true)
  {
   if (df->GetDistance(&pos1.x) < 2.5 * r1)
   {
    for (int i = 0; i < 3; ++i)
     pos1[i] = my_rand(dom->emin[i], dom->emax[i]);
   }
   else break;
   if (exit_loop++ > 10000)
   {
    stat.time = 0.0;
    printf("ERROR GENERATING THE SEEDS \n");
    return stat;
   }
  }
  Insert(dom, pack, fidQ, pos1, r1, &stat); // first

  vector3d dir(my_rand(-1.0, 1.0), my_rand(-1.0, 1.0), my_rand(-1.0, 1.0));
  dir = vector3d::Normalize(dir);
  const double r2 = GetNextRadius(ng, lprs);
  const vector3d pos2 = pos1 + dir * (r2 + pack->s[0].r);
  Insert(dom, pack, fidQ, pos2, r2, &stat); // second

  std::vector<vector3d> spoints;
  const double r3 = GetNextRadius(ng, lprs);
  FindIntersectionPoints(contactThreshold, &pos1[0], r1 + r3,
                         &pos2[0], r2 + r3, spoints, 72);
  const vector3d pos3 = spoints[rand() % 72];
  Insert(dom, pack, fidQ, pos3, r3, &stat); // third
 }
 /*************** END: INSERT THE FIRST THREE SPHERES ****************/

 if (stat.nspheres < 3) return SpherePackStat();

 while (!fidQ.empty())
 {
  bool same_front = (fid == fidQ.front());
  fid = fidQ.front();
  const double new_rad(GetNextRadius(ng, lprs));

  if (!same_front)
  {
   const Sphere s = pack->s[fid];
   const double size = (s.r + 2.1 * new_rad);
   const double minx = s.x - size;
   const double miny = s.y - size;
   const double minz = s.z - size;
   const double maxx = s.x + size;
   const double maxy = s.y + size;
   const double maxz = s.z + size;

   closerFronts.clear();
   int32 dminx, dminy, dminz, dmaxx, dmaxy, dmaxz;
   dminx = static_cast<int32>(floor((minx - dom->emin[0]) * dom->invsvoxel[0]));
   dminy = static_cast<int32>(floor((miny - dom->emin[1]) * dom->invsvoxel[1]));
   dminz = static_cast<int32>(floor((minz - dom->emin[2]) * dom->invsvoxel[2]));
   dmaxx = static_cast<int32>(floor((maxx - dom->emin[0]) * dom->invsvoxel[0]));
   dmaxy = static_cast<int32>(floor((maxy - dom->emin[1]) * dom->invsvoxel[1]));
   dmaxz = static_cast<int32>(floor((maxz - dom->emin[2]) * dom->invsvoxel[2]));
   for (int32 i = dminx; i <= dmaxx; ++i)
    for (int32 j = dminy; j <= dmaxy; ++j)
     for (int32 k = dminz; k <= dmaxz; ++k)
     {
      int32 cc = (i + j * dom->nvoxel[0] + k * dom->nvoxel[0] * dom->nvoxel[1]);
      if (0 <= cc && cc < (int32) dom->totalVoxels) // NOTE: valid voxel id
      {
       std::set<uint32>::iterator tid;
       for (tid = dom->voxels[cc].begin(); tid != dom->voxels[cc].end(); ++tid)
       {
        //check if the AABB the neighbor collides the the front AABB
        const double nminx = pack->s[*tid].x - pack->s[*tid].r;
        const double nminy = pack->s[*tid].y - pack->s[*tid].r;
        const double nminz = pack->s[*tid].z - pack->s[*tid].r;
        const double nmaxx = pack->s[*tid].x + pack->s[*tid].r;
        const double nmaxy = pack->s[*tid].y + pack->s[*tid].r;
        const double nmaxz = pack->s[*tid].z + pack->s[*tid].r;

        bool collides(false);
        if (minx >= nminx && maxx <= nmaxx &&
            miny >= nminy && maxy <= nmaxy &&
            minz >= nminz && maxz <= nmaxz)
         collides = true;
        else
        {
         // Exit with no intersection if separated along an axis
         if (nmaxx < minx || nminx > maxx) collides = false;
         else
         {
          if (nmaxy < miny || nminy > maxy) collides = false;
          else
          {
           if (nmaxz < minz || nminz > maxz) collides = false;
           else collides = true;
          }
         }
        }
        if (collides)
         closerFronts.insert(*tid);
       }
      }
     }
  }
  else
   closerFronts.insert(stat.nspheres - 1);

  std::vector<vector3d> points;
  if (!closerFronts.empty())
  {
   std::set<uint32>::iterator ii, beforeLast = closerFronts.end(); --beforeLast;
   for (ii = closerFronts.begin(); ii != beforeLast; ++ii)
   {
    if (fid == *ii) continue;
    std::set<uint32>::iterator jj = ii; ++jj;
    for (; jj != closerFronts.end(); ++jj)
    {
     if (fid == *jj) continue;
     const vector3d c3(pack->s[*jj].x, pack->s[*jj].y, pack->s[*jj].z);
     FindIntersectionPoints(contactThreshold,
      &pack->s[fid].x, pack->s[fid].r + new_rad,
      &pack->s[*ii].x, pack->s[*ii].r + new_rad,
      c3             , pack->s[*jj].r + new_rad, points);
    }
   }
  }

  const uint16 n_points = (uint16) points.size();
  uint16 valids(0);
  vector3d best(MAX_DISTANCE, MAX_DISTANCE, MAX_DISTANCE);

  for (int i = 0; i < n_points; ++i)
  {
   const vector3d p = points[i];
   if (df->GetDistance(&p.x) >= new_rad)
   {
    bool discarded(false);
    std::set<uint32>::iterator it;
    for (it = closerFronts.begin(); it != closerFronts.end(); ++it)
    {
     const Sphere s = pack->s[*it];
     const vector3d dd(p.x - s.x, p.y - s.y, p.z - s.z);
     const double sq_dist(vector3d::DotProduct(dd, dd));
     const double sq_sum = (s.r + new_rad) * (s.r + new_rad);
     if (sq_sum - sq_dist > sq_sum * COL_THRES)
     {
      discarded = true;
      break;
     }
    }
    if (!discarded)
    {
     if (valids == 0) best = p;
     else
     {
      const vector3d first(pack->s[0].x, pack->s[0].y, pack->s[0].z);
      if (vector3d::DotProduct(first - p, first - p) < vector3d::DotProduct(first - best, first - best))
       best = p;
     }
     valids++;
    }
   }
  }

  if (valids < 2) fidQ.pop();

  if (valids > 0)
  {
   Insert(dom, pack, fidQ, best, new_rad, &stat);
#ifdef _HANDLES_REJECTION_
   while (!lnrs.empty())
   {
    lprs.push(lnrs.front());
    lnrs.pop();
   }
#endif
  }
#ifdef _HANDLES_REJECTION_
  else lnrs.push(new_rad);
#endif
 }
 stat.time = GetTime() - stat.time;
 return stat;
}

void SpherePackStatistics(
 Container* const df, Grid3d* const dom, SpherePack* const pack,
 SpherePackStat* const stat)
{
 stat->vf = 0.0;
 const size_t n_spheres(pack->s.size());
 for (size_t i = 0; i < n_spheres; ++i) stat->vf += pow(pack->s[i].r, 3.0);
 stat->vf = stat->vf * (4.0 * M_PI / 3.0);
 double cdom_vol(1.0);
 df->GetVolume(cdom_vol);
 stat->vf /= cdom_vol;
 stat->porosity = 1.0 - stat->vf;

 stat->cn.reserve(n_spheres);
 stat->cn.resize(n_spheres);
 stat->mcn = 0.0;

 for (uint32 pid = 0; pid < n_spheres; ++pid)
 {
  std::set<uint32> neighbors;
  const Sphere s = pack->s[pid];
  const double size = 1.5 * s.r; // NOTE: a box twice the size of the sphere
  const double minx = s.x - size;
  const double miny = s.y - size;
  const double minz = s.z - size;
  const double maxx = s.x + size;
  const double maxy = s.y + size;
  const double maxz = s.z + size;
  int32 dminx, dminy, dminz, dmaxx, dmaxy, dmaxz;
  dminx = static_cast<int32>(floor((minx - dom->emin[0]) * dom->invsvoxel[0]));
  dminy = static_cast<int32>(floor((miny - dom->emin[1]) * dom->invsvoxel[1]));
  dminz = static_cast<int32>(floor((minz - dom->emin[2]) * dom->invsvoxel[2]));
  dmaxx = static_cast<int32>(floor((maxx - dom->emin[0]) * dom->invsvoxel[0]));
  dmaxy = static_cast<int32>(floor((maxy - dom->emin[1]) * dom->invsvoxel[1]));
  dmaxz = static_cast<int32>(floor((maxz - dom->emin[2]) * dom->invsvoxel[2]));
  for (int32 i = dminx; i <= dmaxx; ++i)
   for (int32 j = dminy; j <= dmaxy; ++j)
    for (int32 k = dminz; k <= dmaxz; ++k)
    {
     int32 cc = (i + j * dom->nvoxel[0] + k * dom->nvoxel[0] * dom->nvoxel[1]);
      if (0 <= cc && cc < (int32) dom->totalVoxels) // NOTE: valid voxel id
      {
       std::set<uint32>::iterator tid;
       for (tid = dom->voxels[cc].begin(); tid != dom->voxels[cc].end(); ++tid)
       {
        if (*tid == pid) continue; // NOTE: If it is the same sphere
        const Sphere neig = pack->s[*tid];

        //check if the AABB the neighbor collides the the current sphere's AABB
        const double nminx = neig.x - neig.r;
        const double nminy = neig.y - neig.r;
        const double nminz = neig.z - neig.r;
        const double nmaxx = neig.x + neig.r;
        const double nmaxy = neig.y + neig.r;
        const double nmaxz = neig.z + neig.r;

        const double fminx = s.x - s.r;
        const double fminy = s.y - s.r;
        const double fminz = s.z - s.r;
        const double fmaxx = s.x + s.r;
        const double fmaxy = s.y + s.r;
        const double fmaxz = s.z + s.r;

        bool collides(false);
        if (fminx >= nminx && fmaxx <= nmaxx &&
            fminy >= nminy && fmaxy <= nmaxy &&
            fminz >= nminz && fmaxz <= nmaxz)
         collides = true;
        else
        {
         // Exit with no intersection if separated along an axis
         if (nmaxx < fminx || nminx > fmaxx) collides = false;
         else
         {
          if (nmaxy < fminy || nminy > fmaxy) collides = false;
          else
          {
           if (nmaxz < fminz || nminz > fmaxz) collides = false;
           else collides = true;
          }
         }
        }
        if (collides)
        {
         vector3d dd(s.x - neig.x, s.y - neig.y, s.z - neig.z);
         double sq_dist = vector3d::DotProduct(dd, dd);
         double sq_sum = (neig.r + s.r) * (neig.r + s.r);
         if (sq_dist - sq_sum <= CNTC_THRES * sq_sum)
          neighbors.insert(*tid);
        }
       }
      }
    }
  stat->cn[pid] = (uint16) neighbors.size();
  stat->mcn += (double) stat->cn[pid];
 }
 stat->mcn /= (double) n_spheres;
}

}