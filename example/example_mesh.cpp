#include "gen_pack.h"
#include "gen_usdf.h"
#include "mesh.h"

/**
 * Creates a sphere triangle mesh.
 */
void* CreateSphereMesh(double radius, int slices, int stacks)
{
 std::vector<double> v; // vertices
 for (int j = 0; j <= slices; ++j)
 {
  double theta = static_cast<double>(j) * 2.0 * M_PI / static_cast<double>(slices);
  double sinTheta = sin(theta);
  double cosTheta = cos(theta);
  for (int i = 0; i <= stacks; ++i)
  {
   double phi = static_cast<double>(i) * M_PI / static_cast<double>(stacks);
   double sinPhi = sin(phi);
   double cosPhi = cos(phi);
   v.push_back(/*x*/radius * cosTheta * sinPhi);
   v.push_back(/*y*/radius * cosPhi);
   v.push_back(/*z*/radius * sinTheta * sinPhi);
  }
 }
 auto Index = [] (int i, int j, int slices) {return i + j * (slices + 1);};
 std::vector<int> idx; // indices
 for (int j = 0; j < slices; ++j)
 {
  for (int i = 0; i < stacks; ++i)
  {
   idx.push_back(Index(i + 1, j + 1, stacks));
   idx.push_back(Index(i, j, stacks));
   idx.push_back(Index(i, j + 1, stacks));

   idx.push_back(Index(i + 1, j, stacks));
   idx.push_back(Index(i, j, stacks));
   idx.push_back(Index(i + 1, j + 1, stacks));
  }
 }
 return PG::Mesh3DCreate(&v[0], v.size()/3, &idx[0], idx.size());
}


int main(int argc, char **argv)
{
 void* mesh = CreateSphereMesh(2.0, 50, 50);

 double rmin(0.02), rmax(0.05);
 PG::NG* ng = new PG::UniformNG(rmin, rmax);
 const int npoints[] = {100, 100, 100};
 PG::Container* container = new PG::USDF(mesh, rmax, npoints);

 PG::Grid3d dom;
 PG::SpherePack* pack = new PG::SpherePack();
 PG::SpherePackStat result = PG::GenerateSpherePack(container, ng, &dom, pack);

 delete pack;
 delete ng;
 delete container;
 PG::Mesh3DDestroy(mesh);

 return 0;
}
