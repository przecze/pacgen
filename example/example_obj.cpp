#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "src/gen_pack.h"
#include "src/gen_usdf.h"
#include "src/mesh.h"

/**
 * A very simple obj reader.
 */
void* CreateMeshFromFile(const char* filename)
{
 std::ifstream in(filename, std::ios::in);
 if (!in) return nullptr;

 std::vector<double> v; // vertices
 std::vector<int> idx; // indices
 std::string line;
 while (getline(in, line))
 {
  if (line.substr(0, 2) == "v ") // vertex
  {
   std::istringstream s(line.substr(2));
   double x,y,z;
   s >> x;
   s >> y;
   s >> z;
   v.push_back(x);
   v.push_back(y);
   v.push_back(z);
  }
  else if (line.substr(0, 2) == "f ") // face
  {
   int a, b, c;
   sscanf(line.c_str(), "f %d %d %d", &a, &b, &c);
   idx.push_back(a-1);
   idx.push_back(b-1);
   idx.push_back(c-1);
  }
 }

 return PG::Mesh3DCreate(&v[0], v.size()/3, &idx[0], idx.size());
}

int main(int argc, char **argv)
{
 void* mesh = CreateMeshFromFile("capsule.obj");

 double rmin(0.02), rmax(0.06);
 PG::NG* ng = new PG::UniformNG(rmin, rmax);
 const int npoints[] = {160, 88, 88};
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