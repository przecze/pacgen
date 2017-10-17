#include <fstream>
#include <iomanip>
#include "src/gen_pack.h"

void ExportTxt(const char* filename, PG::SpherePack* const pack)
{
 remove(filename);
 std::ofstream myfile;
 myfile.open(filename);
 myfile << std::setprecision(15);
 const size_t n_spheres(pack->s.size());
 for (PG::Sphere s : pack->s)
 {
  myfile << s.x << " " << s.y << " " << s.z << " " << s.r << "\n";
 }
 myfile.close();
}

int main(int argc, char **argv)
{
 PG::Container* container = new PG::Box({-1.0, -1.0, -1.0},
                                        { 1.0,  1.0,  1.0});
 double rmin(0.02), rmax(0.05);
 PG::NG* ng = new PG::UniformNG(rmin, rmax);

 PG::Grid3d dom;
 PG::SpherePack* pack = new PG::SpherePack();
 PG::SpherePackStat result = PG::GenerateSpherePack(container, ng, &dom, pack);
 ExportTxt("box.txt", pack);

 delete pack;
 delete ng;
 delete container;

 return 0;
}