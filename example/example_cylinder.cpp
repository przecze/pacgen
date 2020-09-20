#include <fstream>
#include <iomanip>

#include "gen_pack.h"

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
 PG::Container* container = new PG::Cylinder({0.0, 0.0, 0.0},
                                             {0.0, 5.0, 0.0},
                                             1.0);
 double rmin(0.01), rmax(0.03), rmaxprob(0.5);
 PG::NG* ng = new PG::BernoulliNG(rmin, rmax, rmaxprob);

 PG::Grid3d dom;
 PG::SpherePack* pack = new PG::SpherePack();
 PG::SpherePackStat result = PG::GenerateSpherePack(container, ng, &dom, pack);
 ExportTxt("cylinder.txt", pack);

 delete pack;
 delete ng;
 delete container;

 return 0;
}
