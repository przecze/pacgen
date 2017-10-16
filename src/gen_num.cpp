#include "gen_num.h"
#include "time.h"
#include <sstream>
#include <algorithm>

namespace PG
{

NG::NG(const NG &n) :
 m_min(n.m_min), m_max(n.m_max), m_id(n.m_id)
{
 time_t now;
 time(&now);
 srand((unsigned int)now);
}

NG::NG(double min, double max) :
 m_min(min), m_max(max), m_id("NN")
{
 time_t now;
 time(&now);
 srand((unsigned int)now);
}

NG::~NG()
{}

void NG::Reset()
{
 time_t now;
 time(&now);
 srand((unsigned int)now);
}

/***************************************************************************/

UniformNG::UniformNG(const UniformNG &n) :
 NG(n)
{}

UniformNG::UniformNG(double min, double max) :
 NG(min, max)
{ m_id = ("UF"); }

UniformNG::~UniformNG()
{
}

double UniformNG::GetValue()
{
 return m_min + ((double) rand() / (double) RAND_MAX) * (m_max - m_min);
}

std::string UniformNG::ToString()
{
 std::ostringstream s;
 s << "UF_" << m_min << ";" << m_max;
 return s.str();
}

/***************************************************************************/

BernoulliNG::BernoulliNG(const BernoulliNG &n) :
 NG(n), m_max_probability(n.m_max_probability)
{ m_id = ("BN"); }

BernoulliNG::BernoulliNG(double a, double b, double b_probability) :
 NG(a, b), m_max_probability(b_probability)
{ m_id = ("BN"); }

BernoulliNG::~BernoulliNG()
{}

double BernoulliNG::GetValue()
{
 return ((double) rand() / (double) RAND_MAX < 1.0 - m_max_probability) ? m_min : m_max;
}

std::string BernoulliNG::ToString()
{
 std::ostringstream s;
 s << "BN_" << m_min << ";" << m_max << ";" << m_max_probability;
 return s.str();
}

/***************************************************************************/

GaussianNG::GaussianNG(const GaussianNG &n) :
 NG(n), m_mean(n.m_mean), m_sd(n.m_sd)
{ m_id = ("GS"); }

GaussianNG::GaussianNG(double min, double max, double mean, double sd) :
 NG(min, max), m_mean(mean), m_sd(sd)
{ m_id = ("GS"); }

GaussianNG::~GaussianNG()
{}

double GaussianNG::BoxMuller(double mean, double stddev)
{
#if 0
 return mean + generateGaussianNoise(1.0) * stddev;
#endif
 static double n2 = 0.0;
 static int n2_cached = 0;
 if ( !n2_cached )
 {
  double x, y, r;
  do
  {
   x = 2.0*rand()/RAND_MAX - 1;
   y = 2.0*rand()/RAND_MAX - 1;
   r = x*x + y*y;
  } while (r == 0.0 || r > 1.0);
  {
   double d = sqrt(-2.0*log(r)/r);
   double n1 = x*d;
   n2 = y*d;
   double result = n1*stddev + mean;
   n2_cached = 1;
   return result;
  }
 }
 else
 {
  n2_cached = 0;
  return n2*stddev + mean;
 }
}

double GaussianNG::GetValue()
{
 double val = BoxMuller(m_mean, m_sd);
 while( val < m_min || val > m_max )
  val = BoxMuller(m_mean, m_sd);
 return val;
}

std::string GaussianNG::ToString()
{
 std::ostringstream s;
 s << "GS_" << m_min << ";" << m_max << ";" << m_mean << ";" << m_sd;
 return s.str();
}

/***************************************************************************/

GeneralNG::GeneralNG(const double n[], const double p[], size_t size) :
 NG(n[0], n[size - 1])
{
 m_id = ("GG");
 m_size = size;
 m_n    = new double[m_size];
 for(size_t i = 0; i < m_size; ++i) m_n[i] = n[i];
 m_p    = new double[m_size];
 for(size_t i = 0; i < m_size; ++i) m_p[i] = p[i];
 m_cp   = new double[m_size + 1];
 m_cp[0]= 0.0;
 for(size_t i = 0; i < m_size; ++i) m_cp[i+1] = p[i] + m_cp[i];
 if ( fabs(m_cp[m_size] - 1.0) >= 1.e-6 ) printf("The probabilities do not sum 1.0 \n");
}

GeneralNG::GeneralNG(const GeneralNG &n) :
 NG(n)
{
 m_id = ("GG");
 m_size = n.m_size;
 m_n    = new double[m_size];
 for(size_t i = 0; i < m_size; ++i) m_n[i] = n.m_n[i];
 m_p    = new double[m_size];
 for(size_t i = 0; i < m_size; ++i) m_p[i] = n.m_p[i];
 m_cp   = new double[m_size + 1];
 m_cp[0]= 0.0;
 for(size_t i = 0; i < m_size; ++i) m_cp[i+1] = n.m_p[i] + m_cp[i];
 if ( fabs(m_cp[m_size] - 1.0) >= 1.e-6 ) printf("The probabilities do not sum 1.0 \n");
}

GeneralNG::~GeneralNG()
{
 delete[] m_n;
 delete[] m_p;
 delete[] m_cp;
}

double GeneralNG::GetValue()
{
 double r = (double) rand() / (double) RAND_MAX;
 for(size_t i = 0; i < m_size; ++i)
  if ( r >= m_cp[i] && r < m_cp[i+1] ) return m_n[i];
 return m_n[m_size - 1];
}

std::string GeneralNG::ToString()
{
 std::ostringstream s;
 return s.str();
}

/***************************************************************************/

/**
 * n has 'size' elements defining the lower and upper limits of each range
 * p has 'size-1' frequencies
 */
RangesNG::RangesNG(const double n[], const double p[], size_t size) :
 NG(n[0], n[size - 1])
{
 m_id = ("RG");
 m_size = size;
 m_n    = new double[m_size];
 for(size_t i = 0; i < m_size; ++i) m_n[i] = n[i];
 m_p    = new double[m_size - 1];
 for(size_t i = 0; i < m_size - 1; ++i) m_p[i] = p[i];
 m_cp   = new double[m_size];
 m_cp[0]= 0.0;
 for(size_t i = 0; i < m_size - 1; ++i) m_cp[i+1] = p[i] + m_cp[i];
 if ( fabs(m_cp[m_size-1] - 1.0) >= 1.e-6 ) printf("The probabilities do not sum 1.0 \n");
}

RangesNG::RangesNG(const RangesNG &n) :
 NG(n)
{
 m_id = ("RG");
 m_size = n.m_size;
 m_n    = new double[m_size];
 for(size_t i = 0; i < m_size; ++i) m_n[i] = n.m_n[i];
 m_p    = new double[m_size - 1];
 for(size_t i = 0; i < m_size - 1; ++i) m_p[i] = n.m_p[i];
 m_cp   = new double[m_size];
 m_cp[0]= 0.0;
 for(size_t i = 0; i < m_size - 1; ++i) m_cp[i+1] = n.m_p[i] + m_cp[i];
 if ( fabs(m_cp[m_size-1] - 1.0) >= 1.e-6 ) printf("The probabilities do not sum 1.0 \n");
}

RangesNG::~RangesNG()
{
 delete[] m_n;
 delete[] m_p;
 delete[] m_cp;
}

double RangesNG::GetValue()
{
 static double gn;
 static int gn_cached = 0;
 if ( !gn_cached ) gn = (double) rand() / (double) RAND_MAX;
 else gn_cached = 0;
 for(size_t i = 0; i < m_size-1; ++i)
  if ( gn >= m_cp[i] && gn < m_cp[i+1] )
  {
   gn = ((double) rand() / (double) RAND_MAX);
   gn_cached = 1;
   return m_n[i] + gn * (m_n[i + 1] - m_n[i]);
  }
 gn = ((double) rand() / (double) RAND_MAX);
 gn_cached = 1;
 return m_n[m_size - 2] + gn * (m_n[m_size - 1] - m_n[m_size - 2]);
}

std::string RangesNG::ToString()
{
 std::ostringstream s;
 return s.str();
}

}