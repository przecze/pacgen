#ifndef _GEN_NUM_H_
#define _GEN_NUM_H_

#include <string>

namespace PG
{

class NG
{
public:
 double m_min, m_max;
 std::string m_id;
public:
 NG(const NG &n);
 NG(double min, double max);
 virtual ~NG();
 void Reset();
 virtual double GetValue() = 0;
 virtual std::string ToString() = 0;
};

class UniformNG : public NG
{
public:
 UniformNG(const UniformNG &n);
 UniformNG(double min, double max);
 ~UniformNG();
 double GetValue();
 std::string ToString();
};

class BernoulliNG : public NG
{
public:
 double m_max_probability;
 BernoulliNG(const BernoulliNG &n);
 BernoulliNG(double min, double max, double max_probability);
 ~BernoulliNG();
 double GetValue();
 std::string ToString();
};

class GaussianNG : public NG
{
public:
 double m_mean, m_sd;
 GaussianNG(const GaussianNG &n);
 GaussianNG(double min, double max, double mean, double sd);
 ~GaussianNG();
 double GetValue();
 std::string ToString();
private:
 double BoxMuller(double mean, double stddev);
};

class GeneralNG : public NG
{
public:
 double* m_n;
 double* m_p;
 double* m_cp;
 size_t  m_size;
 GeneralNG(const double n[], const double p[], size_t size);
 GeneralNG(const GeneralNG &n);
 ~GeneralNG();
 double GetValue();
 std::string ToString();
};

class RangesNG : public NG
{
public:
 double* m_n;
 double* m_p;
 double* m_cp;
 size_t  m_size;
 RangesNG(const double n[], const double p[], size_t size);
 RangesNG(const RangesNG &n);
 ~RangesNG();
 double GetValue();
 std::string ToString();
};

}

#endif /* _GEN_NUM_H_ */