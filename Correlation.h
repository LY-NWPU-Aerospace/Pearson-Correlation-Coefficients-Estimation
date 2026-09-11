#pragma once
#include <vector>
#include <set>
#include <string>

#define DO_NOT_EXIST_VALUE -801.403

#ifdef DLL_EXPORT
#define DLL_API __declspec(dllexport)
#else
#define DLL_API __declspec(dllimport)
#endif

class Correlation;

int eliminateLinearVariables(std::vector<std::vector<double>> in,
	std::vector<std::vector<double>>& out, std::vector<int>& outIndex,
	std::vector<int>& eliminatedIndex, double tolerance = 1e-3);

class DomainCenter
{
public:
	bool success;
	int num_iters;
	double final_obj;
	double time_cost;
	double condNum;
	double minEigValue;
	double determinant;
	std::vector<double> obj;
	std::vector<std::vector<double>> corr;
	std::vector<double> var;
	std::vector<double> eigValue;
	std::vector<std::vector<double>> eigVector;
};

struct IdentifiedResult
{
	std::vector<double> indenfied_d;
	std::vector<int> indenfied_index;

	std::vector<std::vector<int>> indenfied_factor_Matrix;
	std::vector<std::vector<double>> indenfied_d_Matrix;
};

struct Metrics
{
	double r_b;
	double eps_E_t;
	double d_r_b;
	double d_eps_E_t;
	std::vector<double> direction_eps;
	std::vector<double> direction_rb;
};

struct CorrResult
{
	DomainCenter chebyshev;
	DomainCenter analytic;
	Metrics paras;
	IdentifiedResult identRes;
};

struct MarginalDisInfo
{
	int Type = 0;				// 0-use PDF, 1-use CDF
	std::vector<double> value;
	std::vector<double> pdf;
	std::vector<double> cdf;
};

struct OptOption
{
	bool active = false;				// use opt or not
	int numThreads = 1;					// number of threads for opt DE algorithm
	int popSize = 1e3;					// population size
	int maxIterations = 1e3;			// maximum number of iterations
	std::string saveDir = "\\optRes";	// saved path
};

class DLL_API Correlation_API
{
public:
	Correlation_API();
	~Correlation_API();
	// set the corr matrix
	bool setCorrlationMatrix(std::vector<std::vector<double>> corr, Correlation* curCoef = nullptr);

	void setTolerance(double _value);

	CorrResult runEstimation(std::set<int> variableIndex, OptOption _option);

protected:

	Correlation* core;

	int curIter;
	double tolerance;
};

class DLL_API GaussianCopula
{
public:
	// correlation coefficient matrix
	std::vector<std::vector<double>> corr;

	// marginal distribution for variables
	std::vector<MarginalDisInfo> marDis;

	// sampling according to marDis and corr
	bool sampling(long _size, std::vector<std::vector<double>>& _samples);

};


#undef DO_NOT_EXIST_VALUE