#pragma once
#include <vector>
#include <set>
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

struct DikinEllipsoid
{
public:
	std::vector<std::vector<double>> H;			// Hessian matrix
	std::vector<std::vector<double>> H_i;		// Inverse of H
	std::vector<std::vector<double>> eigVec;	// Eigen vectors, row vector
	std::vector<double> x_center;				// The center of the ellipsoid
	std::vector<double> x_lower;				// The lower bounds of the ellipsoid
	std::vector<double> x_upper;				// The upper bounds of the ellipsoid
	std::vector<double> length;					// The axial length of the ellipsoid
};

struct DikinCenter
{
	bool success;
	int num_iters;
	int opt_type;
	double final_obj;
	double time_cost;
	std::vector<double> obj;
	std::vector<std::vector<double>> corr;
	std::vector<double> var;

	std::vector<double> eigValue;
	double condNum;
	double minEigValue;
	double determinant;

};

struct ErrorBenefit
{
	double maxErrorR_Euclidean;
	double ZaBoundsR_Euclidean;
	double maxError_Riemannian;
	double ZaBounds_Riemannian;
	double minBenefit;
	double maxBeneift;
	double r_b;
};

struct ErrorImprove
{
	double supError;
	double improvement;
	std::vector<double> MaxCorrelated;
	std::vector<double> MinCorrelated;
};

struct CorrResult
{
	DikinCenter chebyshev;
	DikinCenter analytic;
	DikinEllipsoid dikin_base;
	DikinEllipsoid dikin_inflation;
	ErrorBenefit paras;

	std::vector<std::vector<double>> corr;
	std::vector<std::vector<double>> longest;
};

struct MarginalDisInfo
{
	int Type = 0;				// 0-use PDF, 1-use CDF
	std::vector<double> value;
	std::vector<double> pdf;
	std::vector<double> cdf;
};

class DLL_API Correlation_API
{
public:
	Correlation_API();
	~Correlation_API();
	// set the corr matrix
	bool setCorrlationMatrix(std::vector<std::vector<double>> corr, Correlation* curCoef = nullptr);

	// set the tolerance of the convergence
	void setTolerance(double _value);

	// estimate the bounds of the unknown elements of the matrix
	CorrResult runEstimation();

	std::vector<std::vector<double>> getBoundryMatrix(std::vector<double> direction);

	double getBoundryLength(std::vector<double> direction);

	double getBoundryLength(std::set<int> index);

	ErrorImprove getImproveIndex(std::set<int> variableIndex);

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