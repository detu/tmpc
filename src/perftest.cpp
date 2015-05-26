#include <CyberMotion.hpp>

#include <Eigen/Dense>

#include <vector>
#include <ctime>
#include <iostream>

int main(int argc, char * argv[])
{
	CyberMotion cms;

// 	Eigen::MatrixXd x(cms.getStateDim(), nT);
// 	x.fill(0.);
// 
// 	for (size_t i = 0; i < nT; ++i)
// 		cms.getDefaultAxesPosition(x.col(i).data());
// 
// 	Eigen::MatrixXd u(cms.getInputDim(), nT);
// 	u.fill(0.);
// 
// 	Eigen::VectorXd y(cms.getOutputDim());
// 	Eigen::MatrixXd C(cms.getOutputDim(), cms.getStateDim());
// 	Eigen::MatrixXd D(cms.getOutputDim(), cms.getInputDim());
	const unsigned nX = cms.getStateDim(), nU = cms.getInputDim(), nY = cms.getOutputDim();

	{
		const size_t n_iter = 100000;

		std::vector<double> x(nX, 0);
		std::vector<double> u(nU, 0);
		std::vector<double> y(nY);
		std::vector<double> C(nY * nX);
		std::vector<double> D(nY * nU);

		double * px = x.data();
		double * pu = u.data();
		double * py = y.data();
		double * pC = C.data();
		double * pD = D.data();

		cms.getDefaultAxesPosition(px);

		std::cout << "Running single-vector test..." << std::endl;
		const clock_t t0 = clock();
		for (size_t i = 0; i < n_iter; ++i)
			cms.Output(px, pu, py, pC, pD);
		const clock_t t1 = clock();
		const double elapsed_time = (double)(t1 - t0) / CLOCKS_PER_SEC;

		std::cout << "Time elapsed: " << elapsed_time << "s, time per step: " << elapsed_time / n_iter << std::endl;
	}

	{
		const size_t n_iter = 100;
		const unsigned nT = 1000;

		std::vector<double> x(nX * nT, 0);
		std::vector<double> u(nU * nT, 0);
		std::vector<double> y(nY * nT);
		std::vector<double> C(nY * nX * nT);
		std::vector<double> D(nY * nU * nT);

		double * px = x.data();
		double * pu = u.data();
		double * py = y.data();
		double * pC = C.data();
		double * pD = D.data();

		for (size_t j = 0; j < nT; ++j)
			cms.getDefaultAxesPosition(px + j * nX);

		std::cout << "Running multiple-vector test..." << std::endl;
		const clock_t t0 = clock();
		for (size_t i = 0; i < n_iter; ++i)
			for (size_t j = 0; j < nT; ++j)
				cms.Output(px + j * nX, pu + j * nU, py + j * nY, pC + j * nY * nX, pD + j * nY * nU);
		const clock_t t1 = clock();
		const double elapsed_time = (double)(t1 - t0) / CLOCKS_PER_SEC;

		std::cout << "Time elapsed: " << elapsed_time << "s, time per step: " << elapsed_time / (n_iter * nT) << std::endl;
	}
}