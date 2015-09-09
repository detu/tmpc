#include <CyberMotion.hpp>

#include <Eigen/Dense>

#include <vector>
#include <ctime>
#include <iostream>

int main(int argc, char * argv[])
{
	mpmc::CyberMotion cms;

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

		Eigen::VectorXd x = Eigen::VectorXd::Zero(nX);
		Eigen::VectorXd u = Eigen::VectorXd::Zero(nU);
		Eigen::VectorXd y = Eigen::VectorXd::Zero(nY);
		Eigen::MatrixXd C = Eigen::MatrixXd::Zero(nY, nX);
		Eigen::MatrixXd D = Eigen::MatrixXd::Zero(nY, nU);

		x.topRows(cms.getNumberOfAxes()) = cms.getDefaultAxesPosition();

		std::cout << "Running single-vector test..." << std::endl;
		const clock_t t0 = clock();
		for (size_t i = 0; i < n_iter; ++i)
			cms.Output(x, u, y, C, D);
		const clock_t t1 = clock();
		const double elapsed_time = (double)(t1 - t0) / CLOCKS_PER_SEC;

		std::cout << "Time elapsed: " << elapsed_time << "s, time per step: " << elapsed_time / n_iter << std::endl;
	}

	{
		const size_t n_iter = 100;
		const unsigned nT = 1000;

		Eigen::MatrixXd x = Eigen::MatrixXd::Zero(nX, nT);
		Eigen::MatrixXd u = Eigen::MatrixXd::Zero(nU, nT);
		Eigen::MatrixXd y = Eigen::MatrixXd::Zero(nY, nT);
		Eigen::MatrixXd C = Eigen::MatrixXd::Zero(nY, nX * nT);
		Eigen::MatrixXd D = Eigen::MatrixXd::Zero(nY, nU * nT);

		for (size_t j = 0; j < nT; ++j)
			x.col(j).topRows(cms.getNumberOfAxes()) = cms.getDefaultAxesPosition();

		std::cout << "Running multiple-vector test..." << std::endl;
		const clock_t t0 = clock();
		for (size_t i = 0; i < n_iter; ++i)
			for (size_t j = 0; j < nT; ++j)
			{
				Eigen::VectorXd yj(nY);
				Eigen::MatrixXd Cj(nY, nX);
				Eigen::MatrixXd Dj(nY, nU);
				cms.Output(x.col(j), u.col(j), yj, Cj, Dj);

				y.col(j) = yj;
				C.middleCols(nX * j, nX) = Cj;
				D.middleCols(nU * j, nU) = Dj;
			}

		const clock_t t1 = clock();
		const double elapsed_time = (double)(t1 - t0) / CLOCKS_PER_SEC;

		std::cout << "Time elapsed: " << elapsed_time << "s, time per step: " << elapsed_time / (n_iter * nT) << std::endl;
	}
}