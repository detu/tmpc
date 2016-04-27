#include <casadi_interface/GeneratedFunction.hpp>

#include <../src/cybermotion_generated.h>

#include <gtest/gtest.h>
//#define EXPECT_TRUE(X) assert(X)

#include <Eigen/Dense>

#include <iostream>
#include <array>

TEST(test_cybermotion_generated, ode_test)
{
	casadi_interface::GeneratedFunction ode {CASADI_GENERATED_FUNCTION_INTERFACE(cybermotion_ode)};

	EXPECT_EQ(ode.n_in(), 2);
	EXPECT_EQ(ode.n_out(), 2);
	EXPECT_EQ(ode.name(), "cybermotion_ode");

	const unsigned NX = 16;
	const unsigned NU = 8;
	unsigned const NZ = NX + NU;
	const unsigned NP = 0;

	Eigen::Matrix<double, NX, 1> x0;
	x0 << 4.8078, 0.1218, -1.5319, 0.4760, 0.0006, 0.1396, -0.0005, 0.7991, 0, 0, 0, 0, 0, 0, 0, 0;

	Eigen::Matrix<double, NU, 1> u0;
	u0 << 0, 0, 0, 0, 0, 0, 0, 0;

	Eigen::Matrix<double, NX + NU, 1> z0;
	z0 << x0, u0;

	Eigen::Matrix<double, NP, 1> p;

	Eigen::Matrix<double, NX, 1> xdot, xdot_expected;
	Eigen::Matrix<double, NX, NZ> ode_jac, ode_jac_expected;
	Eigen::Matrix<double, NX, NX> A_expected;
	Eigen::Matrix<double, NX, NU> B_expected;

	ode({z0.data(), p.data()}, {xdot.data(), ode_jac.data()});

	/*
	xdot_expected <<
			3.508514747852236049e+00,
			-1.585945328501758400e+00,
			-1.764017709158380853e+01,
			-6.556499999999999828e-02,
			4.682499999999999857e-02,
			9.365500000000001601e-02,
			1.404950000000000088e-01,
			-1.448502944637109913e+01,
			-1.799685562253899107e+01,
			-2.062522132807076858e+01,
			1.000000000000000000e+00,
			2.000000000000000000e+00,
			3.000000000000000000e+00;

	A_expected <<
			-2.263100810676232921e-02,-5.534941559225998314e-02,-2.049643162203821159e-02,4.178874199675344236e+00,-2.061383816039205996e+01,3.410203609512969791e+01,1.841307647932706004e+01,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,-7.755979641141214076e-01,-3.865568766282102553e-01,7.172153584796815418e-01,
			5.062833560940012412e-02,2.229875907023290763e-02,1.895601843988546642e-02,3.218942756800190708e+00,5.257577158912633664e+00,-3.090790925509604392e+01,-1.797459961094623893e+01,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,6.800903710270614688e+00,1.191172357751306254e+00,6.022681908433789033e+00,
			6.535237970533956317e-03,-2.216721991173725356e-02,3.322490365294036439e-04,-5.386078804066578130e+01,-2.722779507772512986e+00,-1.269941895292862633e+01,-7.166891079292938116e+01,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,3.273503843206973518e-01,6.215973597752966207e-01,1.613652790671884718e+00,
			-4.685000000000000275e-02,-9.364999999999999714e-02,-1.405000000000000138e-01,0.000000000000000000e+00,-5.000000000000000278e-02,-1.000000000000000056e-01,-1.499999999999999944e-01,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			4.682999999999999940e-01,-1.405000000000000138e-01,9.364999999999999714e-02,5.000000000000000278e-02,0.000000000000000000e+00,1.499999999999999944e-01,-1.000000000000000056e-01,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			1.405000000000000138e-01,4.682999999999999940e-01,-4.685000000000000275e-02,1.000000000000000056e-01,-1.499999999999999944e-01,0.000000000000000000e+00,5.000000000000000278e-02,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			-9.364999999999999714e-02,4.685000000000000275e-02,4.682999999999999940e-01,1.499999999999999944e-01,1.000000000000000056e-01,-5.000000000000000278e-02,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,-9.328615448764224993e+00,-5.745398244127795095e+00,8.495427402413231022e+00,1.614349158821019592e-01,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,-5.019901754601805877e+00,-1.321900087885137909e+00,-1.132263292073623751e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,-1.178500424491876508e+01,-2.186257078612726890e+00,-6.465069340539971732e-01,-1.857584053029406768e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,-1.321900087885137465e+00,-4.121667227798377908e+00,1.775680651130504994e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,5.760571408587503939e+00,-8.687463367666453662e-01,1.082087917065181593e+01,6.560521865859402446e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,-1.132263292073623973e+00,1.775680651130505217e+00,-5.213537935891887187e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,1.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,1.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,1.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00;

	B_expected <<
			2.841837172418348531e-03,-1.712441215785584156e-02,9.100620991552598946e-03,-9.408152338971555173e-03,-7.744037888739661347e-03,1.077150679659897789e-02,1.262076153907665038e-03,5.272113269609198828e-03,
			-2.155764523418365625e-03,-1.293448722788948652e-02,-5.877652777871859077e-04,1.271208738279471645e-02,4.586030814371005147e-03,5.112958220081481597e-03,6.338400379705978575e-03,-1.687156512605725300e-02,
			-6.250197705115737370e-03,2.232493556741287077e-03,1.543069811056716444e-02,-6.197248195433956164e-03,-1.504035832208307569e-02,-1.441161900325400175e-02,-9.022282803759381492e-03,3.590778076182021827e-03,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			-9.678357977215825900e-03,-7.041247046053312683e-03,-5.862770653072353795e-03,-2.515503211485183304e-04,-6.821459077441941517e-04,-1.013215971818371254e-03,-6.773932840708953712e-03,-5.698421736997943098e-03,
			2.418488824071584730e-03,1.003147260244785163e-03,6.362186589894427285e-03,1.692468737547406373e-04,-9.963191755352077475e-03,-8.548720175291079498e-03,-7.356115778925206329e-03,-6.444474569208749962e-03,
			-6.930358379648532624e-04,-7.029547326158982340e-03,5.015027623537191524e-03,-9.995402839888660809e-03,5.191213794790835653e-04,-5.088494547410207271e-03,3.734322712765715443e-05,-5.098699563038698937e-03,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00;
	*/

	/*
	std::cout << "delta(xdot) = " << std::endl << xdot - xdot_expected << std::endl;
	std::cout << "delta(A) = " << std::endl << A - A_expected << std::endl;
	std::cout << "delta(B) = " << std::endl << B - B_expected << std::endl;
	*/

	std::cout << "xdot = " << std::endl << xdot << std::endl;
	std::cout << "ode_jac = " << std::endl << ode_jac << std::endl;

	/*
	EXPECT_TRUE(xdot.isApprox(xdot_expected));
	EXPECT_TRUE(A.isApprox(A_expected));
	EXPECT_TRUE(B.isApprox(B_expected));
	*/
}

TEST(test_cybermotion_generated, output_test)
{
	casadi_interface::GeneratedFunction output {CASADI_GENERATED_FUNCTION_INTERFACE(cybermotion_output)};

	EXPECT_EQ(output.n_in(), 2);
	EXPECT_EQ(output.n_out(), 2);
	EXPECT_EQ(output.name(), "cybermotion_output");

	const unsigned NX = 16;
	const unsigned NU = 8;
	unsigned const NZ = NX + NU;
	const unsigned NP = 0;
	const unsigned NY = 9;

	Eigen::Matrix<double, NX, 1> x0;
	x0 << 4.8078, 0.1218, -1.5319, 0.4760, 0.0006, 0.1396, -0.0005, 0.7991, 0, 0, 0, 0, 0, 0, 0, 0;

	Eigen::Matrix<double, NU, 1> u0;
	u0 << 0, 0, 0, 0, 0, 0, 0, 0;

	Eigen::Matrix<double, NX + NU, 1> z0;
	z0 << x0, u0;

	Eigen::Matrix<double, NP, 1> p;

	Eigen::Matrix<double, NY, 1> y, y_expected;
	Eigen::Matrix<double, NY, NZ> jac, jac_expected;
	Eigen::Matrix<double, NY, NX> C_expected;
	Eigen::Matrix<double, NY, NU> D_expected;

	output({z0.data(), p.data()}, {y.data(), jac.data()});

	/*
	xdot_expected <<
			3.508514747852236049e+00,
			-1.585945328501758400e+00,
			-1.764017709158380853e+01,
			-6.556499999999999828e-02,
			4.682499999999999857e-02,
			9.365500000000001601e-02,
			1.404950000000000088e-01,
			-1.448502944637109913e+01,
			-1.799685562253899107e+01,
			-2.062522132807076858e+01,
			1.000000000000000000e+00,
			2.000000000000000000e+00,
			3.000000000000000000e+00;

	A_expected <<
			-2.263100810676232921e-02,-5.534941559225998314e-02,-2.049643162203821159e-02,4.178874199675344236e+00,-2.061383816039205996e+01,3.410203609512969791e+01,1.841307647932706004e+01,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,-7.755979641141214076e-01,-3.865568766282102553e-01,7.172153584796815418e-01,
			5.062833560940012412e-02,2.229875907023290763e-02,1.895601843988546642e-02,3.218942756800190708e+00,5.257577158912633664e+00,-3.090790925509604392e+01,-1.797459961094623893e+01,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,6.800903710270614688e+00,1.191172357751306254e+00,6.022681908433789033e+00,
			6.535237970533956317e-03,-2.216721991173725356e-02,3.322490365294036439e-04,-5.386078804066578130e+01,-2.722779507772512986e+00,-1.269941895292862633e+01,-7.166891079292938116e+01,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,3.273503843206973518e-01,6.215973597752966207e-01,1.613652790671884718e+00,
			-4.685000000000000275e-02,-9.364999999999999714e-02,-1.405000000000000138e-01,0.000000000000000000e+00,-5.000000000000000278e-02,-1.000000000000000056e-01,-1.499999999999999944e-01,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			4.682999999999999940e-01,-1.405000000000000138e-01,9.364999999999999714e-02,5.000000000000000278e-02,0.000000000000000000e+00,1.499999999999999944e-01,-1.000000000000000056e-01,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			1.405000000000000138e-01,4.682999999999999940e-01,-4.685000000000000275e-02,1.000000000000000056e-01,-1.499999999999999944e-01,0.000000000000000000e+00,5.000000000000000278e-02,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			-9.364999999999999714e-02,4.685000000000000275e-02,4.682999999999999940e-01,1.499999999999999944e-01,1.000000000000000056e-01,-5.000000000000000278e-02,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,-9.328615448764224993e+00,-5.745398244127795095e+00,8.495427402413231022e+00,1.614349158821019592e-01,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,-5.019901754601805877e+00,-1.321900087885137909e+00,-1.132263292073623751e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,-1.178500424491876508e+01,-2.186257078612726890e+00,-6.465069340539971732e-01,-1.857584053029406768e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,-1.321900087885137465e+00,-4.121667227798377908e+00,1.775680651130504994e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,5.760571408587503939e+00,-8.687463367666453662e-01,1.082087917065181593e+01,6.560521865859402446e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,-1.132263292073623973e+00,1.775680651130505217e+00,-5.213537935891887187e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,1.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,1.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,1.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00;

	B_expected <<
			2.841837172418348531e-03,-1.712441215785584156e-02,9.100620991552598946e-03,-9.408152338971555173e-03,-7.744037888739661347e-03,1.077150679659897789e-02,1.262076153907665038e-03,5.272113269609198828e-03,
			-2.155764523418365625e-03,-1.293448722788948652e-02,-5.877652777871859077e-04,1.271208738279471645e-02,4.586030814371005147e-03,5.112958220081481597e-03,6.338400379705978575e-03,-1.687156512605725300e-02,
			-6.250197705115737370e-03,2.232493556741287077e-03,1.543069811056716444e-02,-6.197248195433956164e-03,-1.504035832208307569e-02,-1.441161900325400175e-02,-9.022282803759381492e-03,3.590778076182021827e-03,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			-9.678357977215825900e-03,-7.041247046053312683e-03,-5.862770653072353795e-03,-2.515503211485183304e-04,-6.821459077441941517e-04,-1.013215971818371254e-03,-6.773932840708953712e-03,-5.698421736997943098e-03,
			2.418488824071584730e-03,1.003147260244785163e-03,6.362186589894427285e-03,1.692468737547406373e-04,-9.963191755352077475e-03,-8.548720175291079498e-03,-7.356115778925206329e-03,-6.444474569208749962e-03,
			-6.930358379648532624e-04,-7.029547326158982340e-03,5.015027623537191524e-03,-9.995402839888660809e-03,5.191213794790835653e-04,-5.088494547410207271e-03,3.734322712765715443e-05,-5.098699563038698937e-03,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00;
	*/

	/*
	std::cout << "delta(xdot) = " << std::endl << xdot - xdot_expected << std::endl;
	std::cout << "delta(A) = " << std::endl << A - A_expected << std::endl;
	std::cout << "delta(B) = " << std::endl << B - B_expected << std::endl;
	*/

	std::cout << "y = " << std::endl << y << std::endl;
	std::cout << "jac = " << std::endl << jac << std::endl;

	/*
	EXPECT_TRUE(xdot.isApprox(xdot_expected));
	EXPECT_TRUE(A.isApprox(A_expected));
	EXPECT_TRUE(B.isApprox(B_expected));
	*/
}
