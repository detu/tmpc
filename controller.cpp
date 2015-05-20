// S-function implementation of an MPC motion cueing controller.

#define S_FUNCTION_LEVEL 2
#define S_FUNCTION_NAME controller

#define NUM_INPUTS          2
/* Input Port  0 */
#define IN_PORT_0_NAME      y_ref
#define INPUT_DIMS_0_COL    1
#define INPUT_0_DTYPE       real_T
#define INPUT_0_COMPLEX     COMPLEX_NO
#define IN_0_FRAME_BASED    FRAME_NO
#define IN_0_BUS_BASED      0
#define IN_0_BUS_NAME       
#define IN_0_DIMS           1-D
#define INPUT_0_FEEDTHROUGH 1
#define IN_0_ISSIGNED        0
#define IN_0_WORDLENGTH      8
#define IN_0_FIXPOINTSCALING 1
#define IN_0_FRACTIONLENGTH  9
#define IN_0_BIAS            0
#define IN_0_SLOPE           0.125
/* Input Port  1 */
#define IN_PORT_1_NAME      x
#define INPUT_DIMS_1_COL    1
#define INPUT_1_DTYPE       real_T
#define INPUT_1_COMPLEX     COMPLEX_NO
#define IN_1_FRAME_BASED    FRAME_NO
#define IN_1_BUS_BASED      0
#define IN_1_BUS_NAME       
#define IN_1_DIMS           1-D
#define INPUT_1_FEEDTHROUGH 1
#define IN_1_ISSIGNED        0
#define IN_1_WORDLENGTH      8
#define IN_1_FIXPOINTSCALING 1
#define IN_1_FRACTIONLENGTH  9
#define IN_1_BIAS            0
#define IN_1_SLOPE           0.125

#define NUM_OUTPUTS          1
/* Output Port  0 */
#define OUT_PORT_0_NAME      u
#define OUTPUT_DIMS_0_COL    1
#define OUTPUT_0_DTYPE       real_T
#define OUTPUT_0_COMPLEX     COMPLEX_NO
#define OUT_0_FRAME_BASED    FRAME_NO
#define OUT_0_BUS_BASED      0
#define OUT_0_BUS_NAME       
#define OUT_0_DIMS           1-D
#define OUT_0_ISSIGNED        1
#define OUT_0_WORDLENGTH      8
#define OUT_0_FIXPOINTSCALING 1
#define OUT_0_FRACTIONLENGTH  3
#define OUT_0_BIAS            0
#define OUT_0_SLOPE           0.125

#define NPARAMS              0

#define SAMPLE_TIME_0        0.05
#define NUM_DISC_STATES      0
#define DISC_STATES_IC       [0]
#define NUM_CONT_STATES      0
#define CONT_STATES_IC       [0]

#define SFUNWIZ_GENERATE_TLC 0
#define SOURCEFILES "__SFB__"
#define PANELINDEX           6
#define USE_SIMSTRUCT        1
#define SHOW_COMPILE_STEPS   1                   
#define CREATE_DEBUG_MEXFILE 1
#define SAVE_CODE_ONLY       1
#define SFUNWIZ_REVISION     3.0

#include "simstruc.h"

#include "MotionPlatformX.hpp"

#include <qpDUNES.h>

/*
extern "C"
{
    return_t qpDUNES_setup(	qpData_t* const qpData,
						uint_t nI,
						uint_t nX,
						uint_t nU,
						uint_t* nD,
						qpOptions_t* options
						);
    
    return_t qpDUNES_cleanup(	qpData_t* const qpData
						);
}
 */

#include <memory>
#include <fstream>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrix;

// Motion platform to be used.
MotionPlatformX platform;

// Number of control intervals
const unsigned Nt = 1;

// Levenberg-Marquardt parameter
const double levenbergMarquardt = 0.01;

// Log stream.
std::ofstream log_stream;

extern void S_Outputs_wrapper(const real_T *y_ref,
			const real_T *x,
			real_T *u,
			const real_T *xD,
			SimStruct *S);
extern void S_Update_wrapper(const real_T *y_ref,
			const real_T *x,
			const real_T *u,
			real_T *xD,
			SimStruct *S);

Eigen::MatrixXd OutputWeightingMatrix()
{
	Eigen::MatrixXd w(6, 6);
	w.fill(0.);
	w.diagonal()[0] = 1;
	w.diagonal()[1] = 1;
	w.diagonal()[2] = 1;
	w.diagonal()[3] = 10;
	w.diagonal()[4] = 10;
	w.diagonal()[5] = 10;

	return w.transpose() * w;
}

/*====================*
 * S-function methods *
 *====================*/
/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *   Setup sizes of the various vectors.
 */
static void mdlInitializeSizes(SimStruct *S)
{
    DECL_AND_INIT_DIMSINFO(inputDimsInfo);
    DECL_AND_INIT_DIMSINFO(outputDimsInfo);
    ssSetNumSFcnParams(S, NPARAMS);
     if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
	 return; /* Parameter mismatch will be reported by Simulink */
     }

    ssSetNumContStates(S, NUM_CONT_STATES);
    ssSetNumDiscStates(S, NUM_DISC_STATES);

    if (!ssSetNumInputPorts(S, NUM_INPUTS)) return;
    /*Input Port 0 */
    ssSetInputPortWidth(S,  0, platform.getOutputDim()); /* */
    ssSetInputPortDataType(S, 0, SS_DOUBLE);
    ssSetInputPortComplexSignal(S,  0, INPUT_0_COMPLEX);
    ssSetInputPortDirectFeedThrough(S, 0, INPUT_0_FEEDTHROUGH);
    ssSetInputPortRequiredContiguous(S, 0, 1); /*direct input signal access*/

    /*Input Port 1 */
    ssSetInputPortWidth(S,  1, platform.getStateDim());
    ssSetInputPortDataType(S, 1, SS_DOUBLE);
    ssSetInputPortComplexSignal(S, 1, INPUT_1_COMPLEX);
    ssSetInputPortDirectFeedThrough(S, 1, INPUT_1_FEEDTHROUGH);
    ssSetInputPortRequiredContiguous(S, 1, 1); /*direct input signal access*/


    if (!ssSetNumOutputPorts(S, NUM_OUTPUTS)) return;
    ssSetOutputPortWidth(S, 0, platform.getInputDim());
    ssSetOutputPortDataType(S, 0, SS_DOUBLE);
    ssSetOutputPortComplexSignal(S, 0, OUTPUT_0_COMPLEX);
    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    
    // One PWork vector for storing a pointer to qpData_t.
    ssSetNumPWork(S, 1);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);

    ssSetSimulinkVersionGeneratedIn(S, "8.4");

    /* Take care when specifying exception free code - see sfuntmpl_doc.c */
    ssSetOptions(S, (/*SS_OPTION_EXCEPTION_FREE_CODE |*/ SS_OPTION_WORKS_WITH_CODE_REUSE));
}

# define MDL_SET_INPUT_PORT_FRAME_DATA
static void mdlSetInputPortFrameData(SimStruct  *S, 
                                     int_T      port,
                                     Frame_T    frameData)
{
    ssSetInputPortFrameData(S, port, frameData);
}
/* Function: mdlInitializeSampleTimes =========================================
 * Abstract:
 *    Specifiy  the sample time.
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, SAMPLE_TIME_0);
    ssSetOffsetTime(S, 0, 0.0);
}
#define MDL_INITIALIZE_CONDITIONS
 /* Function: mdlInitializeConditions ========================================
  * Abstract:
  *    Initialize the states
  */
 static void mdlInitializeConditions(SimStruct *S)
 {
	 /*
	 Initialize at default position, 0 velocity, 0 input.
	 */

	 log_stream.open("controller.log");
	 log_stream << "mdlInitializeConditions()" << std::endl;

	 using namespace Eigen;

	 // Get platform properties.
	 const auto Nq = platform.getNumberOfAxes();
	 const auto Nu = platform.getInputDim();
	 const auto Nx = platform.getStateDim();
	 const auto Ny = platform.getOutputDim();
	 const auto Nz = Nx + Nu;

	 VectorXd q_min(Nq), q_max(Nq), v_min(Nq), v_max(Nq), u_min(Nq), u_max(Nq);
	 platform.getAxesLimits(q_min.data(), q_max.data(), v_min.data(), v_max.data(), u_min.data(), u_max.data());

	 // State vector: x = [default_pos; 0]
	 Eigen::VectorXd x(Nx);
	 x.fill(0.);
	 platform.getDefaultAxesPosition(x.data());

	 // Input vector: u = 0
	 VectorXd u(Nu);
	 u.fill(0.);

	 // Output vector and derivatives at default initial state
	 VectorXd y(Ny);
	 RowMajorMatrix C(Ny, Nx);
	 RowMajorMatrix D(Ny, Nu);
	 platform.Output(x.data(), u.data(), y.data(), C.data(), D.data());

	 // Output weighting matrix
	 const auto W = OutputWeightingMatrix();

	 // G = [C, D]
	 RowMajorMatrix G(Ny, Nx + Nu);
	 G << C, D;

	 // H = G^T W G + \mu I
	 // Adding Levenberg-Marquardt term to make H positive-definite.
	 RowMajorMatrix H = G.transpose() * W * G + levenbergMarquardt * MatrixXd::Identity(Nz, Nz);		

	 // z = [x, u]
	 VectorXd z(Nx + Nu);
	 z << x, u;

	 // g = 2 * (y_bar - y_hat) * W * G
	 //   = 0 (s.t. y_bar = y_hat)
	 VectorXd g(Nz);
	 g.fill(0.);

	 // K = [A, B];
	 // x_{k+1} = K * z_k + c_k
	 MatrixXd A(Nx, Nx);
	 A << MatrixXd::Identity(Nq, Nq), SAMPLE_TIME_0 * MatrixXd::Identity(Nq, Nq),
		  MatrixXd::Zero(Nq, Nq), MatrixXd::Identity(Nq, Nq);

	 MatrixXd B(Nx, Nu);
	 B << SAMPLE_TIME_0 * SAMPLE_TIME_0 / 2. * MatrixXd::Identity(Nq, Nq),
		  SAMPLE_TIME_0 * MatrixXd::Identity(Nq, Nq);

	 RowMajorMatrix K(Nx, Nx + Nu);
	 K << A, B;

	 // c = 0
	 VectorXd c(Nx);
	 c.fill(0.);

	 // Initializing arrays.
	 // H_ stores Nt row-major matrices of size Nz x Nz and 1 matrix of size Nx x Nx.
	 std::vector<double> H_(Nz * Nz * Nt + Nx * Nx);

	 // g_ stores Nt vectors of size Nz and 1 vector of size Nx
	 std::vector<double> g_(Nz * Nt + Nx);

	 // C_ stores Nt row-major matrices of size Nx x Nz
	 std::vector<double> C_(Nx * Nz * Nt);

	 // c_ stores Nt vectors of size Nx
	 std::vector<double> c_(Nx * Nt);

	 // z_min stores Nt vectors of size Nz and 1 vector of size Nx
	 std::vector<double> z_min(Nz * Nt + Nx);

	 // z_max stores Nt vectors of size Nz and 1 vector of size Nx
	 std::vector<double> z_max(Nz * Nt + Nx);

	 for (unsigned i = 0; i < Nt; ++i)
	 {
		 Map<RowMajorMatrix>(H_.data() + Nz * Nz * i, Nz, Nz) = H;
		 Map<VectorXd>(g_.data() + Nz * i, Nz) = g;
		 Map<RowMajorMatrix>(C_.data() + Nx * Nz * i, Nx, Nz) = K;
		 Map<VectorXd>(c_.data() + Nx * i, Nx, 1) = c;
		 
		 Map<VectorXd>(z_min.data() + Nz * i, Nz, 1) << q_min, v_min, u_min;
		 Map<VectorXd>(z_max.data() + Nz * i, Nz, 1) << q_max, v_max, u_max;
	 }

	 Map<RowMajorMatrix>(H_.data() + Nz * Nz * Nt, Nx, Nx) = levenbergMarquardt * MatrixXd::Identity(Nx, Nx);
 	 Map<VectorXd>(g_.data() + Nz * Nt, Nx).fill(0.);
 	 Map<VectorXd>(z_min.data() + Nz * Nt, Nx) << q_min, v_min;
 	 Map<VectorXd>(z_max.data() + Nz * Nt, Nx) << q_max, v_max;
	 
     auto qpData = reinterpret_cast<qpData_t *>(ssGetPWorkValue(S, 0));
     
     /** Initial MPC data setup: components not given here are set to zero (if applicable)
	 *      instead of passing g, D, zLow, zUpp, one can also just pass NULL pointers (0) */

	 log_stream << "H = " << std::endl;
	 for (unsigned i = 0; i < H_.size(); ++i)
	 {
		 log_stream << H_[i] << '\t';

		 if (i < Nz * Nz * Nt)
		 {
			 if ((i + 1) % Nz == 0)
				 log_stream << std::endl;

			 if ((i + 1) % (Nz * Nz) == 0)
				 log_stream << std::endl;
		 }
		 else
		 {
			 if ((i - Nz * Nz * Nt + 1) % Nx == 0)
				 log_stream << std::endl;

			 if ((i - Nz * Nz * Nt + 1) % (Nx * Nx) == 0)
				 log_stream << std::endl;
		 }
	 }

	 log_stream << "g = " << std::endl;
	 for (unsigned i = 0; i < g_.size(); ++i)
	 {
		 log_stream << g_[i] << '\t';

		 if (i < Nz * Nt)
		 {
			 if ((i + 1) % Nz == 0)
				 log_stream << std::endl;
		 }
		 else
		 {
			 if ((i - Nz * Nt + 1) % Nx == 0)
				 log_stream << std::endl;
		 }
	 }
	 log_stream << std::endl;

	 log_stream << "C = " << std::endl;
	 for (unsigned i = 0; i < C_.size(); ++i)
	 {
		 log_stream << C_[i] << '\t';

  		 if ((i + 1) % Nz == 0)
			log_stream << std::endl;

		 if ((i + 1) % (Nz * Nx) == 0)
			log_stream << std::endl;
	 }

	 log_stream << "c = " << std::endl;
	 for (unsigned i = 0; i < c_.size(); ++i)
	 {
		 log_stream << c_[i] << '\t';

		 if ((i + 1) % Nx == 0)
			 log_stream << std::endl;
	 }
	 log_stream << std::endl;

	 log_stream << "z_min = " << std::endl;
	 for (unsigned i = 0; i < z_min.size(); ++i)
	 {
		 log_stream << z_min[i] << '\t';

		 if ((i + 1) % Nz == 0)
			 log_stream << std::endl;
	 }
	 log_stream << std::endl << std::endl;

	 log_stream << "z_max = " << std::endl;
	 for (unsigned i = 0; i < z_max.size(); ++i)
	 {
		 log_stream << z_max[i] << '\t';

		 if ((i + 1) % Nz == 0)
			 log_stream << std::endl;
	 }
	 log_stream << std::endl << std::endl;
	 
	 return_t statusFlag = qpDUNES_init(qpData, H_.data(), g_.data(), C_.data(), c_.data(), z_min.data(), z_max.data(), 0, 0, 0);
	 if (statusFlag != QPDUNES_OK) 
     {
  		 ssSetErrorStatus(S, "qpDUNES_init() failed.");
		 return;
	 }
 }

#define MDL_START  /* Change to #undef to remove function */
#if defined(MDL_START) 
/* Function: mdlStart =======================================================
* Abstract:
*    This function is called once at start of model execution. If you
*    have states that should be initialized once, this is the place
*    to do it.
*/
static void mdlStart(SimStruct *S)
{
    /** Set qpDUNES options */
	qpOptions_t qpOptions 			= qpDUNES_setupDefaultOptions();
	qpOptions.maxIter    			= 100;
	qpOptions.printLevel 			= 2;
	qpOptions.stationarityTolerance = 1.e-6;
    
    /** Allocate data for qpDUNES and set options */
    auto qpData = std::make_unique<qpData_t>();
    
    unsigned int* nD = 0;	  			/* number of affine constraints */
	return_t statusFlag = qpDUNES_setup(qpData.get(), Nt, platform.getStateDim(), platform.getInputDim(), nD, &qpOptions);
	if (statusFlag != QPDUNES_OK)
	{
		ssSetErrorStatus(S, "qpDUNES_setup() failed.");
		return;
	}
    
    ssSetPWorkValue(S, 0, qpData.release());
}
#endif /*  MDL_START */

#define MDL_SET_INPUT_PORT_DATA_TYPE
static void mdlSetInputPortDataType(SimStruct *S, int port, DTypeId dType)
{
    ssSetInputPortDataType( S, 0, dType);
}
#define MDL_SET_OUTPUT_PORT_DATA_TYPE
static void mdlSetOutputPortDataType(SimStruct *S, int port, DTypeId dType)
{
    ssSetOutputPortDataType(S, 0, dType);
}

#define MDL_SET_DEFAULT_PORT_DATA_TYPES
static void mdlSetDefaultPortDataTypes(SimStruct *S)
{
  ssSetInputPortDataType( S, 0, SS_DOUBLE);
  ssSetOutputPortDataType(S, 0, SS_DOUBLE);
}

/* Function: mdlOutputs =======================================================
 *
*/
static void mdlOutputs(SimStruct *S, int_T tid)
{
	using namespace Eigen;

	log_stream << "mdlOutputs()" << std::endl;

    const real_T   *y_ref  = (const real_T*) ssGetInputPortSignal(S,0);
    const real_T   *x  = (const real_T*) ssGetInputPortSignal(S,1);
    real_T         *u  = (real_T *)ssGetOutputPortRealSignal(S,0);
	
	auto qpData = reinterpret_cast<qpData_t *>(ssGetPWorkValue(S, 0));

	// Get platform properties.
	const auto Nq = platform.getNumberOfAxes();
	const auto Nu = platform.getInputDim();
	const auto Nx = platform.getStateDim();
	const auto Ny = platform.getOutputDim();
	const auto Nz = Nx + Nu;

	VectorXd q_min(Nq), q_max(Nq), v_min(Nq), v_max(Nq), u_min(Nq), u_max(Nq);
	platform.getAxesLimits(q_min.data(), q_max.data(), v_min.data(), v_max.data(), u_min.data(), u_max.data());

	Map<const VectorXd> map_x(x, Nx);
	Map<const VectorXd> map_y_ref(y_ref, Ny);

	/** (1) embed current initial value */
	VectorXd z0_min(Nz), z0_max(Nz);
	z0_min << map_x, u_min;
	z0_max << map_x, u_max;

	log_stream << "z0_min = " << std::endl;
	for (unsigned i = 0; i < z0_min.size(); ++i)
		log_stream << z0_min[i] << '\t';
	log_stream << std::endl << std::endl;

	log_stream << "z0_max = " << std::endl;
	for (unsigned i = 0; i < z0_max.size(); ++i)
		log_stream << z0_max[i] << '\t';
	log_stream << std::endl << std::endl;

	return_t statusFlag = qpDUNES_updateIntervalData(qpData, qpData->intervals[0], 0, 0, 0, 0, z0_min.data(), z0_max.data(), 0, 0, 0, 0);
	if (statusFlag != QPDUNES_OK) 
	{
		ssSetErrorStatus(S, "Initial value embedding failed (qpDUNES_updateIntervalData()).");
		return;
	}

	/** (2) solve QP */
	statusFlag = qpDUNES_solve(qpData);
	if (statusFlag != QPDUNES_SUCC_OPTIMAL_SOLUTION_FOUND) 
	{
		ssSetErrorStatus(S, "QP solution failed (qpDUNES_solve()).");
		return;
	}

	/** (3) obtain primal and dual optimal solution */
	// z_opt contains Nt vectors of size Nz and 1 vector of size Nx.
	std::vector<double> z_opt(Nz * Nt + Nx);
	qpDUNES_getPrimalSol(qpData, z_opt.data());
	//qpDUNES_getDualSol(&qpData, lambdaOpt, muOpt);

	// Send optimal input for time 0 to output.
	std::copy_n(z_opt.begin() + Nx, Nu, u);
}

#define MDL_UPDATE  /* Change to #undef to remove function */
/* Function: mdlUpdate ======================================================
   * Abstract:
   *    This function is called once for every major integration time step.
   *    Discrete states are typically updated here, but this function is useful
   *    for performing any tasks that should only take place once per
   *    integration step.
   */
  static void mdlUpdate(SimStruct *S, int_T tid)
  {
    real_T         *xD  = ssGetDiscStates(S);
    const real_T   *y_ref  = (const real_T*) ssGetInputPortSignal(S,0);
    const real_T   *x  = (const real_T*) ssGetInputPortSignal(S,1);
    real_T        *u  = (real_T *)ssGetOutputPortRealSignal(S,0);

    S_Update_wrapper(y_ref, x, u,  xD, S);
}


/* Function: mdlTerminate =====================================================
 * Abstract:
 *    In this function, you should perform any actions that are necessary
 *    at the termination of a simulation.  For example, if memory was
 *    allocated in mdlStart, this is the place to free it.
 */
static void mdlTerminate(SimStruct *S)
{
    auto qpData = reinterpret_cast<qpData_t *>(ssGetPWorkValue(S, 0));
     
    /** cleanup of allocated data */
	qpDUNES_cleanup(qpData);
    
    delete qpData;
}

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif


