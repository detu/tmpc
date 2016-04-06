// S-function implementation of an MPC motion cueing controller.

#define S_FUNCTION_LEVEL 2

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

#define NPARAMS              11

#define SAMPLE_TIME_0        0.012
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
#include "CyberMotion1D.hpp"
#include "MotionPlatformModelPredictiveController.hpp"

#include <memory>
#include <fstream>

// Motion platform to use
auto platform = std::make_shared<mpmc::CyberMotion1D>();
const unsigned n_axes = 8;

// Log stream.
std::ofstream log_stream;

// Error status for Simulink
std::string error_status;

mpmc::MotionPlatformModelPredictiveController * getController(SimStruct * S)
{
	return reinterpret_cast<mpmc::MotionPlatformModelPredictiveController *>(ssGetPWorkValue(S, 0));
}

void setWorkingPointInitialized(SimStruct * S, bool val)
{
	ssSetIWorkValue(S, 0, val ? 1 : 0);
}

bool getWorkingPointInitialized(SimStruct * S)
{
	return ssGetIWorkValue(S, 0) != 0;
}

void setController(SimStruct * S, mpmc::MotionPlatformModelPredictiveController * c)
{
	ssSetPWorkValue(S, 0, c);
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
	ssSetInputPortWidth(S,  0, platform->getOutputDim()); /* */
	ssSetInputPortDataType(S, 0, SS_DOUBLE);
	ssSetInputPortComplexSignal(S,  0, INPUT_0_COMPLEX);
	ssSetInputPortDirectFeedThrough(S, 0, INPUT_0_FEEDTHROUGH);
	ssSetInputPortRequiredContiguous(S, 0, 1); /*direct input signal access*/

	/*Input Port 1 */
	ssSetInputPortWidth(S, 1, 2 * n_axes);
	ssSetInputPortDataType(S, 1, SS_DOUBLE);
	ssSetInputPortComplexSignal(S, 1, INPUT_1_COMPLEX);
	ssSetInputPortDirectFeedThrough(S, 1, INPUT_1_FEEDTHROUGH);
	ssSetInputPortRequiredContiguous(S, 1, 1); /*direct input signal access*/


	if (!ssSetNumOutputPorts(S, NUM_OUTPUTS)) return;
	ssSetOutputPortWidth(S, 0, n_axes);
	ssSetOutputPortDataType(S, 0, SS_DOUBLE);
	ssSetOutputPortComplexSignal(S, 0, OUTPUT_0_COMPLEX);
	ssSetNumSampleTimes(S, 1);
	ssSetNumRWork(S, 0);
	ssSetNumIWork(S, 1);	// Reserve 1 IWork for the "Working Point Initialized" flag.
	
	// One PWork vector for storing a pointer to MPC_Controller.
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
 *    Specify  the sample time.
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
	 Reset the "Working Point Initialized" flag.
	 */
	setWorkingPointInitialized(S, false);
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
	using Eigen::Map;
	using Eigen::VectorXd;

	/** Initialize MPC_Controller */

	// 
	const mxArray * mx_n_intervals = ssGetSFcnParam(S, 8);
	if (!(mx_n_intervals && mxGetNumberOfElements(mx_n_intervals) == 1))
		throw std::runtime_error("nIntervals must be a scalar.");

	double n_intervals = 0;
	if (!((n_intervals = mxGetScalar(mx_n_intervals)) == std::floor(n_intervals) && n_intervals > 0))
		throw std::runtime_error("nIntervals must be a positive integer.");

	auto controller = std::make_unique<mpmc::MotionPlatformModelPredictiveController>(platform, SAMPLE_TIME_0, static_cast<unsigned>(n_intervals));
	controller->setLevenbergMarquardt(0.01);
	
	// Set motion limits from block parameters.
	const mxArray * mx_qMin = ssGetSFcnParam(S, 0);
	if(!mx_qMin || mxGetNumberOfElements(mx_qMin) != n_axes)
		throw std::runtime_error("Invalid number of elements for parameter qMin");
	
	const mxArray * mx_qMax = ssGetSFcnParam(S, 1);
	if(!mx_qMax || mxGetNumberOfElements(mx_qMax) != n_axes)
		throw std::runtime_error("Invalid number of elements for parameter qMax");
	
	const mxArray * mx_vMin = ssGetSFcnParam(S, 2);
	if(!mx_vMin || mxGetNumberOfElements(mx_vMin) != n_axes)
		throw std::runtime_error("Invalid number of elements for parameter vMin");
	
	const mxArray * mx_vMax = ssGetSFcnParam(S, 3);
	if(!mx_vMax || mxGetNumberOfElements(mx_vMax) != n_axes)
		throw std::runtime_error("Invalid number of elements for parameter vMax");
	
	const mxArray * mx_uMin = ssGetSFcnParam(S, 4);
	if(!mx_uMin || mxGetNumberOfElements(mx_uMin) != n_axes)
		throw std::runtime_error("Invalid number of elements for parameter uMin");
	
	const mxArray * mx_uMax = ssGetSFcnParam(S, 5);
	if(!mx_uMax || mxGetNumberOfElements(mx_uMax) != n_axes)
		throw std::runtime_error("Invalid number of elements for parameter uMax");
	
	VectorXd full_x_min(2 * n_axes), full_x_max(2 * n_axes);
	full_x_min << Map<VectorXd>(mxGetPr(mx_qMin), n_axes), Map<VectorXd>(mxGetPr(mx_vMin), n_axes);
	full_x_max << Map<VectorXd>(mxGetPr(mx_qMax), n_axes), Map<VectorXd>(mxGetPr(mx_vMax), n_axes);
	
	const Map<VectorXd> full_u_min(mxGetPr(mx_uMin), n_axes);
	const Map<VectorXd> full_u_max(mxGetPr(mx_uMax), n_axes);

	VectorXd x_min(controller->nX()), x_max(controller->nX()), u_min(controller->nU()), u_max(controller->nU());
	x_min << full_x_min.row(0), full_x_min.row(n_axes);
	x_max << full_x_max.row(0), full_x_max.row(n_axes);
	u_min << full_u_min.row(0);
	u_max << full_u_max.row(0);
	
	controller->setXMin(x_min);
	controller->setXMax(x_max);
	controller->setUMin(u_min);
	controller->setUMax(u_max);

	// Setting the washout position and washout factor from block parameters.
	const mxArray * mx_washoutPos = ssGetSFcnParam(S, 6);
	if (!mx_washoutPos || mxGetNumberOfElements(mx_washoutPos) != n_axes)
		throw std::runtime_error("Invalid number of elements for parameter washoutPos");

	const mxArray * mx_washoutFactor = ssGetSFcnParam(S, 7);
	if (mxGetNumberOfElements(mx_washoutFactor) != 1)
		throw std::runtime_error("washoutFactor must be a scalar");

	controller->setWashoutPosition(Map<VectorXd>(mxGetPr(mx_washoutPos), n_axes).row(0));
	controller->setWashoutFactor(mxGetScalar(mx_washoutFactor));

	// Initialize weighting matrices.
	const mxArray * mx_errorWeight = ssGetSFcnParam(S, 9);
	if (mxGetNumberOfElements(mx_errorWeight) != controller->nY())
		throw std::runtime_error("Invalid number of elements for parameter errorWeight");
	controller->setErrorWeight(Map<VectorXd>(mxGetPr(mx_errorWeight), mxGetNumberOfElements(mx_errorWeight)));

	// Initialize terminal state bounds.
	const mxArray * mx_finalVelocityZero = ssGetSFcnParam(S, 10);
	if (mxGetNumberOfElements(mx_finalVelocityZero) != 1)
		throw std::runtime_error("finalVelocityZero must be a scalar");

	if (mxGetScalar(mx_finalVelocityZero) != 0.)
	{
		VectorXd term_x_min = x_min, term_x_max = x_max;
		term_x_min.bottomRows(1).fill(0.);
		term_x_max.bottomRows(1).fill(0.);

		controller->setTerminalXMin(term_x_min);
		controller->setTerminalXMax(term_x_max);
	}
	else
	{
		controller->setTerminalXMin(x_min);
		controller->setTerminalXMax(x_max);
	}
		
	std::ostringstream os;
	os << "Controller limits set to:" << std::endl
		<< "xMin =\t" << controller->getXMin().transpose() << std::endl
		<< "xMax =\t" << controller->getXMax().transpose() << std::endl
		<< "terminalXMin =\t" << controller->getTerminalXMin().transpose() << std::endl
		<< "terminalXMax =\t" << controller->getTerminalXMax().transpose() << std::endl
		<< "uMin =\t" << controller->getUMin().transpose() << std::endl
		<< "uMax =\t" << controller->getUMax().transpose() << std::endl
		<< "Washout position:" << std::endl
		<< "washoutPos =\t" << controller->getWashoutPosition().transpose() << std::endl
		<< "washoutFactor =\t" << controller->getWashoutFactor() << std::endl
		<< "nIntervals =\t" << controller->getNumberOfIntervals() << std::endl
		<< "errorWeight =\t" << controller->getErrorWeight().transpose() << std::endl;
	
	ssPrintf(os.str().c_str());

	setController(S, controller.release());
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

	Map<const VectorXd> x((const real_T*)ssGetInputPortSignal(S, 1), ssGetCurrentInputPortWidth(S, 1));
	Map<VectorXd> u((real_T *)ssGetOutputPortRealSignal(S, 0), ssGetCurrentOutputPortWidth(S, 0));

	auto * controller = getController(S);

	// Set platform's current state.
	platform->setFullState(x);

	// Make a "reduced" state vector, containing only axis 0 position and acceleration.
	VectorXd x0(controller->nX());
	x0 << x.row(0), x.row(n_axes);

	Map<const VectorXd> y0_ref((const real_T*)ssGetInputPortSignal(S, 0), ssGetCurrentInputPortWidth(S, 0));
	//Map<const VectorXd> x((const real_T*)ssGetInputPortSignal(S, 1), ssGetCurrentInputPortWidth(S, 1));
	//Map<VectorXd> u((real_T *)ssGetOutputPortRealSignal(S, 0), ssGetCurrentOutputPortWidth(S, 0));

	// Initialize working point if needed.
	if (!getWorkingPointInitialized(S))
	{
		controller->InitWorkingPoint(x0);
		setWorkingPointInitialized(S, true);

		std::ostringstream os;
		os << "mdlOutputs(): MPC controller working point initialized to " << x0.transpose() << std::endl;
		ssPrintf(os.str().c_str());
	}

	// Initialize new reference.
	// Assume constant reference output for all prediction horizon.
	MatrixXd y_ref(y0_ref.size(), controller->getNumberOfIntervals());
	for (unsigned i = 0; i < controller->getNumberOfIntervals(); ++i)
		y_ref.col(i) = y0_ref;

	try
	{
		controller->setReference(y_ref);
		controller->EmbedInitialValue(x0);
		
		/** Solve QP */
		controller->Solve();

// 		controller->PrintQP_MATLAB(log_stream);
// 		log_stream << std::flush;

		// Copy u[0] to output.
		u.setConstant(0.);
		u.row(0) = controller->getWorkingU(0);

		// Prepare for the next step.
		controller->PrepareForNext();
	}
	catch (const camels::CondensingSolverSolveException& e)
	{
		{
			std::ofstream os("failed_qp.m");
			controller->PrintQP_MATLAB(os);
			controller->PrintWorkingPoint_MATLAB(os, "wp");
			e.getCondensedQP().Print_MATLAB("cond_qp", os);
		}

		error_status = e.what();
		ssSetErrorStatus(S, error_status.c_str());
	}
	catch (const std::runtime_error& e)
	{
		{
			std::ofstream os("failed_qp.m");
			controller->PrintQP_MATLAB(os);
		}

		error_status = e.what();
		ssSetErrorStatus(S, error_status.c_str());
	}
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
}

/* Function: mdlTerminate =====================================================
 * Abstract:
 *    In this function, you should perform any actions that are necessary
 *    at the termination of a simulation.  For example, if memory was
 *    allocated in mdlStart, this is the place to free it.
 */
static void mdlTerminate(SimStruct *S)
{
	delete getController(S);
}

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif


