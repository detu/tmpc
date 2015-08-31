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

#define NPARAMS              8

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
#include "CyberMotion.hpp"
#include "MPC_Controller.hpp"

#include <memory>
#include <fstream>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrix;

// Motion platform to use
auto platform = std::make_shared<PLATFORM_TYPE>();

// Number of prediction intervals
const unsigned Np = N_PREDICTION;

// Number of control intervals
#ifndef N_CONTROL
	#define N_CONTROL N_PREDICTION
#endif
const unsigned Nc = N_CONTROL;

// Log stream.
std::ofstream log_stream;

// Error status for Simulink
std::string error_status;

mpmc::MPC_Controller * getController(SimStruct * S)
{
	return reinterpret_cast<mpmc::MPC_Controller *>(ssGetPWorkValue(S, 0));
}

void setController(SimStruct * S, mpmc::MPC_Controller * c)
{
	ssSetPWorkValue(S, 0, c);
}

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
    ssSetInputPortWidth(S,  0, platform->getOutputDim()); /* */
    ssSetInputPortDataType(S, 0, SS_DOUBLE);
    ssSetInputPortComplexSignal(S,  0, INPUT_0_COMPLEX);
    ssSetInputPortDirectFeedThrough(S, 0, INPUT_0_FEEDTHROUGH);
    ssSetInputPortRequiredContiguous(S, 0, 1); /*direct input signal access*/

    /*Input Port 1 */
    ssSetInputPortWidth(S,  1, platform->getStateDim());
    ssSetInputPortDataType(S, 1, SS_DOUBLE);
    ssSetInputPortComplexSignal(S, 1, INPUT_1_COMPLEX);
    ssSetInputPortDirectFeedThrough(S, 1, INPUT_1_FEEDTHROUGH);
    ssSetInputPortRequiredContiguous(S, 1, 1); /*direct input signal access*/


    if (!ssSetNumOutputPorts(S, NUM_OUTPUTS)) return;
    ssSetOutputPortWidth(S, 0, platform->getInputDim());
    ssSetOutputPortDataType(S, 0, SS_DOUBLE);
    ssSetOutputPortComplexSignal(S, 0, OUTPUT_0_COMPLEX);
    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    
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
	 Initialize at current position, 0 velocity, 0 input.
	 */

	 log_stream.open("controller.log");
	 log_stream << "mdlInitializeConditions()" << std::endl;

	 try
	 {
		 auto const controller = getController(S);

		 for (unsigned i = 0; i < controller->getNumberOfIntervals(); ++i)
		 {
			 if (i < Nc)
				 controller->W(i) = OutputWeightingMatrix();
			 else
				 controller->W(i).setZero();
		 }

		 Eigen::VectorXd x0(controller->getStateDim());
		 x0.fill(0.);
		 platform->getDefaultAxesPosition(x0.data());

		 std::ostringstream msg;
		 msg << "mdlInitializeConditions(): initializing MPC controller at working point " << x0.transpose() << std::endl;
		 mexPrintf(msg.str().c_str());

		 controller->InitWorkingPoint(x0);
		 controller->PrintQP_MATLAB(log_stream);
		 log_stream << std::flush;
	 }	 
	 catch (const std::runtime_error& e)
	 {
		 error_status = e.what();
		 ssSetErrorStatus(S, error_status.c_str());
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
    using Eigen::Map;
    using Eigen::VectorXd;
    
    /** Initialize MPC_Controller */
	auto controller = std::make_unique<mpmc::MPC_Controller>(platform, SAMPLE_TIME_0, Np);
	controller->setLevenbergMarquardt(0.01);
    controller->setWashoutFactor(5);
    
    // Set motion limits from block parameters.
    const auto n_axes = platform->getNumberOfAxes();
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
    
    VectorXd x_min(2 * n_axes), x_max(2 * n_axes);
    x_min << Map<VectorXd>(mxGetPr(mx_qMin), n_axes), Map<VectorXd>(mxGetPr(mx_vMin), n_axes);
    x_max << Map<VectorXd>(mxGetPr(mx_qMax), n_axes), Map<VectorXd>(mxGetPr(mx_vMax), n_axes);
    
    const Map<VectorXd> u_min(mxGetPr(mx_uMin), n_axes);
    const Map<VectorXd> u_max(mxGetPr(mx_uMax), n_axes);
    
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

	controller->setWashoutPosition(Map<VectorXd>(mxGetPr(mx_washoutPos), n_axes));
	controller->setWashoutFactor(mxGetScalar(mx_washoutFactor));
    
    std::ostringstream os;
	os << "Controller limits set to:" << std::endl
		<< "xMin = " << controller->getXMin().transpose() << std::endl
		<< "xMax = " << controller->getXMax().transpose() << std::endl
		<< "uMin = " << controller->getUMin().transpose() << std::endl
		<< "uMax = " << controller->getUMax().transpose() << std::endl
		<< "Washout position:" << std::endl
		<< "washoutPos = " << controller->getWashoutPosition().transpose() << std::endl
		<< "washoutFactor = " << controller->getWashoutFactor() << std::endl;
    
    mexPrintf(os.str().c_str());

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

	log_stream << "mdlOutputs()" << std::endl;

	const real_T * px = (const real_T*)ssGetInputPortSignal(S, 1);
	real_T * pu = (real_T *)ssGetOutputPortRealSignal(S, 0);

	Map<const VectorXd> y0_ref((const real_T*)ssGetInputPortSignal(S, 0), ssGetCurrentInputPortWidth(S, 0));
	//Map<const VectorXd> x((const real_T*)ssGetInputPortSignal(S, 1), ssGetCurrentInputPortWidth(S, 1));
	//Map<VectorXd> u((real_T *)ssGetOutputPortRealSignal(S, 0), ssGetCurrentOutputPortWidth(S, 0));

	auto * controller = getController(S);

	// Initialize new reference.
	// Assume constant reference output for all prediction horizon.
	MatrixXd y_ref(y0_ref.size(), controller->getNumberOfIntervals());
	for (unsigned i = 0; i < controller->getNumberOfIntervals(); ++i)
		y_ref.col(i) = y0_ref;

	try
	{
        controller->SetReference(y_ref.data());
		controller->EmbedInitialValue(px);
        
		/** Solve QP */
		controller->Solve();

// 		controller->PrintQP_MATLAB(log_stream);
// 		log_stream << std::flush;

		// Copy u[0] to output.
		controller->getWorkingU(0, pu);

		// Prepare for the next step.
		controller->PrepareForNext();
	}
	catch (const std::runtime_error& e)
	{
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


