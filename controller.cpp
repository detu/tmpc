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

#include <qpDUNES.h>
#include <Eigen/Dense>

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
#include <array>
#include <vector>

class System
{
public:
	unsigned getInputDim() const
	{
		return _inputDim;
	}

	unsigned getStateDim() const
	{
		return _stateDim;
	}

	unsigned getOutputDim() const
	{
		return _outputDim;
	}

protected:
	System(unsigned nu, unsigned nx, unsigned ny)
		: _inputDim(nu), _stateDim(nx), _outputDim(ny)
	{
	}

private:
	unsigned _inputDim;
	unsigned _stateDim;
	unsigned _outputDim;
};

class MotionLimits
{
public:
	MotionLimits(double q_min, double q_max, double v_min, double v_max, double u_min, double u_max)
		: _qMin(q_min), _qMax(q_max), _vMin(v_min), _vMax(v_max), _uMin(u_min), _uMax(u_max)
	{
		if (!(q_min < q_max && v_min < v_max && u_min < u_max))
			throw std::invalid_argument("MotionLimits::MotionLimits(): invalid motion limits (min must be less than max)");
	}

private:
	double _qMin;
	double _qMax;
	double _vMin;
	double _vMax;
	double _uMin;
	double _uMax;
};

class MotionPlatform : public System
{
public:
	MotionPlatform(const std::vector<MotionLimits> limits)
		: System(static_cast<unsigned>(limits.size()), static_cast<unsigned>(2 * limits.size()), 6)
		, _gravity(0., 0., -9.81)
	{
	}

	unsigned getNumberOfAxes() const
	{
		return static_cast<unsigned>(_axesLimits.size());
	}

	const std::vector<MotionLimits>& getAxesLimits() const
	{
		return _axesLimits;
	}

	const Eigen::Vector3d& getGravity() const
	{
		return _gravity;
	}

	void setGravity(const Eigen::Vector3d& g)
	{
		_gravity = g;
	}

	typedef Eigen::Map<Eigen::VectorXd> VectorMap;
	typedef Eigen::Map<Eigen::MatrixXd, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>> MatrixMap;

	virtual void Output(const VectorMap& x, const VectorMap& u, VectorMap& y, MatrixMap& C, MatrixMap& D) const = 0;

private:
	std::vector<MotionLimits> _axesLimits;
	Eigen::Vector3d _gravity;
};

class MotionPlatformX : public MotionPlatform
{
public:
	MotionPlatformX()
		: MotionPlatform(DefaultAxesLimits())
	{
	}

	void Output(const VectorMap& x, const VectorMap& u, VectorMap& y, MatrixMap& C, MatrixMap& D) const override
	{
		if (!(x.size() == getStateDim() && u.size() == getInputDim() && y.size() == getOutputDim() 
			&& C.rows() == getOutputDim() && C.cols() == getStateDim() && D.rows() == getOutputDim() && D.cols() == getInputDim()))
			throw std::invalid_argument("MotionPlatformX::Output() input argument dimension mismatch");

		// Specific force
		y.block<3, 1>(0, 0) = getGravity();
		y(0) -= u(0);

		// Rotational velocity
		y.block<3, 1>(3, 0).fill(0.);

		// C = dy/dx
		C.fill(0.);

		// D = dy/du
		D.fill(0.);
		D(0, 0) = -1;
	}

private:
	static std::vector<MotionLimits> DefaultAxesLimits()
	{
		std::vector<MotionLimits> limits;
		limits.push_back(MotionLimits(-1., 1., -1., 1., -1., 1.));
		return limits;
	}
};

// Motion platform to be used.
MotionPlatformX platform;

// number of control intervals
const unsigned Nt = 3; 

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
// 	 // G = [dy/dx, dy/du]
// 	 Eigen::Matrix<double, Ny, Nz> G;
// 	 G << 0. << 0. << 1.;
// 
// 	 // Initializing matrices.
// 	 // H stores Nt matrices of size Nz x Nz and 1 matrix of size Nx x Nx.
// 	 std::array<double, Nz * Nz * Nt + Nx * Nx> H;
// 
// 	 // g stores Nt vectors of size Nz and 1 vector of size Nx
// 	 std::array<double, Nz * Nt + Nx> g;
// 
// 	 for (unsigned i = 0; i < Nt; ++i)
// 	 {
// 		 Eigen::Map<Eigen::Matrix<double, Nz, Nz, Eigen::RowMajor>>(H.data() + Nz * Nz * i) = G.transpose() * G;
// 		 Eigen::Map<Eigen::Matrix<double, Nz, 1>>(g.data() + Nz * i) = G.transpose() * G;
// 	 }
// 
//      auto qpData = reinterpret_cast<qpData_t *>(ssGetPWorkValue(S, 0));
//      
//      /** Initial MPC data setup: components not given here are set to zero (if applicable)
// 	 *      instead of passing g, D, zLow, zUpp, one can also just pass NULL pointers (0) */
// 	 
// 	 return_t statusFlag = qpDUNES_init(qpData, H, g, C, c, z_min.data(), z_max.data(), 0, 0, 0);
// 	 if (statusFlag != QPDUNES_OK) 
//      {
//   		 ssSetErrorStatus(S, "qpDUNES_init() failed.");
// 		 return;
// 	 }
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
    const real_T   *y_ref  = (const real_T*) ssGetInputPortSignal(S,0);
    const real_T   *x  = (const real_T*) ssGetInputPortSignal(S,1);
    real_T         *u  = (real_T *)ssGetOutputPortRealSignal(S,0);

    //S_Outputs_wrapper(y_ref, x, u, xD, S);
	u[0] = -y_ref[0];
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


