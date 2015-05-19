

/*
 * Include Files
 *
 */
#include "simstruc.h"


/* %%%-SFUNWIZ_wrapper_includes_Changes_BEGIN --- EDIT HERE TO _END */
#include <math.h>
/* %%%-SFUNWIZ_wrapper_includes_Changes_END --- EDIT HERE TO _BEGIN */
#define u_width 1
#define y_width 1
/*
 * Create external references here.  
 *
 */
/* %%%-SFUNWIZ_wrapper_externs_Changes_BEGIN --- EDIT HERE TO _END */
/* extern double func(double a); */
/* %%%-SFUNWIZ_wrapper_externs_Changes_END --- EDIT HERE TO _BEGIN */

/*
 * Output functions
 *
 */
void S_Outputs_wrapper(const real_T *y_ref,
			const real_T *x,
			real_T *u,
			const real_T *xD,
			SimStruct *S)
{
      u[0] = y_ref[0]; 
}

/*
  * Updates function
  *
  */
void S_Update_wrapper(const real_T *y_ref,
			const real_T *x,
			const real_T *u,
			real_T *xD,
			SimStruct *S)
{
  /* %%%-SFUNWIZ_wrapper_Update_Changes_BEGIN --- EDIT HERE TO _END */
 
/* %%%-SFUNWIZ_wrapper_Update_Changes_END --- EDIT HERE TO _BEGIN */
}
