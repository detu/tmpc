/* This function was automatically generated by CasADi */
#ifdef __cplusplus
extern "C" {
#endif

#ifndef real_t
#define real_t double
#define to_double(x) (double) x
#define to_int(x) (int) x
#endif /* real_t */

int pendulum_ode(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
void pendulum_ode_incref(void);
void pendulum_ode_decref(void);
int pendulum_ode_n_in(void);
int pendulum_ode_n_out(void);
const char* pendulum_ode_name_in(int i);
const char* pendulum_ode_name_out(int i);
const int* pendulum_ode_sparsity_in(int i);
const int* pendulum_ode_sparsity_out(int i);
int pendulum_ode_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
#ifdef __cplusplus
} /* extern "C" */
#endif
