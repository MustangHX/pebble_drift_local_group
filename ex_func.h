#include "global_var.h"
double drift(double argc,double argc2);
double drag_group(double argc,double argc2);
double w_K(double argc);
double density(double argc);
double sound_sp(double argc);
double viscosity(double argc);
double mean_path( double argc);
double vr_estimate(double argc1,double argc2,double *argc3);
int drift_t(PEBBLE *argc1, double argc2, int argc3);
void Init();
void frag();
void Restart(int argc);
void grow_1();
void grow_3();
double grow_3b();
double grow_3b2();
void dust_evolve(double argc);
void check_disk(double argc);
double v_r_gas(double argc);
double height(double argc);
double coag_kernel(double a_pb1,double a_pb2,double delta_v,double rho1, double rho2, double dt0, double rad);
void coagulation(double argc);
double Sigma(double argc);
double yeta(double argc);
double v_K(double argc);
double w_K(double argc);
double vt_gas(double argc);
double density(double argc);
double temperature(double argc);
double viscosity(double argc);
double mean_path(double argc);
double sound_sp(double argc);

