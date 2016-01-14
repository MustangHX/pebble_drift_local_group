#include <stdio.h>
#include "ex_func.h"
#include "global_ex.h"

void check_disk(double r){
	printf("Properties at 1 AU:\n");
	printf("Surf_dens=%g\nTemperature=%f\nMidplane-density=%g\nSound_speed=%g\n\
H/R=%g\n",Sigma(r),temperature(r),density(r),sound_sp(r),height(r)/r/LUNIT);
}
