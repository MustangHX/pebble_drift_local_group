#include "ex_func.h"
#include "global_ex.h"
#include "global_var.h"
#include<math.h>
void frag(){
int i,j,k,jl,js;
double size_lim,frac,AREA,dr,delta_mass;
dr=size_ring;
for(i=ring_num-1;i>=0;i--){
	for(j=peb_size_num-1;j>=0;j--){
		size_lim=2.25*mean_path(peb_map[i].rad);		
		if(peb_map[i].size[j] > 0.8*size_lim){
			jl=j;
			js=floor(log10(peb_map[i].size[jl]/0.1/2.0)/size_step);
//			js=jl;
			if(js <0) js=0;
			frac=peb_map[i].size[jl]/size_lim;
			frac=frac*0.9;
			if (frac >1.0) frac=1.0;
			//frac=0.0;
			delta_mass=peb_map[i].mass_out[j]*frac;
			peb_map[i].mass_out[jl]-=delta_mass;
			peb_map[i].mass_out[js]+=delta_mass;
	}
}
}

for(i=0;i<ring_num;i++){
	AREA=M_PI*((peb_map[i].rad+dr/2.0)*(peb_map[i].rad+dr/2.0)-(peb_map[i].rad-dr/2.0)*(peb_map[i].rad-dr/2.0))*LUNIT*LUNIT;
	for(j=0;j<peb_size_num;j++){
		peb_map[i].surf_dens[j]=(peb_map[i].mass_out[j]+peb_map[i].mass_in[j])/AREA;
	}
}
}

