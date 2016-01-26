#include "ex_func.h"
#include "global_ex.h"
#include "global_var.h"
#include<math.h>
#include<stdio.h>
void frag(){
int i,j,k,jl,js;
double size_lim,frac,AREA,dr,delta_mass;
dr=size_ring;
for(i=ring_num-1;i>=0;i--){
	dr=peb_map[i].dr;
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

double coag_kernel(double a_pb1,double a_pb2,double delta_v,double rho1, double rho2, double dt0, double rad){
	double cross_section, path,n_density;
	if(a_pb1 >=a_pb2) n_density=rho2/(rho_peb*4.0*M_PI*a_pb1*a_pb1*a_pb1/3.0);
	else n_density=rho1/(rho_peb*4.0*M_PI*a_pb2*a_pb2*a_pb2/3.0)*0.85;//disk within -h to h
	cross_section=M_PI*(a_pb1*a_pb1+a_pb2*a_pb2);
	cross_section=M_PI*(a_pb1+a_pb2)*(a_pb1+a_pb2);
	//delta_v=fabs(v1-v2);
	path=cross_section*n_density*delta_v*dt0*TUNIT*1e0;
	if(path>1.0)	printf("path=%e\t%e\t%f\n",path,n_density,exp(-1.0*path));
	return exp(-1.0*path);
}
void coagulation(double dt0){
        int i,j,jj,i_new,j_new;
	double a_pb1,a_pb11,a_pb2,a_pb22,a_pb3,vr0,v1,v2,delta_v,AREA,dr,a_max,dt_new,dt1,sub_time1;
	double tau,vol_plus,frac,frac_s,coag_eff,ratio_size,ring_mass_before,ring_mass_after,old_sigma,ratio_sigma=1.0,dust_flow=0.0,mass0;
	double peb_flow[peb_size_num]={0.0};
	
	for (i=ring_num-1;i>=10;i--){
		dr=peb_map[i].dr;
		for(j=peb_size_num-1;j>=1;j--){
			if(peb_map[i].surf_dens[j] < 1e-5) continue;
			for(jj=j-1;jj>=0;jj--){
			a_pb1=peb_map[i].size_med[j];
			a_pb2=peb_map[i].size_med[jj];
			v1=peb_map[i].vr[j]-peb_map[i].vr[jj];
                        v2=peb_map[i].vt[j]-peb_map[i].vt[jj];
			delta_v=sqrt(v1*v1+v2*v2);

			a_pb3=pow(a_pb1*a_pb1*a_pb1+a_pb2*a_pb2*a_pb2,1.0/3.0);
			j_new=floor(log10(a_pb3/size_min)/size_step);
			if(j_new>=peb_size_num) j_new=peb_size_num-1;
			
			frac_s=1.0-coag_kernel(a_pb1,a_pb2,delta_v,peb_map[i].rho[j],peb_map[i].rho[jj],dt0,peb_map[i].rad_med);
			if(frac_s> 0.8) printf("rho=%e\t%e\t%f\t%f\t%f\t%f\t%f\n",peb_map[i].rho[j],peb_map[i].rho[jj],frac_s,peb_map[i].rad_med,peb_map[i].size_med[j],peb_map[i].size_med[jj],delta_v);
			//			frac_s=0.5;
			//peb_map[i].mass_in[j]-=peb_map[i].mass_in[j]*frac_s;
			mass0=peb_map[i].mass_out[j];
			peb_map[i].mass_in[j_new]+=mass0*frac_s*(pow(a_pb2/a_pb1,3)+1.0);
			peb_map[i].mass_out[jj]-=mass0*frac_s*pow(a_pb2/a_pb1,3);
			peb_map[i].mass_out[j]-=mass0*frac_s;
		//	peb_map[i].surf_dens[j]=peb_map[i].mass_out[j]/peb_map[i].AREA;
                //	peb_map[i].rho[j]=peb_map[i].surf_dens[j]/sqrt(2.0*M_PI)/peb_map[i].hei[j];
		//	peb_map[i].surf_dens[jj]=peb_map[i].mass_out[jj]/peb_map[i].AREA;
                //	peb_map[i].rho[jj]=peb_map[i].surf_dens[jj]/sqrt(2.0*M_PI)/peb_map[i].hei[jj];
			}
		
		}
	}

        for(i=ring_num-1;i>=0;i--){
        for(j=0;j<peb_size_num;j++){
                peb_map[i].mass_out[j]+=peb_map[i].mass_in[j];
                peb_map[i].mass_in[j]=0.0;
                peb_map[i].surf_dens[j]=peb_map[i].mass_out[j]/peb_map[i].AREA;
                peb_map[i].rho[j]=peb_map[i].surf_dens[j]/sqrt(2.0*M_PI)/peb_map[i].hei[j];
        }
        }


}
