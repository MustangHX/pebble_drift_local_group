//
//  main.c
//  pebble_size_drift
//
//  Created by Xiao Hu on 5/16/15.
//
//

#include <stdio.h>
#include <math.h>
#include "ex_func.h"
#include "global_ex.h"
#include "global_var.h"
double pp_vr_tau[2]={0.0};
int main(){
	int i,j,k,n,i_new,j_new,offset_time=0,num_step=0,tot_num_step=10000,check;
	double AREA,dr=size_ring,a_pb1,a_pb2,a_max,vol_plus,delta_r,frac,tau,vr0;
	double coag_eff=1.0,tot_mass=0.0;
	FILE *fp,*fp2;
	char outname[256];
	char output_peb[256];
	n=0;
   // group();
	printf("PPPPeEEEBBBB%.16f\n",density(1.0));
	drift(1.0,drag_group(1.0,0.01));
        printf("PPPPeEEEBBBB%.12f\n",drag_group(1.0,0.01));
/*	for(i=0;i<100;i++){
		for(j=0;j<100-1;j++){
		pebble[i][j]=pebble[i][j+1]*velocity[i][j+1]*dt/dr;
		}
	}*/
	Init();
        fp=fopen("size_chart.txt","w");
        for(i=0;i<peb_size_num;i++){
                fprintf(fp,"%f\n",peb_map[0].size[i]);
	        }
        fclose(fp);
	        fp=fopen("rad_chart.txt","w");
        for(i=0;i<ring_num;i++){
                fprintf(fp,"%f\n",peb_map[i].rad);
	        }
        fclose(fp);

	num_step=0;
	while (num_step<tot_num_step)
	{
	if(num_step==0){
        sprintf(outname,"out_sigma%d.txt",num_step);
        fp=fopen(outname,"w");
        for(i=0;i<ring_num;i++){
        for(j=0;j<peb_size_num;j++){
        fprintf(fp,"%e\t",peb_map[i].surf_dens[j]);
	        }
        fprintf(fp,"\n");
	}
        fclose(fp);
	}
	for(i=ring_num-1;i>-1;i--){
	for(j=0;j<peb_size_num;j++){
		peb_map[i].vr[j]=vr_estimate(peb_map[i].rad,peb_map[i].size[j],pp_vr_tau);
		a_pb1=peb_map[i].size[j];
		tau=pp_vr_tau[1];
		vr0=peb_map[i].vr[j];
		delta_r=vr0*dt*TUNIT/LUNIT;
		a_max=2.55*mean_path(peb_map[i].rad);
		vol_plus=1.0*M_PI*a_pb1*a_pb1*sqrt(vr0*vr0+0.25*tau*vr0*tau*vr0)*dt*TUNIT;
		if (a_pb1 > 2.55*mean_path(peb_map[i].rad)) 
		{vol_plus=0.0;
			//printf("STOKES LIMIT\n");
			check=1;
		}
		else check=0;
		a_pb2=pow(((vol_plus*coag_eff*density(peb_map[i].rad)/rho_peb0+4.0/3.0*M_PI*a_pb1*a_pb1*a_pb1)*3.0/4.0/M_PI),0.33333333333333333);
		if(a_pb2 >a_max && check==0){
			a_pb2=2.55*mean_path(peb_map[i].rad);
		}
		i_new=i-floor(delta_r/dr)-1;
		if(i_new<0) i_new=0;
		j_new=floor(log10(a_pb2/0.1)/size_step);
		if(log10(a_pb2/0.1)/size_step-j_new*1.0 >=0.5) j_new+=1;
		if(j_new <0) j_new=0;
		else if(j_new >peb_size_num-1) {
			j_new=peb_size_num-1;
//			printf("maximum size %e\t%e\t%e\n",2.55*mean_path(peb_map[i].rad),vr0,tau);
		}
		frac=delta_r/dr-floor(delta_r/dr);
		peb_map[i_new].mass_in[j_new]+=frac*peb_map[i].mass_out[j];
		if(delta_r/dr>1.0){
		peb_map[i_new+1].mass_in[j_new]+=(1.0-frac)*peb_map[i].mass_out[j];
		peb_map[i].mass_out[j]-=peb_map[i].mass_out[j];
		}
		else {
//				peb_map[i].mass_in[j_new]+=(1.0-frac)*peb_map[i].mass_out[j];
                peb_map[i].mass_out[j]-=frac*peb_map[i].mass_out[j];
		}
	}
	}
	for(i=ring_num-1;i>-1;i--){
        for(j=0;j<peb_size_num;j++){
		peb_map[i].mass_out[j]+=peb_map[i].mass_in[j];
		peb_map[i].mass_in[j]=0.0;
		AREA=M_PI*((peb_map[i].rad+dr/2.0)*(peb_map[i].rad+dr/2.0)-(peb_map[i].rad-dr/2.0)*(peb_map[i].rad-dr/2.0))*LUNIT*LUNIT;
		if(i==ring_num-1 && j>=0 && 0) {
			peb_map[i].mass_out[j]=0.2*AREA*dust_budget[i].surf_dens*exp(-1.0*peb_map[i].size[j]/0.1)*exp(0.0*num_step/100);
		}
		else if(i<ring_num-1 && j<10 && 0){
                        peb_map[i].mass_out[j]=0.1*AREA*dust_budget[i].surf_dens*exp(-1.0*peb_map[i].size[j]/0.1)*exp(0.0*num_step/100);
		}

		else peb_map[i].surf_dens[j]+=1e-200;
		peb_map[i].surf_dens[j]=peb_map[i].mass_out[j]/AREA;
	}
	}
	//frag();
	num_step++;
	if(num_step%10==0){
		tot_mass=0.0;
        sprintf(outname,"out_sigma%d.txt",num_step);
        fp=fopen(outname,"w");
		fp2=fopen("mass_check.txt","a+");
        for(i=0;i<ring_num;i++){
        for(j=0;j<peb_size_num;j++){
		if(peb_map[i].surf_dens[j] < 1e-200) peb_map[i].surf_dens[j]=1e-200;
        fprintf(fp,"%e\t",peb_map[i].surf_dens[j]);
		tot_mass+=peb_map[i].mass_out[j];
	}
        fprintf(fp,"\n");		
	}
		fprintf(fp2,"%e\n",tot_mass);
	fclose(fp);
	fclose(fp2);
	printf("%f finished\r",1.0*num_step/tot_num_step);
	}

	}



return 0;
}
