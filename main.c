//
//  main.c
//  pebble_size_drift
//
//  Created by Xiao Hu on 5/16/15.
//
//
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "ex_func.h"
#include "global_ex.h"
#include "global_var.h"
double pp_vr_tau0[2]={0.0};
int main(argc, argv)
	int argc;
	char *argv[];
{
	int i,j,k,n,i_new,j_new,offset_time=0,num_step=0,tot_num_step=(int)(time_yr*1.0/init_step),check,NbRestart;
	double AREA,dr=size_ring,a_pb1,a_max,vol_plus,delta_r,delta_size,d_size,ratio_size,frac,frac_s,tau,vr0,vol1,vol2,mass_flow_inner;
	double coag_eff=1.0,tot_mass=0.0,out_source=0.0,a_p,r0,dt=init_step,time_sum=0.0,dt2,tot_mass_dust;
	
	FILE *fp,*fp2,*fp3,*fp4;
	char outname[256], outname2[256];
	char output_peb[256];
	int Restarting = 0,grow_method=3,drift_method=1;;
	printf("%s\n",argv);
	for(i=0; i< argc; i++){
		printf("%c",argv[i]);
		if (strchr (argv[i], 's')) {
		Restarting = 1;
		i++;
		NbRestart = atoi(argv[i]);
		}
		if (NbRestart < 0) printf("Incorrect restart number\n");
	}						      
	n=0;
   // group();
	printf("PPPPeEEEBBBB%.16f\n",density(1.0));
	//drift(1.0,drag_group(1.0,0.01));
        printf("PPPPeEEEBBBB%.12f\n",drag_group(1.0,0.01));
	Init2();
	check_disk(1.0);
	if(Restarting == 1){
		Restart(NbRestart);
		time_sum+=NbRestart;
	}
        fp=fopen("size_chart.txt","w");
        for(i=0;i<peb_size_num;i++){
                fprintf(fp,"%f\n",peb_map[0].size_med[i]);
	        }
        fclose(fp);
	        fp=fopen("rad_chart.txt","w");
        for(i=0;i<ring_num;i++){
                fprintf(fp,"%f\n",peb_map[i].rad+peb_map[i].dr/2.0);
	        }
        fclose(fp);
	out_source=0.0;
	for(j=0;j<peb_size_num;j++){
		i=ring_num-1;
		AREA=peb_map[i].AREA;
		if(peb_map[i].size_med[j]<peb_size_lim){
		out_source+=AREA*exp(-1.0*size_slope*size_slope*peb_map[i].size_med[j]);
		}
	}

	fp=fopen("1mm.txt","w");
	num_step=0;
	a_p=0.1;
	r0=10.375;
        fprintf(fp,"%e\t%e\n",r0,a_p);
			
	while(r0>1.0){
	vr0=vr_estimate(r0,a_p,pp_vr_tau0);
	tau=pp_vr_tau0[1];
	r0-=vr0*dt*TUNIT/LUNIT;
	a_max=2.25*mean_path(r0);
	vol_plus=1.0*M_PI*a_p*a_p*sqrt(vr0*vr0+0.25*tau*vr0*tau*vr0)*dt*TUNIT;
	if(a_p > a_max) {check=1; vol_plus=0.0;}
	else check=0;
	a_p=pow(((vol_plus*coag_eff*0.1*density(r0)/rho_peb+4.0/3.0*M_PI*a_p*a_p*a_p)*3.0/4.0/M_PI),0.33333333333333333);
	if(check==0 && a_p>=a_max)       a_p=a_max;
	fprintf(fp,"%e\t%e\n",r0,a_p);
	}
	fclose(fp);

	fp=fopen("gas_vr.txt","w");
	for(i=0;i<ring_num;i++){
	fprintf(fp,"%f\t%g\t%g\n",dust_budget[i].rad_med,vr_gas(dust_budget[i].rad_med),vr_gas(dust_budget[i].rad_med)*TUNIT/LUNIT);
	}
	fclose(fp);
	fp=fopen("gas_dust_density.txt","w");
        for(i=0;i<ring_num;i++){
        fprintf(fp,"%f\t%g\t%g\n",dust_budget[i].rad_med,Sigma(dust_budget[i].rad_med),dust_budget[i].surf_dens[0]);
        }
        fclose(fp);
	stokes_size();
	tau_unity();
	num_step=0;
	if(Restarting == 1){
		num_step=NbRestart;
	}

	//start time sequence
		
	while (time_sum<tot_num_step*1.0)
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
	sprintf(outname2,"dust_sigma%d.txt",(int)time_sum);
	fp3=fopen(outname2,"w");
	for(i=0;i<ring_num;i++){
		tot_mass+=dust_budget[i].mass_out;
		fprintf(fp3,"%e\t%e\n",dust_budget[i].rad,dust_budget[i].surf_dens[0]);
	}
	fclose(fp3);
	}
	if(grow_method==1){grow_1();}
	if(grow_method==3){
		dt2=dt;
		dt=grow_3b_ada(dt2);
		//printf("main dt=%f\tdt2=%f\n",dt,dt2);
		//dt=1.0;
		if(dt<0.0) { 
		printf("Actual time step count:%d\t dt=%f\n",num_step,dt);
		return 0;
		}
	}
	dust_evolve(dt);
	//coagulation(dt);
	mass_flow_inner=0.0;
	for(i=ring_num-1;i>-1;i--){
        for(j=0;j<peb_size_num;j++){
                if(i==0) mass_flow_inner+=peb_map[i].mass_out[j];

		peb_map[i].mass_out[j]+=peb_map[i].mass_in[j];
		peb_map[i].mass_in[j]=0.0;
		//AREA=M_PI*((peb_map[i].rad+dr/2.0)*(peb_map[i].rad+dr/2.0)-(peb_map[i].rad-dr/2.0)*(peb_map[i].rad-dr/2.0))*LUNIT*LUNIT;
		AREA=peb_map[i].AREA;
		
		if(i==ring_num-1 && 1 && 1) {
			//peb_map[i].mass_out[j]=0.2*AREA*dust_budget[i].surf_dens*exp(-1.0*peb_map[i].size[j]/0.1)*exp(0.0*num_step/100);
		if(peb_map[i].size_med[j]<peb_size_lim){
			peb_map[i].mass_out[j]+=1.0*dust_gas*peb_dust*AREA*MUNIT*mdot*dt*exp(-1.0*size_slope*peb_map[i].size_med[j])/out_source;
			}
		}
		else if(i<ring_num-1 && j<10 && 0){
                        peb_map[i].mass_out[j]=0.1*AREA*dust_budget[i].surf_dens[0]*exp(-1.0*size_slope*peb_map[i].size_med[j])*exp(0.0*num_step/100);
		}

		else peb_map[i].surf_dens[j]+=peb_low_lim/10;
		peb_map[i].surf_dens[j]=peb_map[i].mass_out[j]/AREA;
		peb_map[i].rho[j]=peb_map[i].surf_dens[j]/sqrt(2.0*M_PI)/peb_map[i].hei[j];
		//printf("peb_rho=%e\t%e\t",peb_map[i].rho[j],peb_map[i].surf_dens[j]);
	}
	}
	//frag();
	num_step++;
	time_sum+=dt;
	if(time_sum-floor(time_sum)<0.00001 && ((int)(time_sum))%((int)(outp_step))==0){
		tot_mass=0.0;
		tot_mass_dust=0.0;
        sprintf(outname,"out_sigma%d.txt",(int)time_sum);
        fp=fopen(outname,"w");
	sprintf(outname2,"dust_sigma%d.txt",(int)time_sum);
	fp3=fopen(outname2,"w");
	fp4=fopen("mass_flow_inner_ring.txt","a+");
	fp2=fopen("mass_check.txt","a+");
        for(i=0;i<ring_num;i++){
        for(j=0;j<peb_size_num;j++){
		if(peb_map[i].surf_dens[j] < peb_low_lim) peb_map[i].surf_dens[j]=peb_low_lim;
        fprintf(fp,"%e\t",peb_map[i].surf_dens[j]);
		tot_mass+=peb_map[i].mass_out[j];
	}
        fprintf(fp,"\n");		
	}
	for(i=0;i<ring_num;i++){
		tot_mass_dust+=dust_budget[i].mass_out;
	fprintf(fp3,"%e\t%e\n",dust_budget[i].rad,dust_budget[i].surf_dens[0]);
	}
	fprintf(fp2,"%2.20g\t%2.20g\t%2.20g\n",tot_mass,tot_mass_dust,tot_mass+tot_mass_dust);
	fprintf(fp4,"%f\t%2.20g\n",time_sum,mass_flow_inner);
	fclose(fp);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);
	printf("%f finished\r",time_sum/(tot_num_step*1.0));
	printf("Actual time step count:%d\t dt=%f\t time=%f\n",num_step,dt,time_sum);
	}

	}
	
	printf("Actual time step count:%d\t dt=%f\n",num_step,dt);


return 0;
}
