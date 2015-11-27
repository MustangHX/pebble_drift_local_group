#include "global_ex.h"
#include "global_var.h"
#include "ex_func.h"
#include<math.h>
double pp_vr_tau[2]={0.0};
void grow_1(){
int i,j,k,n,i_new,j_new,offset_time=0,num_step=0,tot_num_step=(int)(time_yr*1.0/dt),check,NbRestart;
	double AREA,dr=size_ring,a_pb1,a_max,vol_plus,delta_r,delta_size,d_size,ratio_size,frac,frac_s,tau,vr0,vol1,vol2;
	double coag_eff=1.0,tot_mass=0.0,out_source=0.0,vol_plus_sum,dust_dens_sum,dust_rho_sum,old_dens,a_p,r0;
	double vol_plus_arr[peb_size_num],delta_mass[peb_size_num],a_pb2[peb_size_num],peb_mass[peb_size_num];
for(i=ring_num-1;i>-1;i--){
	vol_plus_sum=1e-20;
	dust_dens_sum=0.0;
	dust_rho_sum=0.0;
	AREA=M_PI*((peb_map[i].rad+dr/2.0)*(peb_map[i].rad+dr/2.0)-(peb_map[i].rad-dr/2.0)*(peb_map[i].rad-dr/2.0))*LUNIT*LUNIT;

	for(j=0;j<peb_size_num;j++){
	peb_map[i].vr[j]=vr_estimate(peb_map[i].rad,peb_map[i].size_med[j],pp_vr_tau);
	a_pb1=peb_map[i].size[j];
	tau=pp_vr_tau[1];
	vr0=peb_map[i].vr[j];
	delta_r=vr0*dt*TUNIT/LUNIT;
        a_max=2.25*mean_path(peb_map[i].rad);
	vol_plus_arr[j]=1.0*M_PI*a_pb1*a_pb1*sqrt(vr0*vr0+0.25*tau*vr0*tau*vr0)*dt*TUNIT;
	if(a_pb1>=a_max) vol_plus_arr[j]=0.0;
	a_pb2[j]=pow(((vol_plus_arr[j]*coag_eff*dust_budget[i].rho[j]/rho_peb0+4.0/3.0*M_PI*a_pb1*a_pb1*a_pb1)*3.0/4.0/M_PI),0.33333333333333333);
	if(a_pb1<a_max && a_pb2[j]>=a_max){
	a_pb2[j]=a_max;
	vol_plus_arr[j]=(pow(a_pb2[j],3.0)-pow(a_pb1,3.0))*4.0*M_PI/3.0*rho_peb0/dust_budget[i].rho[j];
	}
	peb_mass[j]=pow(peb_map[i].size[j],3)*M_PI*4.0/3.0*rho_peb0;
	vol_plus_sum+=vol_plus_arr[j]*peb_map[i].mass_out[j]/peb_mass[j];
	dust_dens_sum+=dust_budget[i].surf_dens[j];
	dust_rho_sum+=dust_budget[i].rho[j];
	a_pb2[j]=pow(((vol_plus_arr[j]*coag_eff*dust_budget[i].rho[j]/rho_peb0+4.0/3.0*M_PI*a_pb1*a_pb1*a_pb1)*3.0/4.0/M_PI),0.33333333333333333);
	}
	for(j=0;j<peb_size_num;j++){
        a_pb1=peb_map[i].size[j];
//	dust_budget[i].surf_dens[j]=vol_plus_arr[j]*peb_map[i].mass_out[j]*dust_dens_sum/vol_plus_sum/peb_mass[j];
	old_dens=dust_budget[i].surf_dens[j];
//	dust_budget[i].rho[j]=vol_plus_arr[j]*peb_map[i].mass_out[j]*dust_rho_sum/vol_plus_sum/peb_mass[j];
	a_pb2[j]=pow(((vol_plus_arr[j]*coag_eff*dust_budget[i].rho[j]/rho_peb0+4.0/3.0*M_PI*a_pb1*a_pb1*a_pb1)*3.0/4.0/M_PI),0.33333333333333333);
	a_max=2.25*mean_path(peb_map[i].rad);
	if(a_pb1<a_max && a_pb2[j]>a_max) a_pb2[j]=a_max;
	delta_mass[j]=(pow(a_pb2[j]/a_pb1,3)-1.0)*peb_map[i].mass_out[j];
/*	if(AREA*dust_budget[i].surf_dens[j]<=delta_mass[j]){
	dust_budget[i].surf_dens[j]=0.0;
	delta_mass[j]=AREA*dust_budget[i].surf_dens[j];
//	printf("%1.12g\t%1.12g\n",delta_mass[j],AREA*dust_budget[i].surf_dens[j]);
	}
	else dust_budget[i].surf_dens[j]-=delta_mass[j]/AREA;
	dust_budget[i].rho[j]=dust_budget[i].rho[j]*dust_budget[i].surf_dens[j]/old_dens;

	a_pb2[j]=pow(delta_mass[j]/peb_map[i].mass_out[j]+1.0,1.0/3.0)*a_pb1;
*/	}

	for(j=0;j<peb_size_num;j++){
/*		peb_map[i].vr[j]=vr_estimate(peb_map[i].rad,peb_map[i].size[j],pp_vr_tau);
		a_pb1=peb_map[i].size[j];
		tau=pp_vr_tau[1];
		vr0=peb_map[i].vr[j];
		delta_r=vr0*dt*TUNIT/LUNIT;
		a_max=2.25*mean_path(peb_map[i].rad);
		vol_plus=1.0*M_PI*a_pb1*a_pb1*sqrt(vr0*vr0+0.25*tau*vr0*tau*vr0)*dt*TUNIT;
		*/
		a_pb1=peb_map[i].size[j+1];
		if (a_pb1 >=a_max) 
		{vol_plus=0.0;
			//printf("STOKES LIMIT\n");
			check=1;
		}
		else check=0;
		//delta_mass=vol_plus*coag_eff*density(peb_map[i].rad);
		a_pb2[j]=pow(((vol_plus_arr[j]*coag_eff*dust_budget[i].rho[j]/rho_peb0+4.0/3.0*M_PI*a_pb1*a_pb1*a_pb1)*3.0/4.0/M_PI),0.33333333333333333);
		if(a_pb2[j] >a_max && check==0){
			a_pb2[j]=2.25*mean_path(peb_map[i].rad);
		}
		delta_mass[j]=(pow(a_pb2[j]/a_pb1,3)-1.0)*peb_map[i].mass_out[j];
		
		delta_r=peb_map[i].vr[j]*dt*TUNIT/LUNIT;
		i_new=i-floor(delta_r/dr)-1;
		if (a_pb1 >=a_max) {check=1;}
		if(i_new<0) i_new=0;
		j_new=floor(log10(a_pb2[j]/0.1)/size_step);
		//if(log10(a_pb2[j]/0.1)/size_step-j_new*1.0 >=0.5) j_new+=1;
		if(j_new==j) j_new=j+1;
		if(check==1) j_new=j;
		if(j_new <0){
			j_new=0;
//			j_new=j_new+1;
		}
		else if(j_new >peb_size_num-1) {
			j_new=peb_size_num-1;
//			printf("maximum size %e\t%e\t%e\n",2.25*mean_path(peb_map[i].rad),vr0,tau);
		}
		d_size=peb_map[i].size[j_new]-peb_map[i].size[j_new-1];
		delta_size=a_pb2[j]-a_pb1;
		vol2=pow(peb_map[i_new].size[j_new]/peb_map[i].size[j],3)*peb_map[i].mass_out[j];
		if(j_new==j){
		vol1=pow(peb_map[i_new].size[j_new]/peb_map[i].size[j],3)*peb_map[i].mass_out[j];
		}
		else{
		vol1=pow(peb_map[i_new].size[j_new-1]/peb_map[i].size[j],3)*peb_map[i].mass_out[j];
		}
		if(j_new==j) frac_s=0.0;
		else{
		frac_s=(peb_map[i].mass_out[j]+delta_mass[j]-vol1)/(vol2-vol1);
		}
		frac=delta_r/dr-1.0*floor(delta_r/dr);
		if(delta_r/dr > 0.09) printf("AAATOOLARGE!\n");
                frac_s=delta_size/d_size-1.0*floor(delta_size/d_size);
		if(delta_size/d_size > 0.2) printf("%f\tAAAATOOOLARGE!\n",delta_size/d_size);
		if(check==1) frac_s=0.0;
//		if(num_step%100==0 && (i>25 && i < 30) && j==peb_size_num-1){
//			printf("i=%d i_new=%d j_new=%d frac=%f frac_s=%f\n",i,i_new,j_new,frac,peb_map[i].mass_out[j_new]);
//		}
		ratio_size=peb_map[i_new].size[j_new]/peb_map[i].size[j];
		peb_map[i_new].mass_in[j_new]+=pow(ratio_size,3)*frac_s*frac*peb_map[i].mass_out[j];
		ratio_size=peb_map[i_new].size[j_new-1]/peb_map[i].size[j];
		if(j_new==j) {
		ratio_size=1.0;
		peb_map[i_new].mass_in[j_new]+=pow(ratio_size,3)*(1.0-frac_s)*frac*peb_map[i].mass_out[j];
		}
		else {
		peb_map[i_new].mass_in[j_new-1]+=pow(ratio_size,3)*(1.0-frac_s)*frac*peb_map[i].mass_out[j];
		}
	//	if(delta_r/dr>1.0 || 1){
		if(j_new >peb_size_num-1) j_new=peb_size_num-1;
		ratio_size=peb_map[i_new+1].size[j_new]/peb_map[i].size[j];
		peb_map[i_new+1].mass_in[j_new]+=pow(ratio_size,3)*frac_s*(1.0-frac)*peb_map[i].mass_out[j];
		if(j_new==j) {
		ratio_size=1.0;
		peb_map[i_new+1].mass_in[j_new]+=pow(ratio_size,3)*(1.0-frac_s)*(1.0-frac)*peb_map[i].mass_out[j];
		}
		else{
		ratio_size=peb_map[i_new+1].size[j_new-1]/peb_map[i].size[j];
        	peb_map[i_new+1].mass_in[j_new-1]+=pow(ratio_size,3)*(1.0-frac_s)*(1.0-frac)*peb_map[i].mass_out[j];
		}
		peb_map[i].mass_out[j]-=peb_map[i].mass_out[j];
	//	}
	//	else {
	//	ratio_size=peb_map[i].size[j_new]/peb_map[i].size[j];
	//	peb_map[i].mass_in[j_new]+=pow(ratio_size,3)*(1.0-frac)*peb_map[i].mass_out[j];
          //      peb_map[i].mass_out[j]-=frac*peb_map[i].mass_out[j];
	//	}
//	if(num_step%100==0 && (i>25 && i < 30) && j==peb_size_num-1){
//	printf("i=%d i_new=%d j_new=%d frac=%f frac_s=%f\n",i,i_new,j_new,frac,peb_map[i].mass_out[j_new]);
//	}
	}
	}

}

void grow_3(){

	int i,j,i_new,j_new;
	double a_pb1,a_pb11,a_pb2,a_pb22,a_pb3,vr0,vr1,vr2,AREA,dr,a_max;
	double tau,vol_plus,frac,frac_s,coag_eff,ratio_size;
	dr=size_ring;
	coag_eff=1.0;
for(i=ring_num-1;i>-1;i--){
AREA=M_PI*((peb_map[i].rad+dr/2.0)*(peb_map[i].rad+dr/2.0)-(peb_map[i].rad-dr/2.0)*(peb_map[i].rad-dr/2.0))*LUNIT*LUNIT;
a_max=2.25*mean_path(peb_map[i].rad+dr/2.0);
for(j=0;j<peb_size_num;j++){
	a_pb1=peb_map[i].size[j];
	vr0=vr_estimate(peb_map[i].rad+dr/2.0,a_pb1,pp_vr_tau);
	tau=pp_vr_tau[1];
	if(a_pb1>a_max) vol_plus=0.0;
	else{
	vol_plus=1.0*M_PI*a_pb1*a_pb1*sqrt(vr0*vr0+0.25*tau*vr0*tau*vr0)*dt*TUNIT;
	}
	a_pb11=pow(((vol_plus*coag_eff*dust_budget[i].rho[j]/rho_peb0+4.0/3.0*M_PI*a_pb1*a_pb1*a_pb1)*3.0/4.0/M_PI),1.0/3.0);
	if(a_pb1<a_max && a_pb11>a_max) a_pb11=a_max;
	a_pb2=peb_map[i].size[j+1];
        vr0=vr_estimate(peb_map[i].rad+dr/2.0,a_pb2,pp_vr_tau);
	tau=pp_vr_tau[1];
	if(a_pb2>a_max) vol_plus=0.0;
	else{
		vol_plus=1.0*M_PI*a_pb2*a_pb2*sqrt(vr0*vr0+0.25*tau*vr0*tau*vr0)*dt*TUNIT;
	}
	a_pb22=pow(((vol_plus*coag_eff*dust_budget[i].rho[j]/rho_peb0+4.0/3.0*M_PI*a_pb2*a_pb2*a_pb2)*3.0/4.0/M_PI),1.0/3.0);
	if(a_pb2<a_max && a_pb22>a_max) a_pb22=a_max;

	a_pb3=(a_pb2-a_pb1)/(a_pb22-a_pb11)*(a_pb2-a_pb11)+a_pb1;
	frac_s=(a_pb2-a_pb3)/(a_pb2-a_pb1);
//	frac_s=(a_pb22-a_pb2)/(a_pb2-a_pb1);
	if(a_pb2>a_max) frac_s=0.0;
	if(frac_s>0.2) printf("%f\t AAAATOOOLLLLARGE\n",frac_s);
	vr1=vr_estimate(peb_map[i].rad+dr,(a_pb1+a_pb2)/2.0,pp_vr_tau);
	vr2=vr_estimate(peb_map[i].rad,(a_pb1+a_pb2)/2.0,pp_vr_tau);
//	x=dr*(dr-vr2*dt)/(dr+vr1*dt-vr2*dt);
//	frac=(dr-x)/dr;
	frac=1.0-(dr-vr2*dt*TUNIT/LUNIT)/(dr+vr1*dt*TUNIT/LUNIT-vr2*dt*TUNIT/LUNIT);
	//frac=vr2*dt*TUNIT/LUNIT/dr;
	j_new=j+1;
	if(a_pb2>a_max || j_new>peb_size_num-1) {j_new=j; frac_s=0.0;}
	i_new=i-1;
	if(i_new<0) i_new=0;
	peb_map[i_new].mass_in[j_new]+=pow(a_pb22/a_pb2,3)*frac*frac_s*peb_map[i].mass_out[j];
	peb_map[i].mass_in[j_new]+=pow(a_pb22/a_pb2,3)*(1.0-frac)*frac_s*peb_map[i].mass_out[j];
	peb_map[i_new].mass_in[j]+=pow(a_pb2/a_pb3,3)*frac*(1.0-frac_s)*peb_map[i].mass_out[j];
	peb_map[i].mass_in[j]+=pow(a_pb2/a_pb3,3)*(1.0-frac)*(1.0-frac_s)*peb_map[i].mass_out[j];
	peb_map[i].mass_out[j]-=peb_map[i].mass_out[j];

/*	ratio_size=peb_map[i_new].size[j_new]/peb_map[i].size[j];
		peb_map[i_new].mass_in[j_new]+=pow(ratio_size,3)*frac_s*frac*peb_map[i].mass_out[j];
		ratio_size=peb_map[i_new].size[j_new-1]/peb_map[i].size[j];
		if(j_new==j) {
		ratio_size=1.0;
		peb_map[i_new].mass_in[j_new]+=pow(ratio_size,3)*(1.0-frac_s)*frac*peb_map[i].mass_out[j];
		}
		else {
		peb_map[i_new].mass_in[j_new-1]+=pow(ratio_size,3)*(1.0-frac_s)*frac*peb_map[i].mass_out[j];
		}
	//	if(delta_r/dr>1.0 || 1){
		if(j_new >peb_size_num-1) j_new=peb_size_num-1;
		ratio_size=peb_map[i_new+1].size[j_new]/peb_map[i].size[j];
		peb_map[i_new+1].mass_in[j_new]+=pow(ratio_size,3)*frac_s*(1.0-frac)*peb_map[i].mass_out[j];
		if(j_new==j) {
		ratio_size=1.0;
		peb_map[i_new+1].mass_in[j_new]+=pow(ratio_size,3)*(1.0-frac_s)*(1.0-frac)*peb_map[i].mass_out[j];
		}
		else{
		ratio_size=peb_map[i_new+1].size[j_new-1]/peb_map[i].size[j];
        	peb_map[i_new+1].mass_in[j_new-1]+=pow(ratio_size,3)*(1.0-frac_s)*(1.0-frac)*peb_map[i].mass_out[j];
		}
		peb_map[i].mass_out[j]-=peb_map[i].mass_out[j];
*/			        
}

}

}

void grow_3b(){

	int i,j,i_new,j_new;
	double a_pb1,a_pb11,a_pb2,a_pb22,a_pb3,vr0,vr1,vr2,AREA,dr,a_max;
	double tau,vol_plus,frac,frac_s,coag_eff,ratio_size,ring_mass_before,ring_mass_after,old_sigma;
	dr=size_ring;
	coag_eff=1.0;
for(i=ring_num-1;i>-1;i--){
AREA=M_PI*((peb_map[i].rad+dr/2.0)*(peb_map[i].rad+dr/2.0)-(peb_map[i].rad-dr/2.0)*(peb_map[i].rad-dr/2.0))*LUNIT*LUNIT;
a_max=2.25*mean_path(peb_map[i].rad+dr/2.0);
ring_mass_before=0.0;
ring_mass_after=0.0;
for(j=0;j<peb_size_num;j++){
	ring_mass_before+=peb_map[i].mass_out[j];
}
for(j=0;j<peb_size_num;j++){
	a_pb1=peb_map[i].size[j];
	vr0=vr_estimate(peb_map[i].rad+dr/2.0,a_pb1,pp_vr_tau);
	tau=pp_vr_tau[1];
	if(a_pb1>a_max) vol_plus=0.0;
	else{
	vol_plus=1.0*M_PI*a_pb1*a_pb1*sqrt(vr0*vr0+0.25*tau*vr0*tau*vr0)*dt*TUNIT;
	}
	a_pb11=pow(((vol_plus*coag_eff*dust_budget[i].rho[j]/rho_peb0+4.0/3.0*M_PI*a_pb1*a_pb1*a_pb1)*3.0/4.0/M_PI),1.0/3.0);
	if(a_pb1<a_max && a_pb11>a_max) a_pb11=a_max;
	a_pb2=peb_map[i].size[j+1];
        vr0=vr_estimate(peb_map[i].rad+dr/2.0,a_pb2,pp_vr_tau);
	tau=pp_vr_tau[1];
	if(a_pb2>a_max) vol_plus=0.0;
	else{
		vol_plus=1.0*M_PI*a_pb2*a_pb2*sqrt(vr0*vr0+0.25*tau*vr0*tau*vr0)*dt*TUNIT;
	}
	a_pb22=pow(((vol_plus*coag_eff*dust_budget[i].rho[j]/rho_peb0+4.0/3.0*M_PI*a_pb2*a_pb2*a_pb2)*3.0/4.0/M_PI),1.0/3.0);
	if(a_pb2<a_max && a_pb22>a_max) a_pb22=a_max;

	a_pb3=(a_pb2-a_pb1)/(a_pb22-a_pb11)*(a_pb2-a_pb11)+a_pb1;
	frac_s=(a_pb2-a_pb3)/(a_pb2-a_pb1);
//	frac_s=(a_pb22-a_pb2)/(a_pb2-a_pb1);
	if(a_pb2>a_max) frac_s=0.0;
	if(frac_s>0.2) printf("%f\t AAAATOOOLLLLARGE\n",frac_s);
	vr1=vr_estimate(peb_map[i].rad+dr,(a_pb1+a_pb2)/2.0,pp_vr_tau);
	vr2=vr_estimate(peb_map[i].rad,(a_pb1+a_pb2)/2.0,pp_vr_tau);
//	x=dr*(dr-vr2*dt)/(dr+vr1*dt-vr2*dt);
//	frac=(dr-x)/dr;
	frac=1.0-(dr-vr2*dt*TUNIT/LUNIT)/(dr+vr1*dt*TUNIT/LUNIT-vr2*dt*TUNIT/LUNIT);
	//frac=vr2*dt*TUNIT/LUNIT/dr;
	j_new=j+1;
	if(a_pb2>a_max || j_new>peb_size_num-1) {j_new=j; frac_s=0.0;}
	i_new=i-1;
	if(i_new<0) i_new=0;
	peb_map[i_new].mass_in[j_new]+=pow(a_pb22/a_pb2,3)*frac*frac_s*peb_map[i].mass_out[j];
	ring_mass_after+=pow(a_pb22/a_pb2,3)*frac*frac_s*peb_map[i].mass_out[j];
	peb_map[i].mass_in[j_new]+=pow(a_pb22/a_pb2,3)*(1.0-frac)*frac_s*peb_map[i].mass_out[j];
	ring_mass_after+=pow(a_pb22/a_pb2,3)*(1.0-frac)*frac_s*peb_map[i].mass_out[j];
	peb_map[i_new].mass_in[j]+=pow(a_pb2/a_pb3,3)*frac*(1.0-frac_s)*peb_map[i].mass_out[j];
	ring_mass_after+=pow(a_pb2/a_pb3,3)*frac*(1.0-frac_s)*peb_map[i].mass_out[j];
	peb_map[i].mass_in[j]+=pow(a_pb2/a_pb3,3)*(1.0-frac)*(1.0-frac_s)*peb_map[i].mass_out[j];
	ring_mass_after+=pow(a_pb2/a_pb3,3)*(1.0-frac)*(1.0-frac_s)*peb_map[i].mass_out[j];
	peb_map[i].mass_out[j]-=peb_map[i].mass_out[j];

/*	ratio_size=peb_map[i_new].size[j_new]/peb_map[i].size[j];
		peb_map[i_new].mass_in[j_new]+=pow(ratio_size,3)*frac_s*frac*peb_map[i].mass_out[j];
		ratio_size=peb_map[i_new].size[j_new-1]/peb_map[i].size[j];
		if(j_new==j) {
		ratio_size=1.0;
		peb_map[i_new].mass_in[j_new]+=pow(ratio_size,3)*(1.0-frac_s)*frac*peb_map[i].mass_out[j];
		}
		else {
		peb_map[i_new].mass_in[j_new-1]+=pow(ratio_size,3)*(1.0-frac_s)*frac*peb_map[i].mass_out[j];
		}
	//	if(delta_r/dr>1.0 || 1){
		if(j_new >peb_size_num-1) j_new=peb_size_num-1;
		ratio_size=peb_map[i_new+1].size[j_new]/peb_map[i].size[j];
		peb_map[i_new+1].mass_in[j_new]+=pow(ratio_size,3)*frac_s*(1.0-frac)*peb_map[i].mass_out[j];
		if(j_new==j) {
		ratio_size=1.0;
		peb_map[i_new+1].mass_in[j_new]+=pow(ratio_size,3)*(1.0-frac_s)*(1.0-frac)*peb_map[i].mass_out[j];
		}
		else{
		ratio_size=peb_map[i_new+1].size[j_new-1]/peb_map[i].size[j];
        	peb_map[i_new+1].mass_in[j_new-1]+=pow(ratio_size,3)*(1.0-frac_s)*(1.0-frac)*peb_map[i].mass_out[j];
		}
		peb_map[i].mass_out[j]-=peb_map[i].mass_out[j];
*/			        
}
for(j=0;j<peb_size_num;j++){
	old_sigma=dust_budget[i].surf_dens[j];
	dust_budget[i].surf_dens[j]-=(ring_mass_after-ring_mass_before)/AREA;
	dust_budget[i].rho[j]=dust_budget[i].rho[j]*dust_budget[i].surf_dens[j]/old_sigma;
	if(dust_budget[i].surf_dens[j]<0.0) printf("OMG\n");
}
}

}
