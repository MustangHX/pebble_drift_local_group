#include "global_var.h"
#include "global_ex.h"
#include "ex_func.h"
#include "math.h"
#include <stdio.h>
double Sigma (double r){
	return 108*pow(r,-0.6);
}
void Init(){
	int i,j;
	double AREA;
/*        for(i=0;i<100;i++){
                size=drag_group(i*1.0+0.5,0.01);
                for(j=0;j<100;j++){		
			if(j==(int)(size*10)){
                               pebble[i][j]=i*1.0*2*M_PI*dr;
		                             }
		                  }
		          }
	fp=fopen("velocity.dat","r");
	if(fp=None){
		for(j=0;j<100;j++){
			size=j*1.0;
			for(i=0;i<100;i++){
				velocity[i][j]
			
}
}
}*/
	
	for(i=0;i<ring_num;i++){
		dust_budget[i].rad=i*size_ring+1.0-size_ring;
		dust_budget[i].dr=size_ring;
		for(j=0;j<peb_size_num;j++){
                dust_budget[i].surf_dens[j]=Sigma(dust_budget[i].rad)*0.1;
		dust_budget[i].rho[j]=density(dust_budget[i].rad)*0.1;
		}
	}
i=0;
	for(i=0;i<ring_num;i++){
  //              printf("SIGMAAAAA%e\t%d\t",peb_map[i].surf_dens,i);
//		printf("%d\n",i);
		peb_map[i].dr=size_ring;
		peb_map[i].rad=i*size_ring+1.0-size_ring;
                peb_map[i].time=0.0;
		for(j=0;j<=peb_size_num;j++){
			peb_map[i].size[j]=0.1*pow(10,j*size_step);
		}
		for(j=0;j<peb_size_num;j++){
			peb_map[i].size_med[j]=0.5*(peb_map[i].size[j]+peb_map[i].size[j+1]);
			if(i==1) printf("SIZE=%fcm\n",peb_map[i].size[j]);
				//drag_group((i+1)*0.25,0.01*(j+1));
			AREA=M_PI*((peb_map[i].rad+size_ring/2.0)*(peb_map[i].rad+size_ring/2.0)-(peb_map[i].rad-size_ring/2.0)*(peb_map[i].rad-size_ring/2.0))*LUNIT*LUNIT;
			if (j<1 && i==ring_num-1|| 1) {
				//peb_map[i].mass_out[j]=0.1*AREA*0.01*(Sigma((i+4)*0.25)*exp(-1.0*peb_map[i].size[j]/0.1)+1e-10);
				peb_map[i].mass_out[j]=1.0*AREA*(dust_budget[i].surf_dens[0]*exp(-1.0*peb_map[i].size_med[j]/0.1)+1e-10);
			}
			else peb_map[i].mass_out[j]=1e-20*AREA;
			peb_map[i].mass_out[j]/=peb_size_num;//comes from dust_budget
			peb_map[i].surf_dens[j]=peb_map[i].mass_out[j]/AREA;
			peb_map[i].mass_in[j]=0.0;
		}
  //              printf("SIGMAAAAA%e\t%d\n",peb_map[i].surf_dens,i);
	}
	for(i=0;i<ring_num;i++){
	for(j=0;j<peb_size_num;j++){
//		printf("%d%d %e\t",i,j,peb_map[i].surf_dens[j]);
	}
//	printf("\n");
	}
}

void Restart(int rnum){
	double AREA,dens;
	int i,j,k;
	FILE *fp;
	char name[256];
	printf("RESTART=%d\n",rnum);
        sprintf(name,"out_sigma%d.txt",rnum);
	printf("%s\n",name);
	fp=fopen(name,"r");
        for(i=0;i<ring_num;i++){
	AREA=M_PI*((peb_map[i].rad+size_ring/2.0)*(peb_map[i].rad+size_ring/2.0)-(peb_map[i].rad-size_ring/2.0)*(peb_map[i].rad-size_ring/2.0))*LUNIT*LUNIT;

        for(j=0;j<peb_size_num;j++){
		fscanf(fp,"%lf",&dens);
		peb_map[i].surf_dens[j]=dens;
		peb_map[i].mass_out[j]=AREA*peb_map[i].surf_dens[j];
		peb_map[i].mass_in[j]=0.0;
	}
	}
}
