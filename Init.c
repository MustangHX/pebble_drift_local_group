#include "global_var.h"
#include "global_ex.h"
#include "ex_func.h"
#include "math.h"
#include <stdio.h>
double Sigma (double r){
	return 108*pow(r,-0.6);
}
void Init(){
	int i,j,ii;
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
		dust_budget[i].rad=i*size_ring+1.0;
		dust_budget[i].dr=size_ring;
		dust_budget[i].surf_dens=Sigma(dust_budget[i].rad)*0.01;
	}
i=0;
	for(i=0;i<ring_num;i++){
  //              printf("SIGMAAAAA%e\t%d\t",peb_map[i].surf_dens,i);
//		printf("%d\n",i);
		peb_map[i].dr=size_ring;
		peb_map[i].rad=i*size_ring+1.0;
                peb_map[i].time=0.0;
		for(j=0;j<peb_size_num;j++){
			peb_map[i].size[j]=0.1*pow(10,j*size_step);
			if(i==1) printf("SIZE=%fcm\n",peb_map[i].size[j]);
				//drag_group((i+1)*0.25,0.01*(j+1));
			AREA=M_PI*((peb_map[i].rad+0.125)*(peb_map[i].rad+0.125)-(peb_map[i].rad-0.125)*(peb_map[i].rad-0.125))*LUNIT*LUNIT;
			if (j<1 || 1) {
				//peb_map[i].mass_out[j]=0.1*AREA*0.01*(Sigma((i+4)*0.25)*exp(-1.0*peb_map[i].size[j]/0.1)+1e-10);
				peb_map[i].mass_out[j]=10.0*AREA*(dust_budget[i].surf_dens*exp(-1.0*peb_map[i].size[j]/0.1)+1e-10);
			}
			else peb_map[i].mass_out[j]=1e-20*AREA;
			peb_map[i].mass_out[j]/=peb_size_num;
			peb_map[i].surf_dens[j]=peb_map[i].mass_out[j]/AREA;
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
