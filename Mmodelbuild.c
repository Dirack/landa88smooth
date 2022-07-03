/* Model build for a constant velocity layers with smooth interfaces

This program is used after NIP tomography to build velocity model from inversion
*/

#include <rsf.h>

int main(int argc, char* argv[]){

	float *vel; // Velocity model
	float **zi; // Interfaces nodes
	float *sv; // Layer's velocities
	int nsv; // sv vector dimension
	float dsv; // sv vector sampling
	float osv; // sv vector origin
	int nz; // Velocity model z dimension
	float oz, dz; // Velocity model origin and sampling
	int nsz1, nsz2; // Interfaces dimension
	float osz1, osz2, dsz1, dsz2; // Interfaces origin and sampling
	float z; // Depth
	float xx; // Distance
	int i, j; // Loop counter
	int itf; // Interfaces index
	sf_file vspline; // Layer's velocities
	sf_file zspline; // Interfaces nodes
	sf_file model; // velocity model

	sf_init(argc,argv);

	vspline = sf_input("in");
	zspline = sf_input("sz");
	model = sf_output("out");

	if(!sf_histint(vspline,"n1",&nsv)) sf_error("No n1= in input file");
	if(!sf_histfloat(vspline,"d1",&dsv)) sf_error("No d1= in input file");
	if(!sf_histfloat(vspline,"o1",&osv)) sf_error("No o1= in input file");

	if(!sf_histint(zspline,"n1",&nsz1)) sf_error("No n1= in sz file");
	if(!sf_histfloat(zspline,"d1",&dsz1)) sf_error("No d1= in sz file");
	if(!sf_histfloat(zspline,"o1",&osz1)) sf_error("No o1= in sz file");
	if(!sf_histint(zspline,"n2",&nsz2)) sf_error("No n2= in sz file");
	if(!sf_histfloat(zspline,"d2",&dsz2)) sf_error("No d2= in sz file");
	if(!sf_histfloat(zspline,"o2",&osz2)) sf_error("No o2= in sz file");

	if(!sf_getint("nz",&nz)) nz=1;
	if(!sf_getfloat("oz",&oz)) oz=0.;
	if(!sf_getfloat("dz",&dz)) dz=0.01;

	vel = sf_floatalloc(nz*nsz1);
	sv = sf_floatalloc(nsv);
	sf_floatread(sv,nsv,vspline);
	zi = sf_floatalloc2(nsz1,nsz2);
	sf_floatread(zi[0],nsz1*nsz2,zspline);

	/* Calculate velocity function */
	for(j=0;j<nsz1;j++){
		xx = j*dsz1+osz1;
		itf=0;
		/* Calculate interfaces z coordinates */
		for(i=0;i<nz;i++){
			z = i*dz+oz;
			if(z>zi[itf][j] && itf < nsz2-1) itf++;
			if(z>zi[nsz2-1][j]){
				vel[nz*j+i] = sv[nsv-1];
			}else{
				vel[nz*j+i] = sv[itf];
			}
		} /* Loop over depth */
	} /* Loop over distance */


	sf_putint(model,"n1",nz);
	sf_putint(model,"n2",nsz1);
	sf_putfloat(model,"d1",dz);
	sf_putfloat(model,"d2",dsz1);
	sf_putfloat(model,"o1",oz);
	sf_putfloat(model,"o2",osz1);
	sf_floatwrite(vel,nsz1*nz,model);
}
