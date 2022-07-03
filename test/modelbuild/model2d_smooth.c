/* 2-D ENO interpolation of the velocity model */
#include <rsf.h>
#include "velocity.h"

int main(int argc, char* argv[])
{
        sf_file in, out, sz_file;
        float* slow;
	int n[2]={301,1001};
	float o[2]={0.,-2.};
	float d[2]={0.01,0.01};
	int nsv[2];
	float osv[2];
	float dsv[2];
	float *sv;
	float *sz;
	int im;
	int nsz[2];
	float dsz[2];
	float osz[2];
	float vl[3]={1.508,1.69,2.};

        sf_init (argc, argv);
        in  = sf_input ("in");
        out = sf_output("out");
	sz_file = sf_input("sz");

	slow = sf_floatalloc(301*1001);
	for(im=0;im<n[0]*n[1];im++)
		slow[im]=1.5;

        if (!sf_histint(in,"n1",nsv)) sf_error("No n1= in input");    
        if (!sf_histfloat(in,"d1",dsv)) sf_error("No d1= in input");
        if (!sf_histfloat(in,"o1",osv)) sf_error("No o1= in input");
        if (!sf_histint(in,"n2",nsv+1)) sf_error("No n2= in input");    
        if (!sf_histfloat(in,"d2",dsv+1)) sf_error("No d2= in input");
        if (!sf_histfloat(in,"o2",osv+1)) sf_error("No o2= in input");

        if (!sf_histint(sz_file,"n1",nsz)) sf_error("No n1= in input");    
        if (!sf_histfloat(sz_file,"d1",dsz)) sf_error("No d1= in input");
        if (!sf_histfloat(sz_file,"o1",osz)) sf_error("No o1= in input");
        if (!sf_histint(sz_file,"n2",nsz+1)) sf_error("No n2= in input");    
        if (!sf_histfloat(sz_file,"d2",dsz+1)) sf_error("No d2= in input");
        if (!sf_histfloat(sz_file,"o2",osz+1)) sf_error("No o2= in input");

	sv = sf_floatalloc(nsv[0]*nsv[1]);
	sf_floatread(sv,nsv[0]*nsv[1],in);
	sz = sf_floatalloc(nsz[0]*nsz[1]);
	sf_floatread(sz,nsz[0]*nsz[1],sz_file);

	//interpolation
	updateVelocityModel3(slow,n,o,d,sv,nsv,osv,dsv,sz,nsz,osz,dsz,1,vl);

	sf_putint(out,"n1",301);
        sf_putint(out,"n2",1001);
        sf_putfloat(out,"d1",0.01);
        sf_putfloat(out,"d2",0.01);
	sf_putfloat(out,"o1",0.);
	sf_putfloat(out,"o2",-2.);
        sf_floatwrite(slow,301*1001,out);
}

