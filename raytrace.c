/* Ray tracing interface modified to NIP tomography. */
/*
 Copyright (C) 2004 University of Texas at Austin
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <rsf.h>
#include <stdbool.h>
#include "raytrace.h"
#include "grid2.h"
#include "atela.h"
#include "dynamic.h"
#include <math.h>
/*^*/

#ifndef _raytrace_h

#define DEG2RAD SF_PI/180. /* Degrees to radians conversion */
#define ORDER 4 /* Ray tracing interpolation order */
/*^*/

typedef struct RayTrace* raytrace;
/* abstract data type */
/*^*/

#endif

struct RayTrace {
    bool sym;
    int dim, nt;
    float dt, z0;
    grid2 grd2;
};
/* concrete data type */

static void iso_rhs(void* par, float* y, float* f){}

float getVelocityForRaySampleLocation(
		  void *par, /* Raytrace struct */
		  float **traj, /* Ray trajectory */
		  int i /* Sample index in ray trajectory */)
/*< Get velocity from a point in the ray trajectory
Note: The i variable is the sample index of the ray trajectory traj. This
function returns the velocity from the grid in the sample location
(z,x)=(traj[i][0],traj[i][1]).
>*/
{
	raytrace rt;
	float x[2];

	rt = (raytrace) par;

	x[0]=traj[i][0];
	x[1]=traj[i][1];
	
	return sqrtf(1./grid2_vel(rt->grd2,x));
}

static void dyn_iso_rhs(void* par, float dvdn, float *x, float* y, float* f)
/* right-hand side for isotropic raytracing */
{    
	raytrace rt;
	float v;
	
	rt = (raytrace) par;

	v = sqrtf(1./grid2_vel(rt->grd2,x));
	
	f[0]   = v*v*y[1]; /* v^2 p */
	f[1] = (-1./v)*dvdn*y[0]; /* -1/v dv/dn q */
}

void first_deriv( float h /* step */,
		float *zx /* f(x) */,
		float *der /* Derivative */,
		int n /* vectors dim */)
/*< Calculate first derivative numerically >*/
{
	int i;

	sf_deriv_init(n, 6, 0.);
        sf_deriv(zx,der);

        for (i=0; i < n; i++) {
                der[i] /= h;
        }
}


void second_deriv( float h /* passo h */,
		    float *zx /* f(x) */,
		    float *der /* Derivative */,
		    int n /* vectors dim */)
/*< Calculate second derivative numerically >*/
{
	float *firstder;
	firstder = sf_floatalloc(n);
	first_deriv(h,zx,firstder,n);
	first_deriv(h,firstder,der,n);
	free(firstder);
}


float second_derivative(void *par, float *n, float *x, float v)
/*< Calculate the second derivative in the direction normal to ray trajectory >*/
{
	raytrace rt;
	float vpdx, vmdx;
	float dx=0.01;
	float tmp[2];
	float der[3];
	float fx[3];

	rt = (raytrace) par;

	n[0]*=dx;
	n[1]*=dx;

	// ppdx
	tmp[0] = x[0]+n[1];
	tmp[1] = x[1]-n[0];
	fx[2] = sqrtf(1./grid2_vel(rt->grd2,tmp));

	fx[1]=v;

	// pmdx
	tmp[0]=x[0]-n[1];
	tmp[1]=x[1]+n[0];
	fx[0] = sqrtf(1./grid2_vel(rt->grd2,tmp));	

	second_deriv(0.01,fx,der,3);
	//if(der[1]>0.001 || der[1] <-0.001) sf_warning("v=%f %f %f d=%f %f",fx[0],fx[1],fx[2],der[1],(fx[0]-2*fx[1]+fx[2])/(dx*dx));
	//sf_warning("ddnv=%f d=%f",der[1],(fx[0]-2*fx[1]+fx[2])/(dx*dx));
	//der[1]=(fx[0]-2*fx[1]+fx[2])/(dx*dx);
	return der[1];
}

float calculateRNIPWithDynamicRayTracing(
					  void *par, /* Raytrace struct */
					  float dt, /* Time sampling */
					  float nt, /* Number of time samples */
					  float **traj, /* Ray trajectory (z,x) */
					  float v0 /* Near surface velocity */
)
/*< Calculate RNIP with dynamic ray tracing >*/
{
	float v; // velocity
	int it; // Loop counter
	raytrace rt; // Raytrace struct
	float *x; // Sample coordinate (z,x)
	float *n; // Normal vector
	float *dvdn; // Derivative normal to ray direction
	float mod; // tmp variable
	float rnip; // RNIP parameter

    	rt = (raytrace) par;
	x = sf_floatalloc(2);
	n = sf_floatalloc(2);
	dvdn = sf_floatalloc(nt);

	for(it=0;it<nt;it++){
		x[0]=traj[it][0];
		x[1]=traj[it][1];
		n[0] = (traj[it+1][0]-traj[it][0]);
		n[1] = (traj[it+1][1]-traj[it][1]);
		mod = sqrtf(n[0]*n[0]+n[1]*n[1]);
		n[0] /= mod;
		n[1] /= mod;

		v = sqrtf(1./grid2_vel(rt->grd2,x));

		//sf_warning("v=%f",v);
		/* Calculate derivative for each ray sample */
		dvdn[it]=second_derivative(rt,n,x,v);

	} // Loop over ray samples

	//sf_error("oi");
	/* Initial conditions for a point source */
	x[0]=0.; // q=0
	x[1]=1.; // p=1

	/* Fourth order Runge-Kutta dynamic ray tracing */
	sf_dynamic_runge_init(2,nt,dt);
	rnip = sf_dynamic_runge_step(x,rt,dyn_iso_rhs,dvdn,traj,v0);
	sf_dynamic_runge_close();

	/*if(rnip<0.){
		for(it=0;it<nt;it++)
			if(dvdn[it]<0.001 || dvdn[it] > 0.001)
				sf_warning("rnip=%f dvdn[%d]=%f",rnip,it,dvdn[it]);
	}*/
	return rnip;
}


void sortingXinAscendingOrder(
				float *x, /* x vector to sort */
				float *z, /* z(x) vector */
				int n /* Vectors dimension */)
/*< x vector sorting in ascending order >*/
{
	int i; // Loop counter
	float tmpx, tmpz; // Temporary variables
	int k; // Sorting key (number of changes)

	do{
		k=0;
		for(i=1;i<n;i++){
			if(x[i-1]>x[i]){
				tmpx=x[i-1]; tmpz=z[i-1];
				x[i-1]=x[i]; z[i-1]=z[i];
				x[i]=tmpx; z[i]=tmpz;
				k++;
			}
		} // Loop vector samples
	}while(k!=0);
}

int binarySearch(float xx, float *x, int n)
/*< Binary search to get spline index >*/
{
	int ini=0;
	int fin=n-1;
	int meio=0;
	int l=-1;

	do{
		meio=(ini+fin)/2;
		if(xx>x[meio]){
			ini=meio;
			if(xx<x[meio+1])
				l=meio;
		}else{
			fin=meio;
			if(xx>x[meio-1])
				l=meio-1;
		}
	}while(l==-1);
	return l;
}

static int term(void* par, float* y)
/* grid termination */
{
    raytrace rt;
    
    rt = (raytrace) par;
	
    switch (rt->dim) {
		case 2:
			return grid2_term(rt->grd2,y);
		default:
			sf_error("%s: Cannot raytrace with dim=%d",__FILE__,rt->dim);
			return 0;
    }
}

raytrace raytrace_init(int dim            /* dimensionality (2 or 3) */, 
					   bool sym,          /* if symplectic */
					   int nt             /* number of ray tracing steps */, 
					   float dt           /* ray tracing step (in time) */,
					   int* n             /* slowness dimensions [dim] */, 
					   float* o, float* d /* slowness grid [dim] */,
					   float* slow2       /* slowness squared [n3*n2*n1] */, 
					   int order          /* interpolation order */)
/*< Initialize ray tracing object. 
 * Increasing order increases accuracy but
 decreases efficiency. Recommended values: 3 or 4.
 * slow2 can be changed or deallocated after
 raytrace_init.
 >*/
{
    raytrace rt;
    
    rt = (raytrace) sf_alloc (1,sizeof(*rt));
    
    rt->dim = dim;
    rt->sym = sym;
    rt->nt = nt;
    rt->dt = dt;
    rt->z0 = o[0];
    
    switch (dim) {
		case 2:
			rt->grd2 = grid2_init (n[0], o[0], d[0], 
								   n[1], o[1], d[1],
								   slow2, order);
			break;
		default:
			sf_error("%s: Cannot raytrace with dim=%d",__FILE__,dim);
    }
	
    return rt;
}

void raytrace_close (raytrace rt)
/*< Free internal storage >*/
{
    switch (rt->dim) {
		case 2:
			grid2_close (rt->grd2);
			break;
    }
    free (rt);
}

int trace_ray (raytrace rt  /* ray tracing object */, 
			   float* x     /* point location {z,y,x} [dim] */, 
			   float* p     /* ray parameter vector [dim] */, 
			   float** traj /* output ray trajectory [nt+1,dim] */)
/*< Trace a ray.
 * Values of x and p are changed inside the function.
 * The trajectory traj is stored as follows:
 {z0,y0,z1,y1,z2,y2,...} in 2-D
 {z0,y0,x0,z1,y1,x1,...} in 3-D
 * Vector p points in the direction of the ray. 
 The length of the vector is not important.
 Example initialization:
 p[0] = cos(a); p[1] = sin(a) in 2-D, a is between 0 and 2*pi radians
 p[0] = cos(b); p[1] = sin(b)*cos(a); p[2] = sin(b)*sin(a) in 3-D
 b is inclination between 0 and   pi radians
 a is azimuth     between 0 and 2*pi radians
 * The output code for it = trace_ray(...)
 it=0 - ray traced to the end without leaving the grid
 it>0 - ray exited at the top of the grid
 it<0 - ray exited at the side or bottom of the grid
 * The total traveltime along the ray is 
 nt*dt if (it = 0); abs(it)*dt otherwise 
 >*/
{
    int i, dim, it=0, nt;
    float y[6], s2;
	
    dim = rt->dim;
    nt = rt->nt;
	
    if (!rt->sym) {
		switch (dim) {
			case 2:
				s2 = grid2_vel(rt->grd2,x);
				break;
			default:
				s2 = 0.;
				sf_error("%s: Cannot raytrace with dim=%d",__FILE__,dim);
		}
		
		for (i=0; i < dim; i++) {
			y[i] = x[i];
			y[i+dim] = p[i]*sqrtf(s2);
		}
		
		sf_runge_init(2*dim, nt, rt->dt);
		it = sf_ode23_step (y, rt,iso_rhs,term,traj);
		sf_runge_close();
		
		for (i=0; i < dim; i++) {
			x[i] = y[i];
			p[i] = y[i+dim];
		}
    } else {
		switch (dim) {
			case 2:
				it = atela_step (dim, nt, rt->dt, true, x, p, 
								 rt->grd2, 
								 grid2_vgrad, grid2_vel, grid2_term, traj);
				break;
			default:
				sf_error("%s: cannot handle %d dimensions",__FILE__,rt->dim);
				break;
		}
    }
    
    if (it > 0 && x[0] > rt->z0) {
		return (-it); /* exit through the side or bottom */
    } else {
		return it;
    }
}

/* 	$Id$	 */
