/*
	 velocity_lib.c (c)
	 
	 Purpose: Functions to update velocity model.
	 	 
	 Version 1.0
	 
	 Site: https://dirack.github.io
	 
	 Programmer: Rodolfo A. C. Neves (Dirack) 19/09/2021

	 Email:  rodolfo_profissional@hotmail.com

	 License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

*/

#include <stdio.h>
#include <stdlib.h>
#include <rsf.h>
#include "velocity_lib.h"
#include "raytrace.h"

void calcInterfacesZcoord(      float *zi, /* Interfaces depth coordinates */
                                int nint, /* Number of interfaces */
                                float xs, /* x coordinate */
                                int si, /* Spline index */
                                float **coef /* Cubic spline coefficients */)
/*< Calculate depth coordinates of the interfaces
 * Note: This function calculates interfaces depth coordinates and stores it
 * in the zi vector.
  >*/
{
        int i; // Loop counter

        for(i=0;i<nint;i++)
                zi[i] = coef[i][si*4+0]*xs*xs*xs+coef[i][si*4+1]*xs*xs+coef[i][si*4+2]*xs+coef[i][si*4+3];
}

//TODO: correct "coefficients" in the function name
void calculateSplineCoeficients(int n, /* Vectors (x,y) dimension */
                                float* x, /* x coordinates */
                                float* y, /* y coordinates */
                                float* coef /* Spline coefficients */)
/*< Function to calculate natural cubic spline coefficients

Note: It Receives n points and two vectors x and y with n dimension.
It returns a coefficients vector with 4 coefficients for each of the
n-1 natural cubic splines, coef[(n-1)*4].

IMPORTANT: The number of points must be equal or major than 3 (n>3)
and x vector must be in crescent order.

>*/
{

        float s2[n]; // Second derivatives matrix
        int i, ip1, ip2, im1, m; // Loop counter
        float hb, ha, deltaa, deltab, t; // temporary variables
        float e[n-2]; // hi's vector
        float dp[n-2]; // main diagonal

        /* Vectors dimension must be major than 3 */
        if(n<3) sf_error("Vectors dimension n must be major than 3\n");

        /* x vector must be in crescent order */
        for(i=1;i<n;i++){
                if(x[i-1]>x[i]) sf_error("Vector x should be in ascending order\n");
        }

        /* Simetric tridiagonal linear system build */
        ha = x[1]-x[0]; deltaa = (y[1]-y[0])/ha; m=n-2;
        for(i=0;i<m;i++){
                ip1 = i+1; ip2 = i+2;
                hb = x[ip2]-x[ip1];
                deltab = (y[ip2]-y[ip1])/hb;
                e[i] = hb; dp[i] = 2*(ha+hb);
                s2[ip1] = 6*(deltab-deltaa);
                ha=hb; deltaa=deltab;
        }

        /* Gauss elimination */
        for(i=1;i<m;i++){
                ip1=i+1; im1=i-1;
                t = e[im1]/dp[im1];
                dp[i] = dp[i]-t*e[im1];
                s2[ip1] = s2[ip1]-t*s2[i];
        }

        /* Retroactive substitutive solution */
        s2[m]=s2[m]/dp[m-1];
        for(i=m-1;i>0;i--){
                ip1=i+1; im1=i-1;
                s2[i]=(s2[i]-e[im1]*s2[ip1])/dp[im1];
        }
        s2[0]=0; s2[n-1]=0;
        /* Calculate spline coefficients */
        for(i=0;i<n-1;i++){
                ha = x[i+1]-x[i];
                coef[0+i*4] = (s2[i+1]-s2[i])/(6*ha);
                coef[1+i*4] = s2[i]/2;
                coef[2+i*4] = (y[i+1]-y[i])/ha-(s2[i+1]+2*s2[i])*(ha/6);
                coef[3+i*4] = y[i];
        }
}

void interfaceInterpolationFromNipSources(float **s, /* NIP sources */
                                          int ns, /* Number of NIP sources */
                                          float *sz, /* Spline nodepoints */
                                          int *nsz, /* Number of nodepoints */
                                          float *osz, /* Nodepoints origin */
                                          float *dsz, /* Nodepoints sampling */
					  int itf /* Interface index */)
/*< Use NIP sources location to draw interfaces 
Note: If the velocity model is correct the NIP sources location coincides with interfaces. So, they can be used to draw interface through cubic spline interpolation.
>*/
{
        int nxs=ns; // Number of NIP sources for each interface
        int nxsz=nsz[0]; // Number of nodepoints for each interface
        int i, im; // Loop counter
        float *tsx, *tsz; // Temporary spline vector
        float *coef; // Coefficients matrix
        float xx, xs; // X coordinate
        float oxs; // Spline's origin
        int l; // Spline index
	float *szz; // Interface tmp vector

        tsx = sf_floatalloc(nxs+2);
        tsz = sf_floatalloc(nxs+2);
        coef = sf_floatalloc(4*(nxs+2-1));
	szz = sf_floatalloc(nsz[0]);
        for(i=0;i<1;i++){
                for(im=0;im<nxs;im++){
                        tsz[im]=s[im][0];
                        tsx[im]=s[im][1];
                }
                sortingXinAscendingOrder(tsx,tsz,nxs);
                calculateSplineCoeficients(nxs,tsx,tsz,coef);
                oxs=tsx[0];
                for(im=0;im<nxsz;im++){
                        xx=im*dsz[0]+osz[0];
                        if(xx<tsx[0]){
				szz[im]=tsz[0];
			}else if(xx>tsx[nxs-1]){
				szz[im]=tsz[nxs-1];
			}else{
                                l = binarySearch(xx,tsx,nxs);
                                oxs=tsx[l];
                                xs=xx-oxs;
                                szz[im]=coef[l*4+0]*xs*xs*xs+coef[l*4+1]*xs*xs+coef[l*4+2]*xs+coef[l*4+3];
                        }
                }
        }

	for(im=0;im<nsz[0];im++)
		sz[itf*nsz[0]+im]=szz[im];
	
	free(tsx);
	free(tsz);
	free(coef);
}

void smooth(float *data, int *n, int nrep, int *rect)
/*< Smooth velocity model using ENO interpolation >*/
{
	int dim = 2;
	int dim1 = 1;
	bool diff[2]={0,0};
	bool box[2]={0,0};
	int s[2];
	sf_triangle tr;
	int n1, n2;
	bool adj=false;
	int i, j, i0, irep;

	n1 = n2 = 1;
	for (i=0; i < dim; i++) {
		if (i <= dim1) {
			s[i] = n1; 
			n1 *= n[i];
		} else {
			n2 *= n[i];
		}
	}

	for (i=0; i <= dim1; i++) {
		if (rect[i] <= 1) continue;
			tr = sf_triangle_init (rect[i],n[i],box[i]);
		for (j=0; j < n1/n[i]; j++) {
			i0 = sf_first_index (i,j,dim1+1,n,s);
			for (irep=0; irep < nrep; irep++) {
				if (adj) {
					sf_smooth (tr,i0,s[i],diff[i],data);
				} else {
					sf_smooth2 (tr,i0,s[i],diff[i],data);
				}
			}
		}
		sf_triangle_close(tr);
	}
}

void updateVelocityModel(
			   float *vel,
			   int *n, /* Velocity model dimension n1=n[0] n2=n[1] */
			   float *o, /* Velocity model axis origin o1=o[0] o2=o[1] */
			   float *d, /* Velocity model sampling d1=d[0] d2=d[1] */
			   float *sv, /* Velocity model disturbance */
			   float *sz, /* Interfaces nodepoints */
			   int *nsz, /* Number of nodes, number of interfaces */
			   float *osz, /* Nodes origin, first interface index */
			   float *dsz, /* Nodes sampling, interface increment */
			   bool first, /* TODO remove this flag */
			   bool base, /* TODO Remove this flag */
			   int itf /* interface index */)
/*< Velocity model update
Note: This function uses a sv (layers velocity) vector and sz (depth interfaces
coordinates) vector to build the depth velocity model. There is nsv constant
velocity layers in the model and nsv-1 interfaces separating them.
These interfaces are described with nsz control points (nodes) in the sz vector and
they are interpolated using natural cubic spline interpolation.
 >*/
{
        int i, j, i1, i2, l=0;
	float *x=NULL, **coef=NULL;
	float zi[1];
	float xx, z;
	float *szz=NULL;
	#ifdef SMOOTH_VEL_MODEL
	int rect[2]={10,10}, nrep=1;
	#endif

	if(first){
		for(i2=0;i2<n[1];i2++){

			for(i1=0;i1<n[0];i1++){
					vel[i2*n[0]+i1]=sv[0];
			}
		}
	}else{

		x = sf_floatalloc(nsz[0]);
		szz = sf_floatalloc(nsz[0]);

		for(i=0;i<nsz[0];i++){
			x[i] = i*dsz[0]+osz[0];
			szz[i] = sz[(itf-1)*nsz[0]+i];
		}

	       /* Calculate coefficients matrix (interfaces interpolation) */
		coef = sf_floatalloc2(4*(nsz[0]-1),1);
		calculateSplineCoeficients(nsz[0],x,szz,coef[0]);
		// TODO remove this 'base' flag
		if(base){
			itf++; sf_warning("sv=%f",sv[itf]);
		}

		/* Calculate velocity function */
		for(j=0;j<n[1];j++){

			xx = d[1]*j+o[1];
			if(xx>x[l+1]) l++;
			/* Calculate interfaces z coordinates */
			calcInterfacesZcoord(zi,1,xx-x[l],l,coef);
			for(i=0;i<n[0];i++){
				z = i*d[0]+o[0];
				if(z>zi[0]){
					vel[n[0]*j+i] = sv[itf];
				}else{
					vel[n[0]*j+i] = sqrt(1./vel[n[0]*j+i]);
				}
			} /* Loop over depth */
		} /* Loop over distance */

		free(coef);
		free(x);
		free(szz);
	}

	#ifdef SMOOTH_VEL_MODEL
	smooth(vel,n,nrep,rect);
	#endif
}

void buildSlownessModelFromVelocityModel(
					 float *vel,
					 int *n, /* Velocity model dimension n1=n[0] n2=n[1] */
			 		 float *o, /* Velocity model axis origin o1=o[0] o2=o[1] */
					 float *d, /* Velocity model sampling d1=d[0] d2=d[1] */
					 float *sv, /* Velociy disturbance */
					 float *sz, /* Interfaces nodepoints */
					 int *nsz, /* Nodes per interface, number of interfaces */
					 float *osz, /* Nodes origin, interfaces first index */
					 float *dsz, /* Nodes sampling, interfaces increment */
					 bool first, /* TODO: remove this flag */
					 bool base, /* TODO: remove this flag */
					 int itf /* Interface index */)
/*< Slowness model build from velocity model
Note: This function is a function wrapper to updateVelocityModel function.
It calls that function to update the velocity model and build the slowness
model matrix using the slowness definition slow=(1.0/(v*v)). 
 >*/
{

	int i, nm; // Loop counters and indexes

	nm =n[0]*n[1];
	updateVelocityModel(vel,n,o,d,sv,sz,nsz,osz,dsz,first,base,itf);

	/* transform velocity to slowness */
	for(i=0;i<nm;i++){
			vel[i] = 1.0/(vel[i]*vel[i]);
	}
}

