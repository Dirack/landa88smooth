/*
* test_twopointrt.c (C)
* 
* Purpose: Two point ray tracing using dynamic ray tracing system.
* 
* Site: https://www.geofisicando.com
* 
* Versão 1.0
* 
* Programmer: Rodolfo A C Neves (Dirack) 02/17/2022
* 
* Email: rodolfo_profissional@hotmail.com
* 
* License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.
*/

/*
It uses a smooth velocity model as input. The purpose is to make the normal ray
converge for the receiver coordinate in the acquisition surface using quantities
calculated through dynamic ray tracing. It will be 'q' and 'BETA' parameters.
*/

#include "Unity/unity.h"
#include "raytrace.h"
#include "forward.h"
#include "rungekutta2.h"
#include <stdio.h>
#include <rsf.h>

sf_file in;
int n[2]={301,1001};
float d[2]={0.01,0.01};
float o[2]={0.,-2.};
float v0=1.508;
float *slow;
float *slow2;

void init()
/*< Load smooth velocity model >*/
{
	int nm=n[0]*n[1];
	int im;
	for(im=0;im<nm;im++)
		slow[im]=1./(slow[im]*slow[im]);
}

void setUp(){}

void tearDown(){}

void test_twoPointRayTracingInConstantVelocityModel()
/*< Test two point ray tracing convergence
Use Two point ray tracing formula for angle perturbation to get the ray
that arrives at m0=2km in acquisition surface.
The ray is a normal ray that arrives at m0=2.5Km at acquisition surface. Use
the dynamic ray tracing to get q and calculate the angle increment dphi to trace
a new ray that should be closer to the receiver location R=2.Km. Repeat this procedure
until m0 converge to receiver position.
>*/
{
        float x[2];
        float p[2];
        raytrace rt;
        float **traj;
        float normalRayAngleRad;
        float a[1]={0.,};
        int nt=10000;
        float dt=0.001;
        int it;
        float t;
        int is;
        float **s;
        float rnip=0.;
        int i;
	float q;
	float y[2];
	float beta;
	float dphi;
	float m0;

        s = sf_floatalloc2(2,1);

        /* initialize ray tracing object */
        rt = raytrace_init(2,true,nt,dt,n,o,d,slow,ORDER);
        traj = sf_floatalloc2(2,nt+1);
        normalRayAngleRad = a[0]*DEG2RAD;

        for(i=0;i<7;i++){
                s[0][0] = 1.;
                s[0][1] = 2.5;
                /* Set initial ray point and ray vector */
                setInitialRayPointAndRayVector(s,x,p,0,normalRayAngleRad);
               /* Ray tracing */
                it = trace_ray (rt, x, p, traj);

                rnip = calculateRNIPWithDynamicRayTracing2(rt,dt,it,traj,v0,&q);
		y[0] = traj[it][0]; y[1] = traj[it][1];
		beta = calculateBetaWithRayTrajectory(y,traj,it);
		dphi = (-cosf(3.1415-beta))*(2.-traj[it][1])/(0.44*q);

		m0 = traj[it][1];
		printf("m0=%f R=2. diff=%f\n",m0,fabs(traj[it][1]-2.));
		printf("dphi=%f RNIP=%f BETA=%f\n",dphi,rnip,beta);

        	normalRayAngleRad += dphi;
        }

	TEST_ASSERT_FLOAT_WITHIN(0.01,2.,m0);

        raytrace_close(rt);
        free(traj);
	free(s);
}

void test_twoPointRayTracingInTwoLayersVelocityModel()
/*< *** REPEAT THE PERVIOUS TEST FOR A TWO LAYERS MODEL ***
Test two point ray tracing convergence
Use Two point ray tracing formula for angle perturbation to get the ray
that arrives at m0=2km in acquisition surface.
The ray is a normal ray that arrives at m0=2.5Km at acquisition surface. Use
the dynamic ray tracing to get q and calculate the angle increment dphi to trace
a new ray that should be closer to the receiver location R=2.Km. Repeat this procedure
until m0 converge to receiver position.
>*/
{
        float x[2];
        float p[2];
        raytrace rt;
        float **traj;
        float normalRayAngleRad;
        float a[1]={10,};
        int nt=10000;
        float dt=0.001;
        int it;
        float t;
        int is;
        float **s;
        float rnip=0.;
        int i;
	float q;
	float y[2];
	float beta;
	float dphi;
	float m0;

        s = sf_floatalloc2(2,1);

        /* initialize ray tracing object */
        rt = raytrace_init(2,true,nt,dt,n,o,d,slow,ORDER);
        traj = sf_floatalloc2(2,nt+1);
        normalRayAngleRad = a[0]*DEG2RAD;

        for(i=0;i<5;i++){
                s[0][0] = 1.85;
                s[0][1] = 2.5;
                /* Set initial ray point and ray vector */
                setInitialRayPointAndRayVector(s,x,p,0,normalRayAngleRad);
               /* Ray tracing */
                it = trace_ray (rt, x, p, traj);

                rnip = calculateRNIPWithDynamicRayTracing2(rt,dt,it,traj,v0,&q);
		y[0] = traj[it][0]; y[1] = traj[it][1];
		beta = calculateBetaWithRayTrajectory(y,traj,it);
		dphi = (-cosf(3.1415-beta))*(2.-traj[it][1])/(0.44*q);

		m0 = traj[it][1];
		printf("m0=%f R=2. diff=%f\n",m0,fabs(traj[it][1]-2.));
		printf("dphi=%f RNIP=%f BETA=%f\n",dphi,rnip,beta);

        	normalRayAngleRad += dphi;
        }

	TEST_ASSERT_FLOAT_WITHIN(0.01,2.,m0);

        raytrace_close(rt);
        free(traj);
	free(s);
}


int main(int argc, char* argv[]){

        /* Redirect the stdin to RSF velocity file */
        freopen("modelbuild/vel3.rsf","r",stdin);

        sf_init(argc,argv);
        in = sf_input("in");

        /* Read velocity file */
        slow=sf_floatalloc(n[0]*n[1]);
        sf_floatread(slow,n[0]*n[1],in);

        init();
        UNITY_BEGIN();
        RUN_TEST(test_twoPointRayTracingInConstantVelocityModel);
	RUN_TEST(test_twoPointRayTracingInTwoLayersVelocityModel);
        return UNITY_END();
}

