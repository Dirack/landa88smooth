/*
* test_rungekutta.c (C)
* 
* Purpose: Calculate RNIP using dynamic ray tracing system.
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
It uses a smooth velocity model as input. The purpose is to calculate RNIP parameter
using dynamic ray tracing system. This velocity model has constant velocity layers
and plane interfaces, v={1.508,1.69,2.}
*/

#include "Unity/unity.h"
#include "raytrace.h"
#include "forward.h"
#include "rungekutta.h"
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

void test_getRNIPUsingDynamicRayTracingInConstantVelocityModel()
/*< Get RNIP for the first layer using dynamic ray tracing system
For the first layer, RNIP is the normal ray length
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

        s = sf_floatalloc2(2,1);

        /* initialize ray tracing object */
        rt = raytrace_init(2,true,nt,dt,n,o,d,slow,ORDER);
        traj = sf_floatalloc2(2,nt+1);
        normalRayAngleRad = a[0]*DEG2RAD;

        for(i=0;i<5;i++){
                s[0][0] = 1.;
                s[0][1] = i*0.5+2.5;
                /* Set initial ray point and ray vector */
                setInitialRayPointAndRayVector(s,x,p,0,normalRayAngleRad);
               /* Ray tracing */
                it = trace_ray (rt, x, p, traj);

                rnip = calculateRNIPWithDynamicRayTracing(rt,dt,it,traj,v0);

                TEST_ASSERT_FLOAT_WITHIN(0.1,1.0,rnip);
        }

        raytrace_close(rt);
        free(traj);
	free(s);
}

void test_getRNIPUsingDynamicRayTracingInTwoLayersVelocityModel()
/*< REPEAT THE FIRST TEST FOR THE SECOND INTERFACE
Now, RNIP is transmited through interfaces by dynamic ray system and
should be closer to the RNIP value calculated using Hubral's propagation
and transmission laws
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

        s = sf_floatalloc2(2,1);

        /* initialize ray tracing object */
        rt = raytrace_init(2,true,nt,dt,n,o,d,slow,ORDER);
        traj = sf_floatalloc2(2,nt+1);
        normalRayAngleRad = a[0]*DEG2RAD;
        for(i=0;i<5;i++){
                s[0][0] = 1.85;
                s[0][1] = i*0.5+2.5;
                /* Set initial ray point and ray vector */
                setInitialRayPointAndRayVector(s,x,p,0,normalRayAngleRad);

                /* Ray tracing */
                it = trace_ray (rt, x, p, traj);

                rnip = calculateRNIPWithDynamicRayTracing(rt,dt,it,traj,v0);

                TEST_ASSERT_FLOAT_WITHIN(0.1,1.95,rnip);
        }
        raytrace_close(rt);
        free(traj);
	free(s);
}

void test_getRnipInInterfaceLayerForAFanOfRays()
/*< REPEAT THE SECOND TEST FOR A BENDING RAY
The normal ray is tilted 30 degrees and launched. RNIP from dynamic ray tracing system
should agree with the RNIP value calculated using Hubral laws. This time, we should
consider incident and transmited angles in Hubrals transmission law.
>*/
{
        float x[2];
        float p[2];
        raytrace rt;
        float **traj;
        int nt=10000;
        float dt=0.001;
        float normalRayAngleRad;
        int it;
        float rnip;
        float vv0=1.508;
        float ei, et, rtt, ri;

        x[0] = 2.5;
        x[1] = 5.;

        // Test for 30 degrees ray
        /* initialize ray tracing object */
        rt = raytrace_init(2,true,nt,dt,n,o,d,slow,ORDER);
        traj = sf_floatalloc2(2,nt+1);

        p[0]=-cosf(SF_PI/5.); p[1]=sinf(SF_PI/5.);
        ei=SF_PI/5.;
        et=asin((1.508/1.69)*sin(ei));
        ri=1.508/cos(ei);
        rtt=(1.69/1.508)*ri*((cos(et)*cos(et))/(cos(ei)*cos(ei)));

	/* Ray tracing */
        it = trace_ray (rt, x, p, traj);

        // Test rt after transmission
        rnip = calculateRNIPWithDynamicRayTracing(rt,dt,it,traj,1.508);
        TEST_ASSERT_FLOAT_WITHIN(0.1,rtt+1./cos(et),rnip);

        raytrace_close(rt);
        free(traj);
}

int main(int argc, char* argv[]){

        /* Redirect the stdin to velocity model */
        freopen("modelbuild/vel.rsf","r",stdin);

        sf_init(argc,argv);
        in = sf_input("in");

        /* Read seismic velocity model */
        slow=sf_floatalloc(n[0]*n[1]);
        sf_floatread(slow,n[0]*n[1],in);

        init();
        UNITY_BEGIN();
        RUN_TEST(test_getRNIPUsingDynamicRayTracingInConstantVelocityModel);
	RUN_TEST(test_getRNIPUsingDynamicRayTracingInTwoLayersVelocityModel);
	RUN_TEST(test_getRnipInInterfaceLayerForAFanOfRays);
        return UNITY_END();
}

