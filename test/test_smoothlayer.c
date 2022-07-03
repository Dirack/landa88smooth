/*
* test_smoothlayer.c (C)
* 
* Purpose: Get RNIP using dynamic ray tracing system for a smooth velocity model
* with curved interfaces.
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
It uses a smooth velocity model as input. The purpose is to obtain RNIP using
dynamic ray tracing system.
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

void init()
/*< Load velocity model >*/
{
	int nm=n[0]*n[1];
	int im;
	for(im=0;im<nm;im++)
		slow[im]=1./(slow[im]*slow[im]);
}

void setUp(){}

void tearDown(){}

void test_getRnipIniCurvedInterfaceLayerForAFanOfRays()
/*< Get RNIP using dynamic ray tracing system >*/
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

        x[0] = 1.85;
        x[1] = 5.;

	TEST_IGNORE_MESSAGE("TODO");

        // Test for 30 degrees ray
        /* initialize ray tracing object */
        rt = raytrace_init(2,true,nt,dt,n,o,d,slow,ORDER);
        traj = sf_floatalloc2(2,nt+1);

        p[0]=-cosf(-SF_PI/10000); p[1]=sinf(-SF_PI/10000);
        ei=-SF_PI/10000;
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
        freopen("modelbuild/vel3.rsf","r",stdin);

        sf_init(argc,argv);
        in = sf_input("in");

        /* Read velocity model */
        slow=sf_floatalloc(n[0]*n[1]);
        sf_floatread(slow,n[0]*n[1],in);

        init();
        UNITY_BEGIN();
        RUN_TEST(test_getRnipIniCurvedInterfaceLayerForAFanOfRays);
        return UNITY_END();
}

