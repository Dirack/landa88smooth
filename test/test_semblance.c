#include "Unity/unity.h"
#include "raytrace.h"
#include "forward.h"
#include "datamis.h"
#include "rungekutta.h"
#include "velocity.h"
#include <time.h>
#include <stdio.h>
#include <rsf.h>

int n[2]={301,1001};
float d[2]={0.01,0.01};
float o[2]={0.,-2.};
float v0=1.508;
float *slow;
float ***data;
int data_n[3]={1001,161,482};
float data_o[3]={0.,0.,0.};
float data_d[3]={0.004,0.025,0.00625};
sf_file in; // datacube RSF file
FILE *SEMB;
FILE *RNIP;
FILE *VEL;
float sz[42]={1.254316,1.254316,1.254316,1.254316,1.254316,1.254316,1.254316,1.259718,1.263270,1.208849,1.117164,1.020568,0.947812,0.926439,0.926439,0.926439,0.926439,0.926439,0.926439,0.926439,0.926439,1.959288,1.959288,1.959288,1.959288,1.959288,1.959288,1.959288,1.959288,1.965899,1.969283,1.967570,1.968356,1.974897,1.975764,1.975764,1.975764,1.975764,1.975764,1.975764,1.975764,1.975764};
int nsz[2]={21,2};
float osz[2]={-2,0};
float dsz[2]={0.5,1.};

void init(){
       int i, j;

        for(i=0;i<n[0]*n[1];i++){
                slow[i]=1./(slow[i]*slow[i]);
        }

}

void setUp(){};

void tearDown(){};

void test_semblanceTwoLayersModel()
/*< Test model setup in a constant velocity model: NIP sources position and NIP angles >*/
{
	float **s;
        int ns=1;
        float m0[1]={0.};
        float t0[1]={2.183};
        float a[1]={0.};
	float med[9];
	float average=0.;;
        float BETA;
        float v;
	int nv=50;
	float ov=1.508;
	float dv=0.01;
	int i, j;
	float lambda=1.;
	float vm=1.9;
	float semb=0.;
	float sumAmplitudes, sumAmplitudes2;
	int numSamples=1;
        float x[2];
        float p[2];
        raytrace rt;
        float **traj;
        float normalRayAngleRad;
        int nt=10000;
        float dt=0.001;
        int it;
        float t;
        int is;
        float rnip=0.;
	float sv[3]={1.508,1.69,2.};

        s = sf_floatalloc2(1,ns);

	s[0][0]=1.85;
	s[0][1]=5.;
	m0[0]=s[0][1];
	t0[0]=2.33;
	BETA=0.;

	/* initialize ray tracing object */
	updateVelocityModel(slow, n, o, d, sv, sz, nsz, osz, dsz, false, false,1);
	init();
        //rt = raytrace_init(2,true,nt,dt,n,o,d,slow,ORDER);
        traj = sf_floatalloc2(2,nt+1);
        normalRayAngleRad = a[0]*DEG2RAD;

	//for(i=0;i<n[0];i++)
	//	printf("i=%d v=%f\n",i,slow[i]);

	for(i=0;i<nv;i++){
		v=dv*i+ov;
		sv[1]=v;	
		updateVelocityModel(slow, n, o, d, sv, sz, nsz, osz, dsz, false, false,1);
		init();
		rt = raytrace_init(2,true,nt,dt,n,o,d,slow,ORDER);
		//for(j=0;j<n[0];j++)
		//printf("j=%d v=%f\n",j,slow[j]);

		/* Set initial ray point and ray vector */
                setInitialRayPointAndRayVector(s,x,p,0,normalRayAngleRad);

                /* Ray tracing */
                it = trace_ray (rt, x, p, traj);

		rnip=calculateRNIPWithDynamicRayTracing(rt,dt,it,traj,v0);
		sumAmplitudes=0., sumAmplitudes2=0.;
                numSamples = stackOverCRETimeCurve(rnip,BETA,traj[it][1],2*it*dt,v0,&sumAmplitudes,&sumAmplitudes2,data,data_n,data_o,data_d);
		sf_warning("rnip=%f t=%f\n",rnip,it*dt*2);
		semb = (sumAmplitudes*sumAmplitudes)/(0.004*numSamples*sumAmplitudes2);
		fprintf(VEL,"%f\n",v);
		fprintf(SEMB,"%f\n",semb);
		fprintf(RNIP,"%f\n",rnip);

		raytrace_close(rt);
	}

}


int main(int argc, char* argv[]){

	int i;

        /* Redirect the stdin to datacube file */
        //freopen("data/interpolatedDataCube2.rsf","r",stdin);
        freopen("modelbuild/interpolatedDataCube.rsf","r",stdin);
	VEL = fopen("vel.txt","w");
	SEMB = fopen("semb.txt","w");
	RNIP = fopen("rnip.txt","w");

        sf_init(argc,argv);
        in = sf_input("in");

        /* Read seismic data cube */
        data=sf_floatalloc3(data_n[0],data_n[1],data_n[2]);
        sf_floatread(data[0][0],data_n[0]*data_n[1]*data_n[2],in);

        slow = sf_floatalloc(n[0]*n[1]);
	for(i=0;i<n[0]*n[1];i++)
			slow[i]=1.508;
	init();
	test_semblanceTwoLayersModel();
	fclose(VEL);
	fclose(SEMB);
	fclose(RNIP);
        //return UNITY_END();
}

