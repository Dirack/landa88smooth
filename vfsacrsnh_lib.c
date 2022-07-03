/*
	 vfsacrsnh_lib.c (c)
	 
	 Purpose: 'Mvfsacrsnh.c' library.
	 	 
	 Version 1.0
	 
	 Site: https://dirack.github.io
	 
	 Programmer: Rodolfo A. C. Neves (Dirack) 19/09/2019

	 Email:  rodolfo_profissional@hotmail.com

	 License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

*/

/*
TODO: Modify macro definition in search window for each interface.
Large windows can make the result oscilate a lot and do not converge
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <rsf.h>
/*^*/

#define signal(s) ((s<0)?(-1.):(1.))
/*< Signal function >*/
/*^*/

float getRandomNumberBetween0and1(){
/*< Function to get a random number between 0 and 1 >*/

	return (float)(rand()%1000)/1000;
}

float getVfsaIterationTemperature(int iteration,float dampingFactor,float inicialTemperature){
/*< Temperature function for VFSA algorithm >*/

	return inicialTemperature*expf(-dampingFactor*pow(iteration,0.25));

}

/* TODO: Modify this function for multiple interfaces */
void disturbParameters( float temperature, /* Temperature of this interation in VFSA */
			float *disturbedVel, /* Parameters disturbed vector */
			float *vel,
			float minvel,
			float maxvel,
			float scale, /* Scale to multiply by disturbance */
			int itf)
/*< Disturb parameters from the previous iteration of VFSA
 Note: It receives a parameter vector and distubs it accordingly to 
VFSA disturb parameters step.
 >*/
{

	float u;
	float disturbance;
	int i,j;

	u=getRandomNumberBetween0and1();
				
	disturbance = signal(u - 0.5) * temperature * (pow( (1+temperature),fabs(2*u-1) )-1);
	disturbance *= scale;

	disturbedVel[itf] = vel[itf] + (disturbance);

	if (disturbedVel[itf] >= maxvel)
		disturbedVel[itf] = maxvel - (maxvel-minvel) * getRandomNumberBetween0and1();

	if (disturbedVel[itf] <= minvel)
		disturbedVel[itf] = (maxvel-minvel) * getRandomNumberBetween0and1() + minvel;
}

