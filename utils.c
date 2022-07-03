/*
	 utils.c (c)
	 
	 Purpose: TODO
	 	 
	 Version 1.0
	 
	 Site: https://dirack.github.io
	 
	 Programmer: Rodolfo A. C. Neves (Dirack) 14/12/2021

	 Email:  rodolfo_profissional@hotmail.com

	 License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

*/

#include <rsf.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"

void dumpfloat1(char* l, float* v, int n)
/*< dump vector >*/
{
	int i;

	for(i=0;i<n;i++)
		sf_warning("%s[%d]=%f",l,i,v[i]);
}

