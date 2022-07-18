# Makefile
# 
# Purpose: Run all tests with 'make'.
# 
# Site: https://www.geofisicando.com
# 
# Version 1.0
# 
# Programmer: Rodolfo A C Neves (Dirack) 07/07/2022
# 
# Email: rodolfo_profissional@hotmail.com
# 
# License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

DEP= ../raytrace.c ../grid2.c ../atela.c ../dynamic.c
SRCLANDA= ../Mlanda88_smooth.c ../velocity_lib.c ../tomography.c ../vfsacrsnh_lib.c

all:

landa88_smooth:	landa88_smooth.x GDB1

landa88_smooth.x:
	@echo "Program landa88tomo.x for GDB built"
	gcc -DGDB_DEBUG -g $(SRCLANDA) $(DEP) -I$$RSFSRC/include -I.. -L$$RSFSRC/lib -o $@ -lrsf -lm -lblas

GDB1:	landa88_smooth.x
	echo "Running GDB script for landa88tomo.x"
	gdb -q --command=gdb_scripts/gdb_script1 -ex "r" --args $<

clean:
	rm *.x

.PHONY: clean landa88_smooth.x

help:
	@echo "Use 'make' to run unit tests and intergration test"
	@echo "Use 'make progname' to generate a GDB version of progname for debug"
	@echo "Programs available:"
	@echo "\tlanda88smooth"

