# Makefile
# 
# Purpose: Run all tests with 'make'.
# 
# Site: https://www.geofisicando.com
# 
# Version 1.0
# 
# Programmer: Rodolfo A C Neves (Dirack) 05/10/2020
# 
# Email: rodolfo_profissional@hotmail.com
# 
# License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

DEP= ../raytrace.c ../grid2.c ../atela.c ../dynamic.c ../interface2d.c ../layer2d.c ../model2d.c
TESTS=test_rungekutta.x test_smoothlayer.x test_twopointrt.x test_semblance.x
EXP=test_model1 test_model3

all:	test exp

test:	$(TESTS)

exp:	$(EXP)

test_rungekutta.x:	test_rungekutta.c ../grid2.c ../atela.c ../interface2d.c Unity/unity.c
	gcc $^ -I$$RSFSRC/include -I$$RSFSRC/user/landa88 -L$$RSFSRC/lib -o $@ -lrsf -lm
	./$@ && rm $@

test_smoothlayer.x:	test_smoothlayer.c ../grid2.c ../atela.c ../interface2d.c Unity/unity.c
	gcc $^ -I$$RSFSRC/include -I$$RSFSRC/user/landa88 -L$$RSFSRC/lib -o $@ -lrsf -lm
	./$@ && rm $@

test_twopointrt.x:	test_twopointrt.c ../grid2.c ../atela.c ../interface2d.c Unity/unity.c
	gcc $^ -I$$RSFSRC/include -I$$RSFSRC/user/landa88 -L$$RSFSRC/lib -o $@ -lrsf -lm
	./$@ && rm $@

test_model1:    model1/SConstruct
	scons -C model1 nit=200 inverted-original.view

test_model3:    model3/SConstruct
	scons -C model3 nit=200 inverted-original.view

test_semblance.x:       test_semblance.c ../grid2.c ../atela.c Unity/unity.c
	gcc $^ -I$$RSFSRC/include -I$$RSFSRC/user/landa88 -L$$RSFSRC/lib -o $@ -lrsf -lm -lblas
	./$@ #> times.asc && echo "n1=101 d1=0.0125 o1=0 n2=1 d2=1 o2=0 data_format=ascii_float in=times.asc" >> times.asc && rm $@
	python3 plotsemb.py; python3 plotrnip.py

.PHONY:	$(TESTS)
