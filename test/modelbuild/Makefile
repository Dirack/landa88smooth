# Makefile
# 
# Purpose: Compile the programs to build smooth velocity models.
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

all:	model2d_smooth.x

#TODO
#model1d.x:	model1d.c velocity.h
#	gcc $< -I$$RSFSRC/include -L$$RSFSRC/lib -o $@ -lrsf -lm

#model2d.x:	model2d.c velocity.h
#	gcc $< -I$$RSFSRC/include -L$$RSFSRC/lib -o $@ -lrsf -lm

model2d_smooth.x:	model2d_smooth.c velocity.h
	gcc $< -I$$RSFSRC/include -L$$RSFSRC/lib -o $@ -lrsf -lm -lblas

clean:
	rm *.x
