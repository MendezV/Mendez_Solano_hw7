#Variable: Se definen acá y se usan con $(NOMBRE)
enrg=0.01
all:  estado_$(enrg).pdf

# $@ significa el archivo target
# $^ significa todos los prerrequistos 
sod_test.x:sod_test.c
	echo "Creando archivo .x "
	cc $^ -lm -o $@
estado_$(enrg).dat:sod_test.x
	echo "Creando archivo .dat"
	./$^ $(enrg) 
estado_$(enrg).pdf:estado_$(enrg).dat plot.py
	echo "Creando archivo .pdf  "
	python plot.py $(enrg)
