ifeq ($(OS), Windows_NT)

all: bayesEpi.dll
	# comment
bayesEpi.dll: bayesEpi.cpp
	Rscript --vanilla -e "TMB::compile('bayesEpi.cpp','-O0 -g')"

else

all: bayesEpi.so
	# comment
bayesEpi.so: bayesEpi.cpp
	Rscript --vanilla -e "TMB::compile('bayesEpi.cpp','-O0 -g')"

endif


clean:
	rm -rf *o
