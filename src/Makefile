ifeq ($(OS), Windows_NT)
all: cc.dll

else
all: cc.so

endif

cc.so: cc.cpp
	Rscript --vanilla -e "TMB::compile('cc.cpp','-O0 -g')"

clean:
	rm -rf *o
