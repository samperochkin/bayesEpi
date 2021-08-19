all: cc.so
	# comment

cc.so: cc.cpp
	Rscript --vanilla -e "TMB::compile('cc.cpp','-O0 -g')"
# anothertemplate.so: anothertemplate.cpp
	# Rscript --vanilla -e "TMB::compile('anothertemplate.cpp','-O0 -g')"

clean:
	rm -rf *o
