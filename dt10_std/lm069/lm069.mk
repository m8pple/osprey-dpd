
%/lm069.0.state : %/dmpci.lm069
	(cd $* && ../../../build/dpd lm069)

clean_intermediates :
	rm */*.vtk */dmpchs.* */dmpcls.* */dmpcis.* */dmpcas.* 