#{ R2RVSAO: SAO redshift package

set	r2rvsaobin	 = "r2rvsao$bin(arch)/"
package	r2rvsao, bin=r2rvsaobin$
 
task	bcvcorr,
	emsao,
	eqwidth,
	linespec,
	sumspec,
	velset,
	xcsao		= "r2rvsao$x_rvsao.e"

task	contpars	= "r2rvsao$contpars.par"
task	contsum		= "r2rvsao$contsum.par"
task	relearn		= "r2rvsao$relearn.cl"
task	qplot		= "r2rvsao$qplot.cl"
task	zvel		= "r2rvsao$zvel.cl"

# Write the welcome message
if (motd)
    type r2rvsao$r2rvsao.msg
;
clbye()
