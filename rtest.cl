#{ RVSAO: SAO redshift package

set	rvsaobin	 = "rvsao$bin(arch)/"
set	rvscripts	 = "rvsao$Scripts/"
package	rvsao, bin=rvsaobin$
string	test
 
task	bcvcorr,
	emsao,
	eqwidth,
	linespec,
	pemsao,
	pxcsao,
	sumspec,
	wlrange,
	listspec,
	pix2wl,
	wl2pix,
	velset,
	xcsao		= "rvsao$x_rvsao.e"

task	contpars	= "rvsao$contpars.par"
task	contsum		= "rvsao$contsum.par"
task	emplot		= "rvscripts$emplot.cl"
task	qplot		= "rvscripts$qplot.cl"
task	relearn		= "rvscripts$relearn.cl"
task	rvrelearn	= "rvscripts$rvrelearn.cl"
task	setvel		= "rvscripts$setvel.cl"
task	qplotc		= "rvscripts$qplotc.cl"
task	qplotx		= "rvscripts$qplotx.cl"
task	skyplot		= "rvscripts$skyplot.cl"
task	xcplot		= "rvscripts$xcplot.cl"
task	zvel		= "rvscripts$zvel.cl"
task	pvel		= "rvscripts$pvel.cl"

rvsao.newversion = "RVSAO 2.8.2, May 2, 2014"

# Write the welcome message
if (motd) {
    test = rvsao.newversion
    if (rvsao.newversion == rvsao.version)
        type rvsao$rvsao.msg
    else
        type rvsao$rvupdate.msg
    }
;
clbye()
