# File setvel.cl
# February 8, 2011
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# SETVEL -- Shift a single spectrum to a variety of different radial velocities

procedure setvel (baseimage)

file	baseimage="image0.fits"	{prompt="Original zero-velocity image"}

string	rootout="im"		{prompt="Root for output image name"}
file	velfile="F2cq_em5.velocities"	{prompt="File of velocities to shift"}

begin

real	velocity
int	ifile
file	outfile

	sumspec.spectra = baseimage
	sumspec.specnum = "1"
	sumspec.specband = 1
	sumspec.specdir = ""
	sumspec.compname = "Sample Hectospec spectrum"
	sumspec.compdir = "Specimage"
	sumspec.nspec = 1
	sumspec.save_names = yes
	sumspec.copy_header = yes
	sumspec.normin = 0.
	sumspec.fixbad = "no"
	sumspec.badlines = "badlines.dat"
	sumspec.linedir = "hectospec$lib/"
	sumspec.cont_remove = "no"
	sumspec.cont_split = 1
	sumspec.reject = no
	sumspec.abs_reject = 2.
	sumspec.em_reject = 2.
	sumspec.contout = "no"
	sumspec.cont_plot = no
	sumspec.cont_add = 0.
	sumspec.spec_smooth = 0
	sumspec.st_lambda = 3700.
	sumspec.end_lambda = 9150.
	sumspec.pix_lambda = 1.
	sumspec.npts = 5451
	sumspec.complog = no
	sumspec.interp_mode = "spline3"
	sumspec.normout = 1.
	sumspec.spec_plot = no
	sumspec.spec_int = no
	sumspec.comp_plot = no
	sumspec.comp_int = no
	sumspec.ymin = 0.
	sumspec.ymax = 2.
	sumspec.velcomp = 0.
	sumspec.zcomp = INDEF
	sumspec.svel_corr = "none"
	sumspec.nsmooth = 0
	sumspec.device = "stdgraph"
	sumspec.plotter = "stdplot"
	sumspec.logfiles = "STDOUT"
	sumspec.nsum = 1
	sumspec.debug = no
	sumspec.cursor = ""

	ifile = 0
	list = velfile
	while (fscan (list,velocity)) {
	    ifile = ifile + 1
	    outfile = rootout // "_"
	    if (ifile < 10)
		outfile = outfile // "0"
	    if (ifile < 100)
		outfile = outfile // "0"
	    if (ifile < 1000)
		outfile = outfile // "0"
	    outfile = outfile // i // ".fits"
	    
	    sumspec (compfile=outfile,velcomp=velocity)
	    }

	print (ifile, " files created")
	return
end

# Feb 27 2006	New script

# Feb  8 2011	Change declaration of velocity to real
