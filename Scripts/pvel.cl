# File rvsao/pvel.cl
# November 2, 2012
# By Jessica Mink, Harvard-Smithsonian Center for Astrophysics

# PVEL -- Compute correlation velocities for
#         a file or list of files using rvsao.pxcsao

procedure pvel (spectra)

string spectra=""	{prompt="Spectrum file(s) or list file"}
int nskip=0		{prompt="Number of files in list to skip"}
int ap1=1		{prompt="First spectrum aperture"}
int ap2=1		{prompt="Last spectrum aperture"}
int specband=0		{prompt="Spectrum band if multispec file"}
int specext=0		{prompt="FITS extension for spectra (0 if none)"}
string specdir=""	{prompt="Spectrum directory"}
string templates=""	{prompt="List of template spectra"}
string tempdir=""	{prompt="Template directory"}
bool emis_vel=no	{prompt="Compute emission line velocities (y or n)"}
bool corr_vel=yes	{prompt="Compute cross correlation velocities (y or n)"}
bool plot=no		{prompt="Plot results on display (y or n)"}
bool hard_copy=no	{prompt="Make printer hard copies (y or n)"}
bool curmode=no		{prompt="Wait after plot with cursor (y or n)"}
bool verbose=no		{prompt="Print what's happening (y or n)"}
bool debug=no		{prompt="Print everything that happens (y or n)"}

begin

bool vb, corvel,emvel, plres, hcres, curwait, dellist, first
int ifile, nfiles, nfsk
file listfile,specfile, specpath
string iml, tchar

	vb = verbose
	corvel = corr_vel
	emvel = emis_vel
	plres = plot
	hcres = hard_copy
	curwait = curmode
	iml = spectra
	nfsk = nskip
	first = yes

#  Make list of files from rfn range or file of rfn's

#  If first character of spectra is @, get number of images in list
        tchar = substr (iml,1,1)
        if (tchar == "@") {
            listfile = substr (iml,2,100)
	    dellist = no
            }

#  Else list to file, adding image header extension
	else {
	    listfile = mktemp ("zvel")
	    files (iml,>listfile)
	    dellist = yes
	    } 
        if (vb)
            print ("PVEL: Reading filenames from ",listfile)

#  Set up parameters in archive writing, graphic display, and hard copies

	if (corvel) {
	    if (emvel)
		pxcsao.save_vel=yes
	    pxcsao.specnum = ap1 //"-" // ap2
	    pxcsao.specext = specext
	    pxcsao.specband = specband
	    pxcsao.specdir = specdir
	    pxcsao.templates = templates
	    pxcsao.tempdir = tempdir
	    pxcsao.displot = plres
	    pxcsao.hardcopy = hcres
	    pxcsao.curmode = curwait
	    pxcsao.debug = debug
	    }

	if (emvel) {
	    pemsao.displot = plres
	    pemsao.hardcopy = hcres
	    pemsao.curmode = curwait
	    pemsao.debug = debug
	    }

#  Loop through reduced files in list

	ifile = 0
	list = listfile
	while (fscan (list,specfile) != EOF) {
	    if (ifile < nfsk) {
		ifile = ifile + 1
		next
		}
	    specpath = specdir // specfile
	    if (!access (specpath)) {
		print ("PVEL: Cannot read file ",specpath)
		ifile = ifile + 1
		next
		}
	    if (vb) {
		if (ap2 > ap1)
		    print ("PVEL: Processing ",specfile,"[",ap1,"-",ap2,"]")
		else
		    print ("PVEL: Processing ",specfile)
		}
	    if (ifile < 10) {
		pxcsao.logfiles = "STDOUT,pxcsao00" // ifile // ".log"
		}
	    else if (ifile < 100) {
		pxcsao.logfiles = "STDOUT,pxcsao0" // ifile // ".log"
		}
	    else {
		pxcsao.logfiles = "STDOUT,pxcsao" // ifile // ".log"
		}

#	Find cross-correlation velocity
	    if (corvel)
		pxcsao (specfile)

#	After processing first spectrum file, reset aperture range
	    pxcsao.specnum = "1-" // ap2

#	Find emission line velocity
	    if (emvel)
		pemsao (specfile)

	    ifile = ifile + 1
	    }

	if (dellist)
	    delete (listfile)
end
# Jul 22 2010	New task after ZVEL
# Jul 23 2010	Replace specnum with ap1 and ap2; add nskip for restarting
# Aug 11 2010	Zero-pad logfile number to 3 digits

# Nov  2 2012	Fix procedure name and setting of save_vel
