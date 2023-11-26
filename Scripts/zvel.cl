# File rvsao/zvel.cl
# February 6, 2002
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# ZVEL -- Compute correlation and emission line velocities for
#         a file or list of files using RVSAO package

procedure zvel (imagelist)

string imagelist=""		{prompt="file or files"}
string imagext=".imh"		{prompt="Image header file name extension"}
bool emis_vel=yes		{prompt="Compute emission line velocities (y or n)"}
bool corr_vel=yes		{prompt="Compute cross correlation velocities (y or n)"}
bool plot=yes			{prompt="Plot results on display (y or n)"}
bool hard_copy=no		{prompt="Make printer hard copies (y or n)"}
bool curmode=no			{prompt="Wait after plot with cursor (y or n)"}
bool verbose=no			{prompt="Print what's happening (y or n)"}
bool debug=no			{prompt="Print everything that happens (y or n)"}

begin

bool vb, corvel,emvel, plres, hcres, curwait, dellist
int ifile, nfiles, icom, lname
file listfile,iraffile,irafhead
string iml, tchar,imx, imhl, iml1

	vb = verbose
	corvel = corr_vel
	emvel = emis_vel
	plres = plot
	hcres = hard_copy
	curwait = curmode
	iml = imagelist
	imx = imagext

#  Make list of files from rfn range or file of rfn's

#  If first character of imagelist is @, get number of images in list
        tchar = substr (iml,1,1)
	icom = stridx (',',iml)
        if (tchar == "@") {
            listfile = substr (iml,2,100)
	    dellist = no
            }

#  If comma-delimited list, list to file, adding image header extension
	else if (icom > 0) {
	    imhl = ""
	    while (icom > 0) {
		imhl = imhl // substr (iml,1,icom-1)
		imhl = imhl // imx // ","
		iml = substr (iml,icom+1,100)
		icom = stridx (',',iml)
		print (icom," ",iml)
		}
	    imhl = imhl // iml // imx
	    listfile = mktemp ("zvel")
	    files (imhl,>listfile)
	    dellist = yes
	    }

#  Else list to file, adding image header extension
	else {
	    listfile = mktemp ("zvel")
	    iml = iml // imx
	    files (iml,>listfile)
	    dellist = yes
	    } 
        if (vb)
            print ("ZVEL: Reading filenames from ",listfile)

#  Set up parameters in archive writing, graphic display, and hard copies

	if (corvel) {
	    if (emvel)
		xcsao.save_vel=yes
	    xcsao.displot = plres
	    xcsao.hardcopy = hcres
	    xcsao.curmode = curwait
	    xcsao.debug = debug
	    }

	if (emvel) {
	    emsao.displot = plres
	    emsao.hardcopy = hcres
	    emsao.curmode = curwait
	    emsao.debug = debug
	    }

#  Loop through reduced files in list

	list = listfile
	while (fscan (list,irafhead) != EOF) {
	    if (!access (irafhead)) {
		irafhead = irafhead // imx
		if (!access (irafhead)) {
		    print ("ZVEL: Cannot read file ",irafhead)
		    next
		    }
		lname = strlen (irafhead) - 4
		iraffile = substr (irafhead,1,lname)
		}
	    else
		iraffile = irafhead
	    if (vb)
		print ("ZVEL: Processing ",iraffile)

#	Find cross-correlation velocity
	    if (corvel)
		xcsao (iraffile)

#	Find emission line velocity
	    if (emvel)
		emsao (iraffile)
	    }

	if (dellist)
	    delete (listfile)
end
# May 23 1994	Read using SCAN

# Feb  6 1995	Add extension only if it is not already present

# Feb  6 2002	Fix bugs found testing on 2dF spectra
