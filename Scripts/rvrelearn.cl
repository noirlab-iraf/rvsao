# File rvsao/Scripts/rvrelearn.cl
# By Jessica Mink, Center for Astrophysics
# May 25, 2012

# RVRELEARN Update parameter files for all of the tasks in the RVSAO package

procedure rvrelearn ()

bool verbose=yes {prompt="List tasks as parameters files are relearned (yes or no)"}

begin
    bool vb
    vb = verbose

# SPP Tasks
    relearn ("bcvcorr")
    if (vb)
	print ("RVRELEARN: BCVCORR relearned")
    relearn ("emsao")
    if (vb)
	print ("RVRELEARN: EMSAO relearned")
    relearn ("eqwidth")
    if (vb)
	print ("RVRELEARN: EQWIDTH relearned")
    relearn ("linespec")
    if (vb)
	print ("RVRELEARN: LINESPEC relearned")
    relearn ("pemsao")
    if (vb)
	print ("RVRELEARN: PEMSAO relearned")
    relearn ("pxcsao")
    if (vb)
	print ("RVRELEARN: PXCSAO relearned")
    relearn ("sumspec")
    if (vb)
	print ("RVRELEARN: SUMSPEC relearned")
    relearn ("wlrange")
    if (vb)
	print ("RVRELEARN: WLRANGE relearned")
    relearn ("listspec")
    if (vb)
	print ("RVRELEARN: LISTSPEC relearned")
    relearn ("pix2wl")
    if (vb)
	print ("RVRELEARN: PIX2WL relearned")
    relearn ("wl2pix")
    if (vb)
	print ("RVRELEARN: WL2PIX relearned")
    relearn ("velset")
    if (vb)
	print ("RVRELEARN: VELSET relearned")
    relearn ("xcsao")
    if (vb)
	print ("RVRELEARN: XCSAO relearned")

# CL PSETs
    relearn ("contpars")
    if (vb)
	print ("RVRELEARN: CONTPARS relearned")
    relearn ("contsum")
    if (vb)
	print ("RVRELEARN: CONTSUM relearned")

# CL Tasks
    relearn ("emplot")
    if (vb)
	print ("RVRELEARN: EMPLOT relearned")
    relearn ("qplot")
    if (vb)
	print ("RVRELEARN: QPLOT relearned")
    relearn ("setvel")
    if (vb)
	print ("RVRELEARN: SETVEL relearned")
    relearn ("qplotc")
    if (vb)
	print ("RVRELEARN: QPLOTC relearned")
    relearn ("skyplot")
    if (vb)
	print ("RVRELEARN: SKYPLOT relearned")
    relearn ("xcplot")
    if (vb)
	print ("RVRELEARN: XCPLOT relearned")
    relearn ("zvel")
    if (vb)
	print ("RVRELEARN: ZVEL relearned")
    relearn ("pvel")
    if (vb)
	print ("RVRELEARN: PVEL relearned")
    relearn ("qplotx")
    if (vb)
	print ("RVRELEARN: QPLOTX relearned")

    rvsao.version = rvsao.newversion
    if (vb)
	print ("RVRELEARN: RVSAO parameters updated to ",rvsao.newversion)
end

# Sep 15 2008	New task

# Jun 15 2009	Add WLRANGE

# May 25 2012	Add QPLOTX
