# File rvsao/Scripts/qplotx.cl
# May 25, 2012
# By Jessica Mink, Harvard-Smithsonian Center for Astrophysics

# QPLOTX -- Use XCSAO or EMSAO to plot emission and absorption lines already found
#            in a spectrum file or list of spectra and set wavelength limits

procedure qplotx (spectra)

string	spectra=""	{prompt="RFN or file"}
string	qtask="emsao"	{prompt="program to run (xcsao or emsao)"}
string	velplot="combination"	{prompt="Velocity to plot",
				 enum="combination|emission|correlation"}
real	w1=INDEF	{prompt="Starting wavelength in angstroms INDEF=file"}
real	w2=INDEF	{prompt="Ending wavelength in angstroms INDEF=file"}
int	dispmode=4	{prompt="Display mode (2 or 4-continuum, bad lines)"}

begin

if (qtask == "emsao") {
    emsao (spectra,vel_init="combination",linefit=no,st_lambda=w1,end_lambda=w2,dispmode=dispmode,curmode=yes,save_vel=yes,vel_plot=velplot)
    }
else {
    xcsao (spectra,correlate=no,st_lambda=w1,end_lambda=w2,curmode=yes,dispmode=dispmode,save_vel=yes,vel_plot=velplot)
    }
 
end

# May 25 2012	New script based on QPLOTC
