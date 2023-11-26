# File rvsao/Scripts/qplot.cl
# August 19, 1999
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# QPLOT -- Use EMSAO to plot emission and absorption lines already found
#            in a spectrum file or list of spectra

procedure qplot (spectra)

string	spectra=""	{prompt="RFN or file"}
string	qtask="xcsao"	{prompt="program to run (xcsao or emsao)"}
string	velplot="combination"	{prompt="Velocity to plot",
				 enum="combination|emission|correlation"}

begin

if (qtask == "emsao") {
    emsao (spectra,vel_init="combination",linefit=no,curmode=yes,save_vel=yes,vel_plot=velplot)
    }
else {
    xcsao (spectra,correlate=no,curmode=yes,dispmode=2,save_vel=yes,vel_plot=velplot)
    }
 
end

# Sep 28 1995	New script

# Jul 16 1999	Fix so headers are updated properly
# Aug 19 1999	Add velplot parameter
