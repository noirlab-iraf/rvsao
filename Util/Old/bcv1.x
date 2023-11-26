#*** File rvsao/util/bcv.x
#*** January 3, 1997
#*** By Doug Mink
#*** From G. Torres Fortran code (1989)

#  bcv calculates the correction required to reduce observed (topocentric)
#  radial velocities of a given star to the barycenter of the solar system.
#  it includes correction for the effect of the earth's rotation.
#  the maximum error of this routine is not expected to be larger than 0.6 m/sec

procedure bcv (djd,dlong,dlat,dalt,dra,dec,deq,bcvel,hcvel,gcvel)

double	djd	# Julian date (days)
double	dlong	# geodetic longitude (degrees, west positive)
double	dlat	# geodetic latitude (degrees)
double	dalt	# altitude above sea level (meters)
double	dra	# right ascension of star (hours)
double	dec	# declination of star (degrees)
double	deq	# mean equator and equinox for coordinates
			# e.g., 1950.0
double	bcvel	# barycentric correction (km/s) (returned)
double	hcvel	# heliocentric correction (km/s) (returned)
double	gcvel	# geocentric correction (km/s) (returned)

double	dlongs	# geodetic longitude (radians, west+)
double	dlats	# geodetic latitude (radians)
double	dras	# right ascension of star (radians)
double	decs	# declination of star (radians)

double	dc[3],dcc[3],dprema[3,3],dvelh[3],dvelb[3]
double	dpi,daukm,dctrop,dcbes,dc1900,dct0,dtr,dst
double	deqt,dra2,dec2,dha,darg,tpi
int	k

begin
	dpi = 3.1415926535897932d0
	tpi = 2.d0 * dpi
	daukm = 1.4959787d08
	dctrop = 365.24219572d0
	dcbes = 0.313d0
	dc1900 = 1900.0d0
	dct0 = 2415020.0d0
	dtr = dpi / 180.0d0

#	open (15,'testbcv')
#	write(15,*) 'long:',dlong,' lat:',dlat
#	write(15,*) 'ra:',dra,' dec:',dec,deq

#  calculate local sidereal time

	dlongs = dlong*dtr
	call sidtim (djd,dlongs,dst)

#  precess r.a. and dec. to mean equator and equinox of date (deqt)

	deqt = (djd - dct0 - dcbes)/dctrop + dc1900
	dras = dra * 15.0d0 * dtr
	decs = dec * dtr
	dc[1] = cos(dras) * cos(decs)
	dc(2) = sin(dras) * cos(decs)
	dc(3) =	     sin(decs)

	call pre (deq,deqt,dprema)
	do k = 1, 3 {
	    dcc[k]=dc[1]*dprema[k,1]+dc[2]*dprema[k,2]+dc[3]*dprema[k,3]
	    }

	if (dcc[1] != 0.0d0) {
	    darg = dcc[2] / dcc[1]
	    dra2 = datan (darg)
	    if (dcc[1] < 0.0d0)
		dra2 = dra2 + dpi
	    else if (dcc[2] < 0.0d0)
		dra2 = dra2 + tpi
	    endif
	    }
	else {
	    if (dcc[2] > 0.0d0)
		dra2 = dpi/2.0d0
	    else
		dra2 = 1.5d0*dpi
	    }

	dec2 = dasin (dcc[3])

#  calculate hour angle = local sidereal time - r.a.

	dha = dst - dra2
	dha = dmod (dha + tpi , tpi)

#  calculate observer's geocentric velocity
#  (altitude assumed to be zero)

	dlats = dlat * dtr
	call geovel (dlats,dalt,dec2,-dha,gcvel)

#  calculate components of earth's barycentric veolcity,
#  dvelb(i), i=1,2,3  in units of a.u./s

	call barvel (djd,deqt,dvelh,dvelb)

#  project barycentric velocity to the direction of the star, and
#  convert to km/s

	bcvel = 0.0d0
	hcvel = 0.0d0
	do k=1,3 {
	    bcvel = bcvel + dvelb[k]*dcc[k]*daukm
	    hcvel = hcvel + dvelh[k]*dcc[k]*daukm
	    }

	return
end

#------------------------------------------------------------------------------
#
#       procedure:     sidereal (djd,dlong,dst)
#
#       purpose:	computes the mean local sidereal time.
#
#       input:	  djd   = julian date
#	               dlong = observer's longitude (radians)
#
#       output:	 dst = mean local sidereal time (radians)
#
#       note:	   constants taken from the american ephemeris
#	               and nautical almanac, 1980)
#
#       author:	 g. torres (1989)
#
#------------------------------------------------------------------------------

	procedure sidtim (djd,dlong,dst)

	double	djd	# julian date
	double	dlong	# longitude
	double	dst	# sidereal time (returned)

	double	dtpi, djd0, dst0, dut, dt

#   constants d1,d2,d3 for calculating greenwich mean sidereal time at 0 ut

	double	d1,d2,d3
	data d1/ 1.739935934667999d0  /
	data d2/ 6.283319509909095d02 /
	data d3/ 6.755878646261384d-06/

	double	dpi, df, dct0, dcjul

begin
	dpi = 3.141592653589793d0
	df = 1.00273790934d0
	dct0 = 2415020.0d0
	dcjul = 36525.0d0
	dtpi = 2.0d0 * dpi

	djd0 = dint (djd) + 0.5d0
	if (djd0 > djd)
	    djd0 = djd0 - 1.0d0
	dut = (djd - djd0) * dtpi

	dt = (djd0 - dct0)/dcjul
	dst0 = d1 + d2*dt + d3*dt*dt
	dst0 = mod (dst0,dtpi)
	dst = df*dut + dst0 - dlong
	dst = mod (dst + 2.0d0*dtpi , dtpi)

	return
end

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#
#       procedure:     pre (deq1,deq2,dprema)
#
#       purpose:	calculates the matrix of general precession from
#	               deq1 to deq2.
#
#       input:	  deq1 = initial epoch of mean equator and equinox
#	               deq2 = final epoch of mean equator and equinox
#
#       output:	 dprema = 3 x 3 matrix of general precession
#
#       note:	   the precession angles (dzeta,dzett,and dthet) are
#	               computed from the constants (dc1-dc9) corresponding
#	               to the definitions in the explanatory supplement
#	               to the american ephemeris (1961, p.30f).
#
#       author:	 p. stumpff (ibm-version 1979): astron. astrophys.
#	                suppl. ser. 41, 1 (1980)
#	               m. h. slovak (vax 11/780 implementation 1986)
#	               g. torres (1989)
#
#------------------------------------------------------------------------------

	procedure pre (deq1,deq2,dprema)

	double	dprema(3,3)

	double	dt0, dt, dts, dtc, dzeta, dzett, dthet, dszeta, dczeta
	double	dszett, dczett, dsthet, dcthet, da, db, dc, dd

	double	dcsar, dc1900, dc1m2, dc1, dc2, dc3, dc4, dc5, dc6
	double	dc7, dc8, dc9
	data dcsar/4.848136812d-6/, dc1900/1900.0d0/, dc1m2/0.01d0/,
     1	dc1/2304.25d0/, dc2/1.396d0/, dc3/0.302d0/, dc4/0.018d0/,
     2	dc5/0.791d0/, dc6/2004.683d0/, dc7/-0.853d0/, dc8/-0.426d0/,
     3	dc9/-0.042d0/

begin
	dt0 = (deq1 - dc1900)*dc1m2
	dt = (deq2 - deq1)*dc1m2
	dts = dt * dt
	dtc = dts * dt
	dzeta = ((dc1+dc2*dt0)*dt+dc3*dts+dc4*dtc)*dcsar
	dzett = dzeta + dc5*dts*dcsar
	dthet = ((dc6+dc7*dt0)*dt+dc8*dts+dc9*dtc)*dcsar
	dszeta = sin(dzeta)
	dczeta = cos(dzeta)
	dszett = sin(dzett)
	dczett = cos(dzett)
	dsthet = sin(dthet)
	dcthet = cos(dthet)
	da = dszeta * dszett
	db = dczeta * dszett
	dc = dszeta * dczett
	dd = dczeta * dczett

	dprema[1,1] = dd * dcthet - da
	dprema[1,2] = -1.d0 * dc * dcthet - db
	dprema[1,3] = -1.d0 * dsthet * dczett
	dprema[2,1] = db * dcthet + dc
	dprema[2,2] = -1.d0 * da * dcthet + dd
	dprema[2,3] = -1.d0 * dsthet * dszett
	dprema[3,1] = dczeta * dsthet
	dprema[3,2] = -1.d0 * dszeta * dsthet
	dprema[3,3] = dcthet

	return
end

#------------------------------------------------------------------------------ 
#
#
#       procedure:     geovel (dphi,dh,dec,dha,dvelg)
#
#       purpose:	calculates the correction required to transform
#	               the topocentric radial velocity of a given star
#	               to geocentric.
#	               - the maximum error of this routine is not expected
#	               to be larger than 0.1 m/s.
#
#       input:	  dphi = observer's geodetic latitude (radians)
#	               dh = observer's altitude above sea level (meters)
#
#	               dec = star's declination (radians) for mean
#	                      equator and equinox of date
#	               dha  = hour angle (radians)
#
#       output:	 dvelg = geocentric correction (km/s)
#
#       notes:	  vr = r.w.cos(dec).sin(hour angle), where r =
#	               geocentric radius at observer's latitude and
#	               altitude, and w = 2.pi/t, t = length of sidereal
#	               day (sec).  the hour angle is positive east of
#			the meridian.
#			other relevant equations from e. w. woolard
#	               & g. m. clemence (1966), spherical astronomy,
#			p.45 and p.48
#
#       author:	 g. torres (1989)
#
#------------------------------------------------------------------------------

procedure geovel(dphi,dh,dec,dha,dvelg)

double	dphi,dh,dec,dha,dvelg

double	de2,d1,d2,dr0,da,df,dw,dphig,drh

#   da = earth's equatorial radius (km)
#   df = polar flattening
#   dw = angular rotation rate (2.pi/t)

data da/6378.140d0/, df/0.00335281d0/, dw/7.2921158554d-05/

begin
	de2 = df*(2.0d0 - df)

#  calculate geocentric radius dr0 at sea level (km)

	d1 = 1.0d0 - de2*(2.0d0 - de2)*sin(dphi)**2
	d2 = 1.0d0 - de2*sin(dphi)**2
	dr0 = da * sqrt(d1/d2)

#  calculate geocentric latitude dphig

	d1 = de2*sin(2.0d0*dphi)
	d2 = 2.0d0*d2
	dphig = dphi - datan(d1/d2)

#  calculate geocentric radius drh at altitude dh (km)

	drh = dr0*cos(dphig) + dh/1.0d3*cos(dphi)

#  projected component to star at declination = dec and
#  at hour angle = dha, in units of km/s

	  dvelg = dw * drh * dcos(dec) * dsin(dha)

	return
end

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#
#       procedure:     barvel (dje,deq,dvelh,dvelb)
#
#       purpose:	calculates the heliocentric and barycentric
#	               velocity components of the earth.  the largest
#	               deviations from the jpl-de96 ephemeris are 42 cm/s
#	               for both heliocentric and barycentric velocity
#	               components.
#
#       input:	  dje = julian ephermeris date
#	               deq = epoch of mean equator and mean equinox of dvelh
#	                     and dvelb.  if deq = 0, both vectors are refer-
#	                     red to the mean equator and equinox of dje.
#
#       output:	 dvelh[k] = heliocentric velocity components
#	               dvelb[k] = barycentric velocity components
#	                (dx/dt, dy/dt, dz/dt, k=1,2,3  a.u./s)
#
#       author:	 p. stumpff  (ibm-version 1979): astron. astrophys.
#	                suppl. ser. 41, 1 (1980)
#	               m. h. slovak (vax 11/780 implementation 1986)
#	               g. torres (1989)
#
#------------------------------------------------------------------------------

procedure barvel (dje,deq,dvelh,dvelb)

double	dje,deq
double	dvelh[3]	# heliocentric velocity correction (returned)
double	dvelb[3]	# barycentric velocity correction (returned)

int	k, n
double	dcfel[3,8],dceps[3],dcargs[2,15]
double	dcargm[2,3], e, g
double	cc2pi,ccsec3,cckm,ccmld,ccfdi,t,tsq,a
double	dc1mme,ccsgd,dt,dtsq,dml,dlocal,pertld,pertr,pertrd
double	cosa,sina,esq,param,dparam,twoe,twog,f,sinf,cosf
double	dyhd,dzhd,b,tl,pertpd,pertp,drld,drd,dxhd,phid,psid
double	plon,pomg,pecc,dyahd,dzahd,dyabd,dzabd,deqdat
double	dcsld,dxbd,dybd,dzbd
double	ccamps[5,15],ccsec[3,4]
double	ccpamv[4],sn[4],ccsel[3,17],ccampm[4,3]

double	dc2pi,dc1,dct0,dcjul,dcbes,dctrop,dc1900,dtl,deps
double	phi,pertl

common /barxyz/ dprema,dpsi,d1pdro,dsinls,dcosls,dsinep,
     		dcosep,forbel,sorbel,sinlp,coslp,
     		sinlm,coslm,sigma,ideq
double	dprema[3,3],dpsi,d1pdro,dsinls,dcosls,dsinep,dcosep
double	forbel[7],sorbel[17],sinlp[4],coslp[4],sinlm,coslm,sigma
int	ideq

#  constants dcfel[i,k] of fast-changing elements
#	       i = 1             i = 2           i = 3
data dcfel/ 1.7400353d+00, 6.2833195099091d+02, 5.2796d-06,
	    6.2565836d+00, 6.2830194572674d+02,-2.6180d-06,
	    4.7199666d+00, 8.3997091449254d+03,-1.9780d-05,
	    1.9636505d-01, 8.4334662911720d+03,-5.6044d-05,
	    4.1547339d+00, 5.2993466764997d+01, 5.8845d-06,
	    4.6524223d+00, 2.1354275911213d+01, 5.6797d-06,
	    4.2620486d+00, 7.5025342197656d+00, 5.5317d-06,
	    1.4740694d+00, 3.8377331909193d+00, 5.6093d-06/

#   constants dceps and ccsel[i,k] of slowly changing elements
#	       i = 1        i = 2         i = 3
data dceps/ 4.093198d-01,-2.271110d-04,-2.860401d-08/
data ccsel/ 1.675104d-02,-4.179579d-05,-1.260516d-07,
	    2.220221d-01, 2.809917d-02, 1.852532d-05,
	    1.589963d+00, 3.418075d-02, 1.430200d-05,
	    2.994089d+00, 2.590824d-02, 4.155840d-06,
	    8.155457d-01, 2.486352d-02, 6.836840d-06,
	    1.735614d+00, 1.763719d-02, 6.370440d-06,
	    1.968564d+00, 1.524020d-02,-2.517152d-06,
	    1.282417d+00, 8.703393d-03, 2.289292d-05,
	    2.280820d+00, 1.918010d-02, 4.484520d-06,
	    4.833473d-02, 1.641773d-04,-4.654200d-07,
	    5.589232d-02,-3.455092d-04,-7.388560d-07,
	    4.634443d-02,-2.658234d-05, 7.757000d-08,
	    8.997041d-03, 6.329728d-06,-1.939256d-09,
	    2.284178d-02,-9.941590d-05, 6.787400d-08,
	    4.350267d-02,-6.839749d-05,-2.714956d-07,
	    1.348204d-02, 1.091504d-05, 6.903760d-07,
	    3.106570d-02,-1.665665d-04,-1.590188d-07/

#   constants of the arguments of the short-period perturbations by
#   the planets:  dcargs[i,k]
#	         i = 1             i = 2
data dcargs/ 5.0974222d+00,-7.8604195454652d+02,
	     3.9584962d+00,-5.7533848094674d+02,
	     1.6338070d+00,-1.1506769618935d+03,
	     2.5487111d+00,-3.9302097727326d+02,
	     4.9255514d+00,-5.8849265665348d+02,
	     1.3363463d+00,-5.5076098609303d+02,
	     1.6072053d+00,-5.2237501616674d+02,
	     1.3629480d+00,-1.1790629318198d+03,
	     5.5657014d+00,-1.0977134971135d+03,
	     5.0708205d+00,-1.5774000881978d+02,
	     3.9318944d+00, 5.2963464780000d+01,
	     4.8989497d+00, 3.9809289073258d+01,
	     1.3097446d+00, 7.7540959633708d+01,
	     3.5147141d+00, 7.9618578146517d+01,
	     3.5413158d+00,-5.4868336758022d+02/

#   amplitudes ccamps[n,k] of the short-period perturbations
#	  n = 1        n = 2        n = 3        n = 4        n = 5
data ccamps/
	-2.279594d-5, 1.407414d-5, 8.273188d-6, 1.340565d-5,-2.490817d-7,
	-3.494537d-5, 2.860401d-7, 1.289448d-7, 1.627237d-5,-1.823138d-7,
	 6.593466d-7, 1.322572d-5, 9.258695d-6,-4.674248d-7,-3.646275d-7,
	 1.140767d-5,-2.049792d-5,-4.747930d-6,-2.638763d-6,-1.245408d-7,
	 9.516893d-6,-2.748894d-6,-1.319381d-6,-4.549908d-6,-1.864821d-7,
	 7.310990d-6,-1.924710d-6,-8.772849d-7,-3.334143d-6,-1.745256d-7,
	-2.603449d-6, 7.359472d-6, 3.168357d-6, 1.119056d-6,-1.655307d-7,
	-3.228859d-6, 1.308997d-7, 1.013137d-7, 2.403899d-6,-3.736225d-7,
	 3.442177d-7, 2.671323d-6, 1.832858d-6,-2.394688d-7,-3.478444d-7,
	 8.702406d-6,-8.421214d-6,-1.372341d-6,-1.455234d-6,-4.998479d-8,
	-1.488378d-6,-1.251789d-5, 5.226868d-7,-2.049301d-7, 0.0d0,
	-8.043059d-6,-2.991300d-6, 1.473654d-7,-3.154542d-7, 0.0d0,
	 3.699128d-6,-3.316126d-6, 2.901257d-7, 3.407826d-7, 0.0d0,
	 2.550120d-6,-1.241123d-6, 9.901116d-8, 2.210482d-7, 0.0d0,
	-6.351059d-7, 2.341650d-6, 1.061492d-6, 2.878231d-7, 0.0d0/

#   constants of the secular perturbations in longitude ccsec3 and
#   ccsec[n,k]
#	       n = 1        n = 2          n = 3
data ccsec/  1.289600d-06, 5.550147d-01, 2.076942d+00,
	     3.102810d-05, 4.035027d+00, 3.525565d-01,
	     9.124190d-06, 9.990265d-01, 2.622706d+00,
	     9.793240d-07, 5.508259d+00, 1.559103d+01/

#   constants dcargm[i,k] of the arguments of the perturbations of the
#   motion of the moon
#	                 i = 1             i = 2
data dcargm/  5.1679830d+00, 8.3286911095275d+03,
	      5.4913150d+00,-7.2140632838100d+03,
	      5.9598530d+00, 1.5542754389685d+04/

#   amplitudes ccampm[n,k] of the perturbations of the moon
#	   n = 1         n = 2         n = 3         n = 4
data ccampm/
	 1.097594d-01, 2.896773d-07, 5.450474d-02, 1.438491d-07,
	-2.223581d-02, 5.083103d-08, 1.002548d-02,-2.291823d-08,
	 1.148966d-02, 5.658888d-08, 8.249439d-03, 4.063015d-08/

#  ccpamv = a*m*dl/dt (planets)
data ccpamv/8.326827d-11,1.843484d-11,1.988712d-12,1.881276d-12/

begin

	dc2pi = 6.2831853071796d0
	cc2pi = 6.283185d0
	dc1 = 1.0d0
	dct0 = 2415020.0d0
	dcjul = 36525.0d0
	dcbes = 0.313d0
	dctrop = 365.24219572d0
	dc1900 = 1900.0d0

#  sidereal rate dcsld in longitude, rate ccsgd in mean anomaly
	dcsld = 1.990987d-07
	ccsgd = 1.990969d-07
	ccsec3 = -7.757020d-08

#  dc1mme = 1 - mass(earth+moon)
	dc1mme = 0.99999696d0

#   some constants used in the calculation of the lunar contribution
	cckm = 3.122140d-05
	ccmld = 2.661699d-06
	ccfdi = 2.399485d-07

#  control-parameter ideq, and time-arguments

	ideq = deq
	dt = (dje - dct0) / dcjul
	t = dt
	dtsq = dt * dt
	tsq = dtsq

#  values of all elements for the instant dje

	do k = 1 , 8 {
	    dlocal = dmod (dcfel[1,k]+dt*dcfel[2,k]+dtsq*dcfel[3,k], dc2pi)
	    if (k == 1) dml = dlocal
	    if (k != 1) forbel[k-1] = dlocal
	    }

	deps = dmod (dceps[1]+dt*dceps[2]+dtsq*dceps[3], dc2pi)

	do k = 1 , 17 {
	    sorbel[k] = dmod (ccsel[1,k]+t*ccsel[2,k]+tsq*ccsel[3,k], dc2pi)
	    }

#  secular perturbations in longitude

	do k = 1 , 4 {
	    a = dmod (ccsec[2,k]+t*ccsec[3,k], dc2pi)
	    sn[k] = sin (a)
	    }

#  periodic perturbations of the emb (earth-moon barycenter)

	pertl = ccsec[1,1]	    *sn[1] +ccsec[1,2]*sn[2] +
		(ccsec[1,3]+t*ccsec3)*sn[3] +ccsec[1,4]*sn[4]

	pertld = 0.0
	pertr =  0.0
	pertrd = 0.0

	do k = 1 , 15 {
	    a = dmod (dcargs[1,k]+dt*dcargs[2,k], dc2pi)
	    cosa = cos(a)
	    sina = sin(a)
	    pertl = pertl + ccamps[1,k]*cosa + ccamps[2,k]*sina
	    pertr = pertr + ccamps[3,k]*cosa + ccamps[4,k]*sina
	    if (k < 11) {
		pertld = pertld + (ccamps[2,k]*cosa-ccamps[1,k]*sina)*ccamps[5,k]
		pertrd = pertrd + (ccamps[4,k]*cosa-ccamps[3,k]*sina)*ccamps[5,k]
		}
	    }

#  elliptic part of the motion of the emb

	e = sorbel[1]
	g = forbel[1]
	esq = e * e
	dparam = dc1 - esq
	param = dparam
	twoe = e + e
	twog = g + g
	phi = twoe*((1.0 - esq*(1./8.))*sin(g) + e*(5./8.)*sin(twog)
	      + esq*0.5416667*sin(g+twog) )
	f = g + phi
	sinf = sin(f)
	cosf = cos(f)
	dpsi = dparam/(dc1 + e*cosf)
	phid = twoe*ccsgd*((1.0+esq*1.5)*cosf+e*(1.25-sinf*sinf*0.5))
	psid = ccsgd*e*sinf/sqrt(param)

#  perturbed heliocentric motion of the emb

	d1pdro = (dc1 + pertr)
	drd = d1pdro*(psid+dpsi*pertrd)
	drld = d1pdro*dpsi*(dcsld+phid+pertld)
	dtl = dmod (dml+phi+pertl, dc2pi)
	dsinls = dsin(dtl)
	dcosls = dcos(dtl)
	dxhd = drd*dcosls - drld*dsinls
	dyhd = drd*dsinls + drld*dcosls

#  influence of eccentricity, evection and variation of the geocentric
#  motion of the moon

	pertl =  0.0
	pertld = 0.0
	pertp =  0.0
	pertpd = 0.0

	do k = 1, 3 {
	    a = dmod (dcargm[1,k] + dt*dcargm[2,k], dc2pi)
	    sina = sin(a)
	    cosa = cos(a)
	    pertl   = pertl  + ccampm[1,k]*sina
	    pertld  = pertld + ccampm[2,k]*cosa
	    pertp   = pertp  + ccampm[3,k]*cosa
	    pertpd  = pertpd - ccampm[4,k]*sina
	    }

#   heliocentric motion of the earth

	tl =  forbel[2] + pertl
	sinlm = sin(tl)
	coslm = cos(tl)
	sigma = cckm/(1.0 + pertp)
	a = sigma*(ccmld+pertld)
	b = sigma*pertpd
	dxhd = dxhd + a*sinlm + b*coslm
	dyhd = dyhd - a*coslm + b*sinlm
	dzhd =      - sigma*ccfdi*cos(forbel[3])

#   barycentric motion of the earth

	dxbd = dxhd*dc1mme
	dybd = dyhd*dc1mme
	dzbd = dzhd*dc1mme

	do k = 1 , 4 {
	    plon = forbel[k+3]
	    pomg = sorbel[k+1]
	    pecc = sorbel[k+9]
	    tl = dmod (plon+2.0*pecc*sin(plon-pomg), cc2pi)
	    sinlp[k] = sin(tl)
	    coslp[k] = cos(tl)
	    dxbd = dxbd + ccpamv[k]*(sinlp[k] + pecc*sin(pomg))
	    dybd = dybd - ccpamv[k]*(coslp[k] + pecc*cos(pomg))
	    dzbd = dzbd - ccpamv[k]*sorbel[k+13]*cos(plon-sorbel[k+5])
	    }

#   transition to mean equator of date

	  dcosep = dcos(deps)
	  dsinep = dsin(deps)
	  dyahd = dcosep*dyhd - dsinep*dzhd
	  dzahd = dsinep*dyhd + dcosep*dzhd
	  dyabd = dcosep*dybd - dsinep*dzbd
	  dzabd = dsinep*dybd + dcosep*dzbd
#
	if (ideq == 0) {
	    dvelh[1] = dxhd
	    dvelh[2] = dyahd
	    dvelh[3] = dzahd
	    dvelb[1] = dxbd
	    dvelb[2] = dyabd
	    dvelb[3] = dzabd

	    return
	    }

#   general precession from epoch dje to deq
	deqdat = (dje - dct0 - dcbes)/dctrop + dc1900

	call pre (deqdat,deq,dprema)

	do n = 1 , 3 {
	    dvelh[n] = dxhd*dprema[n,1] + dyahd*dprema[n,2] + 
		       dzahd*dprema[n,3]
	    dvelb[n] = dxbd*dprema[n,1] + dyabd*dprema[n,2] + 
		       dzabd*dprema[n,3]
	    }

	return
end

#        1990	return hcv and bcv

# Apr  1 1992	use altitude when computing velocity correction
# Aug 14 1992	return 0 from julday if year is <= 0

# Apr 13 1994	declare more variables
# Jun 15 1994	fix parameter descriptions
# Nov 16 1994	remove ^l's between procedures
# Nov 16 1994	shorten three 80-character lines

# Jul 13 1995	move julian date to spp procedure
