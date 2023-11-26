c*** File rvsao/Util/bcv.f
c*** January 30, 1997
c*** By G. Torres (1989)
c*** Modified by D. Mink

c  BCV calculates the correction required to reduce observed (topocentric)
c  radial velocities of a given star to the barycenter of the solar system.
c  It includes correction for the effect of the earth's rotation.
c  The maximum error of this routine is not expected to be larger than 0.6 m/sec

	Subroutine BCV (DJD,DLONG,DLAT,DALT,DRA,DEC,DEQ,BCVEL,HCVEL,GCVEL)

	Real*8 DJD
c				Heliocentric Julian date (days)
	Real*8 DLONG
c				Geodetic longitude (degrees, west positive)
	Real*8 DLAT
c				Geodetic latitude (degrees)
	Real*8 DALT
c				Altitude above sea level (meters)
	Real*8 DRA
c				Right ascension of star (hours)
	Real*8 DEC
c				Declination of star (degrees)
	Real*8 DEQ
c				Mean equator and equinox for coordinates
c				e.g., 1950.0
	Real*8 BCVEL
c				Barycentric correction (km/s) (returned)
	Real*8 HCVEL
c				Heliocentric correction (km/s) (returned)
	Real*8 GCVEL
c				Geocentric correction (km/s) (returned)

	Real*8 DLONGS
c				Geodetic longitude (radians, west+)
	Real*8 DLATS
c				Geodetic latitude (radians)
	Real*8 DRAS
c				Right ascension of star (radians)
	Real*8 DECS
c				Declination of star (radians)

	Real*8 DC(3),DCC(3),DPREMA(3,3),DVELH(3),DVELB(3)
	Real*8 DPI,DAUKM,DCTROP,DCBES,DC1900,DCT0,DTR,DST
	Real*8 DEQT,DRA2,DEC2,DHA,DARG
	Integer*4 K
	Data DPI/3.1415926535897932d0/, DAUKM/1.4959787d08/,
     1	     DCTROP/365.24219572d0/, DCBES/0.313d0/, DC1900/1900.0d0/
     2	     DCT0/2415020.0d0/

	DTR = DPI / 180.0d0

	DLONGS = DLONG*DTR
	DLATS = DLAT * DTR
	DRAS = DRA * 15.0d0 * DTR
	DECS = DEC * DTR
c	Open (15,FILE='testbcv')
c	Write(15,*) 'Julian Date:',DJD
c	Write(15,*) 'Long:',DLONGS,' Lat:',DLATS,DALT
c	Write(15,*) 'RA:',DRAS,' Dec:',DECS,DEQ

c  Calculate local sidereal time

	Call SIDTIM (DJD,DLONGS,DST)

c  Precess R.A. and Dec. to mean equator and equinox of date (DEQT)

	DEQT = (DJD - DCT0 - DCBES)/DCTROP + DC1900
	DC(1) = Dcos(DRAS) * Dcos(DECS)
	DC(2) = Dsin(DRAS) * Dcos(DECS)
	DC(3) =	      Dsin(DECS)

	Call PRE (DEQ,DEQT,DPREMA)
	Do 100 K=1,3
	    DCC(K)=DC(1)*DPREMA(K,1)+DC(2)*DPREMA(K,2)+DC(3)*DPREMA(K,3)
100	    Continue

	If (DCC(1) .ne. 0.0d0) Then
	    DARG = DCC(2) / DCC(1)
	    DRA2 = Datan (DARG)
	    If (DCC(1) .lt. 0.0d0) Then
		DRA2 = DRA2 + DPI
	    Elseif (DCC(2) .lt. 0.0d0) Then
		DRA2 = DRA2 + 2.0d0*DPI
	    Endif
	Else
	    If (DCC(2) .gt. 0.0d0) Then
		DRA2 = DPI/2.0d0
	    Else
		DRA2 = 1.5d0*DPI
	    Endif
	Endif

	DEC2 = DASIN (DCC(3))

c  Calculate hour angle = local sidereal time - R.A.

	DHA = DST - DRA2
c	DHA = Dmod (DHA + 2.0d0*DPI , 2.0d0*DPI)

c	Write(15,*) 'RA2=',DRA2,'  DEC2=',DEC2,'  ALT=',DALT
c	Write(15,*) 'ST=',DST,'  HA=',DHA,'  LAT:',DLATS

c  Calculate observer's geocentric velocity
c  (altitude assumed to be zero)

	Call GEOVEL (DLATS,DALT,DEC2,-DHA,GCVEL)

c  Calculate components of earth's barycentric velocity,
c  DVELB(I), I=1,2,3  in units of A.U./S

	Call BARVEL (DJD,DEQT,DVELH,DVELB)

c	Write(15,*) 'BCV: DJD=',DJD,'  DEQT=',DEQT
c	Write(15,*) 'BCV: DVELH=',DVELH
c	Write(15,*) 'BCV: DVELB=',DVELB

c  Project barycentric velocity to the direction of the star, and
c  convert to km/s

	BCVEL = 0.0d0
	HCVEL = 0.0d0
	Do 200 K=1,3
	    BCVEL = BCVEL + DVELB(K)*DCC(K)*DAUKM
	    HCVEL = HCVEL + DVELH(K)*DCC(K)*DAUKM
200	    Continue

c	Close (15)
	Return

	End
c------------------------------------------------------------------------------
c
c       Subroutine:     sidereal (DJD,DLONG,DST)
c
c       PURPOSE:	COMPUTES THE MEAN LOCAL sidereal TIME.
c
c       INPUT:	  DJD   = JULIAN DATE
c		       DLONG = OBSERVER'S LONGITUDE (RADIANS)
c
c       OUTPUT:	 DST = MEAN LOCAL sidereal TIME (RADIANS)
c
c       NOTE:	   CONSTANTS TAKEN FROM THE AMERICAN EPHEMERIS
c		       AND NAUTICAL ALMANAC, 1980)
c
c       AUTHOR:	 G. TORRES (1989)
c
c------------------------------------------------------------------------------

	Subroutine SIDTIM (DJD,DLONG,DST)

	Real*8 DJD
c			Julian Date
	Real*8 DLONG
c			Longitude
	Real*8 DST
c			Sidereal TIme (returned)

	Real*8 DTPI, DJD0, DST0, DUT, DT

	Real*8 D1,D2,D3

	Real*8 DPI, DF, DCT0, DCJUL

	DPI = 3.141592653589793d0
	DTPI = 2.0d0*DPI
	DF = 1.00273790934d0
	DCT0 = 2415020.0d0
	DCJUL = 36525.0d0

c   Constants D1,D2,D3 for calculating Greenwich Mean Sidereal Time at 0 UT
	D1 = 1.739935934667999d0
	D2 = 6.283319509909095d02
	D3 = 6.755878646261384d-06

	DJD0 = IDINT(DJD) + 0.5d0
	If (DJD0.GT.DJD)  DJD0 = DJD0 - 1.0d0
	DUT = (DJD - DJD0)*DTPI

	DT = (DJD0 - DCT0)/DCJUL
	DST0 = D1 + D2*DT + D3*DT*DT
	DST0 = Dmod (DST0,DTPI)
	DST = DF*DUT + DST0 - DLONG
	DST = Dmod (DST + 2.0d0*DTPI , DTPI)

	Return

	End
c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
c
c       Subroutine:     PRE (DEQ1,DEQ2,DPREMA)
c
c       PURPOSE:	CALCULATES THE MATRIX OF GENERAL PRECESSION FROM
c		       DEQ1 TO DEQ2.
c
c       INPUT:	  DEQ1 = INITIAL EPOCH OF MEAN EQUATOR AND EQUINOX
c		       DEQ2 = FINAL EPOCH OF MEAN EQUATOR AND EQUINOX
c
c       OUTPUT:	 DPREMA = 3 X 3 MATRIX OF GENERAL PRECESSION
c
c       NOTE:	   THE PRECESSION ANGLES (DZETA,DZETT,AND DTHET) ARE
c		       COMPUTED FROM THE CONSTANTS (DC1-DC9) CORRESPONDING
c		       TO THE DEFINITIONS IN THE EXPLANATORY SUPPLEMENT
c		       TO THE AMERICAN EPHEMERIS (1961, P.30F).
c
c       AUTHOR:	 P. STUMPFF (IBM-VERSION 1979): ASTRON. ASTROPHYS.
c			SUPPL. SER. 41, 1 (1980)
c		       M. H. SLOVAK (VAX 11/780 IMPLEMENTATION 1986)
c		       G. TORRES (1989)
c
c------------------------------------------------------------------------------

	Subroutine PRE (DEQ1,DEQ2,DPREMA)

	Real*8 DEQ1, DEQ2
	Real*8 DPREMA(3,3)

	Real*8 DT0, DT, DTS, DTC, DZETA, DZETT, DTHET, DSZETA, DCZETA
	Real*8 DSZETT, DCZETT, DSTHET, DCTHET, DA, DB, DC, DD

	Real*8 DCSAR, DC1900, DC1M2, DC1, DC2, DC3, DC4, DC5, DC6
	Real*8 DC7, DC8, DC9

	DCSAR = 4.848136812d-6
	DC1900 = 1900.0d0
	DC1M2 = 0.01d0
	DC1 = 2304.25d0
	DC2 = 1.396d0
	DC3 = 0.302d0
	DC4 = 0.018d0
	DC5 = 0.791d0
	DC6 = 2004.683d0
	DC7 = -0.853d0
	DC8 = -0.426d0
	DC9 = -0.042d0

	DT0 = (DEQ1 - DC1900)*DC1M2
	DT = (DEQ2 - DEQ1)*DC1M2
	DTS = DT * DT
	DTC = DTS * DT
	DZETA = ((DC1+DC2*DT0)*DT+DC3*DTS+DC4*DTC)*DCSAR
	DZETT = DZETA + DC5*DTS*DCSAR
	DTHET = ((DC6+DC7*DT0)*DT+DC8*DTS+DC9*DTC)*DCSAR
	DSZETA = Dsin(DZETA)
	DCZETA = Dcos(DZETA)
	DSZETT = Dsin(DZETT)
	DCZETT = Dcos(DZETT)
	DSTHET = Dsin(DTHET)
	DCTHET = Dcos(DTHET)
	DA = DSZETA * DSZETT
	DB = DCZETA * DSZETT
	DC = DSZETA * DCZETT
	DD = DCZETA * DCZETT

c	Write(15,*) 'DZETA=',DZETA,'  DZETT=',DZETT,'  DTHET=',DTHET
c	Write(15,*) 'CZETA=',DCZETA,'  SZETA=',DSZETA
c	Write(15,*) 'CZETT=',DCZETT,'  SZETT=',DSZETT
c	Write(15,*) 'CTHET=',DCTHET,'  STHET=',DSTHET

	DPREMA(1,1) = DD * DCTHET - DA
	DPREMA(1,2) = -1.d0 * DC * DCTHET - DB
	DPREMA(1,3) = -1.d0 * DSTHET * DCZETT
	DPREMA(2,1) = DB * DCTHET + DC
	DPREMA(2,2) = -1.d0 * DA * DCTHET + DD
	DPREMA(2,3) = -1.d0 * DSTHET * DSZETT
	DPREMA(3,1) = DCZETA * DSTHET
	DPREMA(3,2) = -1.d0 * DSZETA * DSTHET
	DPREMA(3,3) = DCTHET

	Return

	End

c------------------------------------------------------------------------------ 
c
c
c       Subroutine:     GEOVEL (DPHI,DH,DEC,DHA,DVELG)
c
c       PURPOSE:	CALCULATES THE CORRECTION REQUIRED TO TRANSFORM
c		       THE TOPOCENTRIC RADIAL VELOCITY OF A GIVEN STAR
c		       TO GEOCENTRIC.
c		       - THE MAXIMUM ERROR OF THIS ROUTINE IS NOT EXPECTED
c		       TO BE LARGER THAN 0.1 M/S.
c
c       INPUT:	  DPHI = OBSERVER'S GEODETIC LATITUDE (RADIANS)
c		       DH = OBSERVER'S ALTITUDE ABOVE SEA LEVEL (METERS)
c
c		       DEC = STAR'S DECLINATION (RADIANS) FOR MEAN
c			      EQUATOR AND EQUINOX OF DATE
c		       DHA  = HOUR ANGLE (RADIANS)
c
c       OUTPUT:	 DVELG = GEOCENTRIC CORRECTION (KM/S)
c
c       NOTES:	  VR = R.W.COS(DEC).SIN(HOUR ANGLE), WHERE R =
c		       GEOCENTRIC RADIUS AT OBSERVER'S LATITUDE AND
c		       ALTITUDE, AND W = 2.PI/T, T = LENGTH OF sidereal
c		       DAY (SEC).  THE HOUR ANGLE IS POSITIVE EAST OF
c			THE MERIDIAN.
c			OTHER RELEVANT EQUATIONS FROM E. W. WOOLARD
c		       & G. M. CLEMENCE (1966), SPHERICAL ASTRONOMY,
c			P.45 AND P.48
c
c       AUTHOR:	 G. TORRES (1989)
c
c------------------------------------------------------------------------------

	Subroutine GEOVEL(DPHI,DH,DEC,DHA,DVELG)

	Real*8 DPHI,DH,DEC,DHA,DVELG

	Real*8 DE2,D1,D2,DR0,DA,DF,DW,DPHIG,DRH

c  Earth's equatorial radius (KM)
	DA = 6378.140d0

c  Polar flattening
	DF = 0.00335281d0

c  Angular rotation rate (2.PI/T)
	DW = 7.2921158554d-05

	DE2 = DF*(2.0d0 - DF)

c  Calculate geocentric radius DR0 at sea level (KM)
	D1 = 1.0d0 - DE2*(2.0d0 - DE2)*Dsin(DPHI)**2
	D2 = 1.0d0 - DE2*Dsin(DPHI)**2
	DR0 = DA * DSQRT(D1/D2)

c  Calculate geocentric latitude DPHIG
	D1 = DE2*Dsin(2.0d0*DPHI)
	D2 = 2.0d0*D2
	DPHIG = DPHI - DataN(D1/D2)

c  Calculate geocentric radius DRH at altitude DH (KM)
	DRH = DR0*Dcos(DPHIG) + DH/1.0D3*Dcos(DPHI)

c  Projected component to star at declination = DEC and
c  at hour angle = DHA, in units of km/s
	DVELG = DW * DRH * Dcos(DEC) * Dsin(DHA)

	Return

	End
c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
c
c       Subroutine:     BARVEL (DJE,DEQ,DVELH,DVELB)
c
c       PURPOSE:	CALCULATES THE HELIOCENTRIC AND BARYCENTRIC
c		       VELOCITY COMPONENTS OF THE EARTH.  THE LARGEST
c		       DEVIATIONS FROM THE JPL-DE96 EPHEMERIS ARE 42 CM/S
c		       FOR BOTH HELIOCENTRIC AND BARYCENTRIC VELOCITY
c		       COMPONENTS.
c
c       INPUT:	  DJE = JULIAN EPHERMERIS DATE
c		       DEQ = EPOCH OF MEAN EQUATOR AND MEAN EQUINOX OF DVELH
c			     AND DVELB.  If DEQ = 0, BOTH VECTORS ARE REFER-
c			     RED TO THE MEAN EQUATOR AND EQUINOX OF DJE.
c
c       OUTPUT:	 DVELH(K) = HELIOCENTRIC VELOCITY COMPONENTS
c		       DVELB(K) = BARYCENTRIC VELOCITY COMPONENTS
c			(DX/DT, DY/DT, DZ/DT, K=1,2,3  A.U./S)
c
c       AUTHOR:	 P. STUMPFF  (IBM-VERSION 1979): ASTRON. ASTROPHYS.
c			SUPPL. SER. 41, 1 (1980)
c		       M. H. SLOVAK (VAX 11/780 IMPLEMENTATION 1986)
c		       G. TORRES (1989)
c
c------------------------------------------------------------------------------

	Subroutine BARVEL (DJE,DEQ,DVELH,DVELB)

	Real*8 DJE,DEQ
	Real*8 DVELH(3)
c				Heliocentric velocity correction (returned)
	Real*8 DVELB(3)
c				Barycentric velocity correction (returned)

	Integer*4 K, N
	Real*8 DCFEL(3,8),DCEPS(3),DCARGS(2,15)
	Real*8 DCARGM(2,3), E, G
	Real*8 CC2PI,CCSEC3,CCKM,CCMLD,CCFDI,T,TSQ,A
	Real*8 DC1MME,CCSGD,DT,DTSQ,DML,DLOCAL,PERTLD,PERTR,PERTRD
	Real*8 COSA,SINA,ESQ,PARAM,DPARAM,TWOE,TWOG,F,SINF,COSF
	Real*8 DYHD,DZHD,B,TL,PERTPD,PERTP,DRLD,DRD,DXHD,PHID,PSID
	Real*8 PLON,POMG,PECC,DYAHD,DZAHD,DYABD,DZABD,DEQDAT
	Real*8 DCSLD,DXBD,DYBD,DZBD
	Real*8 CCAMPS(5,15),CCSEC(3,4)
	Real*8 CCPAMV(4),SN(4),CCSEL(3,17),CCAMPM(4,3)

	Common /BARXYZ/ DPREMA(3,3),DPSI,D1PDRO,DSINLS,DCOSLS,DSINEP,
     1			DCOSEP,FORBEL(7),SORBEL(17),SINLP(4),COSLP(4),
     2			SINLM,COSLM,SIGMA,IDEQ
	Real*8 DPREMA,DPSI,D1PDRO,DSINLS,DCOSLS,DSINEP,DCOSEP
	Real*8 FORBEL,SORBEL,SINLP,COSLP,SINLM,COSLM,SIGMA
	Integer*4 IDEQ

	Equivalence (SORBEL(1),E),(FORBEL(1),G)

	Real*8 DC2PI,DC1,DCT0,DCJUL,DCBES,DCTROP,DC1900,DTL,DEPS
	Real*8 PHI,PERTL

	Data DC2PI/6.2831853071796d0/, CC2PI/6.283185/,
     1	     DC1  /1.0d0/ , DCT0/2415020.0d0/, DCJUL/36525.0d0/,
     2	     DCBES/0.313d0/, DCTROP/365.24219572d0/, DC1900/1900.0d0/

c  Constants DCFEL(I,K) of fast-changing elements

c		       I = 1	     I = 2	   I = 3

	Data DCFEL/ 1.7400353d+00, 6.2833195099091d+02, 5.2796d-06,
     1	      6.2565836d+00, 6.2830194572674d+02,-2.6180d-06,
     1	      4.7199666d+00, 8.3997091449254d+03,-1.9780d-05,
     1	      1.9636505d-01, 8.4334662911720d+03,-5.6044d-05,
     1	      4.1547339d+00, 5.2993466764997d+01, 5.8845d-06,
     1	      4.6524223d+00, 2.1354275911213d+01, 5.6797d-06,
     1	      4.2620486d+00, 7.5025342197656d+00, 5.5317d-06,
     1	      1.4740694d+00, 3.8377331909193d+00, 5.6093d-06/

c   CONSTANTS DCEPS AND CCSEL(I,K) OF SLOWLY CHANGING ELEMENTS

c		       I = 1	I = 2	 I = 3

	Data DCEPS/ 4.093198d-01,-2.271110d-04,-2.860401d-08/

	Data CCSEL/ 1.675104d-02,-4.179579d-05,-1.260516d-07,
     1	      2.220221d-01, 2.809917d-02, 1.852532d-05,
     1	      1.589963d+00, 3.418075d-02, 1.430200d-05,
     1	      2.994089d+00, 2.590824d-02, 4.155840d-06,
     1	      8.155457d-01, 2.486352d-02, 6.836840d-06,
     1	      1.735614d+00, 1.763719d-02, 6.370440d-06,
     1	      1.968564d+00, 1.524020d-02,-2.517152d-06,
     1	      1.282417d+00, 8.703393d-03, 2.289292d-05,
     1	      2.280820d+00, 1.918010d-02, 4.484520d-06,
     1	      4.833473d-02, 1.641773d-04,-4.654200d-07,
     1	      5.589232d-02,-3.455092d-04,-7.388560d-07,
     1	      4.634443d-02,-2.658234d-05, 7.757000d-08,
     1	      8.997041d-03, 6.329728d-06,-1.939256d-09,
     1	      2.284178d-02,-9.941590d-05, 6.787400d-08,
     1	      4.350267d-02,-6.839749d-05,-2.714956d-07,
     1	      1.348204d-02, 1.091504d-05, 6.903760d-07,
     1	      3.106570d-02,-1.665665d-04,-1.590188d-07/

c   CONSTANTS OF THE ARGUMENTS OF THE SHORT-PERIOD PERTURBATIONS BY
c   THE PLANETS:  DCARGS(I,K)

c			I = 1	     I = 2

	Data DCARGS/ 5.0974222d+00,-7.8604195454652d+02,
     1	       3.9584962d+00,-5.7533848094674d+02,
     1	       1.6338070d+00,-1.1506769618935d+03,
     1	       2.5487111d+00,-3.9302097727326d+02,
     1	       4.9255514d+00,-5.8849265665348d+02,
     1	       1.3363463d+00,-5.5076098609303d+02,
     1	       1.6072053d+00,-5.2237501616674d+02,
     1	       1.3629480d+00,-1.1790629318198d+03,
     1	       5.5657014d+00,-1.0977134971135d+03,
     1	       5.0708205d+00,-1.5774000881978d+02,
     1	       3.9318944d+00, 5.2963464780000d+01,
     1	       4.8989497d+00, 3.9809289073258d+01,
     1	       1.3097446d+00, 7.7540959633708d+01,
     1	       3.5147141d+00, 7.9618578146517d+01,
     1	       3.5413158d+00,-5.4868336758022d+02/

c   AMPLITUDES CCAMPS(N,K) OF THE SHORT-PERIOD PERTURBATIONS

c	  N = 1	N = 2	N = 3	N = 4	N = 5

	Data CCAMPS/
     1 -2.279594d-5, 1.407414d-5, 8.273188d-6, 1.340565d-5,-2.490817d-7,
     1 -3.494537d-5, 2.860401d-7, 1.289448d-7, 1.627237d-5,-1.823138d-7,
     1  6.593466d-7, 1.322572d-5, 9.258695d-6,-4.674248d-7,-3.646275d-7,
     1  1.140767d-5,-2.049792d-5,-4.747930d-6,-2.638763d-6,-1.245408d-7,
     1  9.516893d-6,-2.748894d-6,-1.319381d-6,-4.549908d-6,-1.864821d-7,
     1  7.310990d-6,-1.924710d-6,-8.772849d-7,-3.334143d-6,-1.745256d-7,
     1 -2.603449d-6, 7.359472d-6, 3.168357d-6, 1.119056d-6,-1.655307d-7,
     1 -3.228859d-6, 1.308997d-7, 1.013137d-7, 2.403899d-6,-3.736225d-7,
     1  3.442177d-7, 2.671323d-6, 1.832858d-6,-2.394688d-7,-3.478444d-7,
     1  8.702406d-6,-8.421214d-6,-1.372341d-6,-1.455234d-6,-4.998479d-8,
     1 -1.488378d-6,-1.251789d-5, 5.226868d-7,-2.049301d-7, 0.0d0,
     1 -8.043059d-6,-2.991300d-6, 1.473654d-7,-3.154542d-7, 0.0d0,
     1  3.699128d-6,-3.316126d-6, 2.901257d-7, 3.407826d-7, 0.0d0,
     1  2.550120d-6,-1.241123d-6, 9.901116d-8, 2.210482d-7, 0.0d0,
     1 -6.351059d-7, 2.341650d-6, 1.061492d-6, 2.878231d-7, 0.0d0/

c   CONSTANTS OF THE SECULAR PERTURBATIONS IN LONGITUDE CCSEC3 AND
c   CCSEC(N,K)

c			N = 1	N = 2	  N = 3

	Data CCSEC3/-7.757020d-08/

	Data CCSEC/  1.289600d-06, 5.550147d-01, 2.076942d+00,
     1	       3.102810d-05, 4.035027d+00, 3.525565d-01,
     1	       9.124190d-06, 9.990265d-01, 2.622706d+00,
     1	       9.793240d-07, 5.508259d+00, 1.559103d+01/

c   Sidereal RATE DCSLD IN LONGITUDE, RATE CCSGD IN MEAN ANOMALY

	Data DCSLD/ 1.990987d-07/, CCSGD/ 1.990969d-07/

c   SOME CONSTANTS USED IN THE CALCULATION OF THE LUNAR CONTRIBUTION

	Data CCKM/3.122140d-05/, CCMLD/2.661699d-06/, CCFDI/2.399485d-07/

c   CONSTANTS DCARGM(I,K) OF THE ARGUMENTS OF THE PERTURBATIONS OF THE
c   MOTION OF THE MOON

c			 I = 1	     I = 2

	Data DCARGM/  5.1679830d+00, 8.3286911095275d+03,
     1		5.4913150d+00,-7.2140632838100d+03,
     1		5.9598530d+00, 1.5542754389685d+04/

c   AMPLITUDES CCAMPM(N,K) OF THE PERTURBATIONS OF THE MOON

c	   N = 1	 N = 2	 N = 3	 N = 4

	Data CCAMPM/
     1	 1.097594d-01, 2.896773d-07, 5.450474d-02, 1.438491d-07,
     1	-2.223581d-02, 5.083103d-08, 1.002548d-02,-2.291823d-08,
     1	 1.148966d-02, 5.658888d-08, 8.249439d-03, 4.063015d-08/

c  CCPAMV = A*M*DL/DT (PLANETS); DC1MME = 1 - MASS(EARTH+MOON)

	Data CCPAMV/8.326827d-11,1.843484d-11,1.988712d-12,1.881276d-12/,
     1	     DC1MME/0.99999696d0/

c  Program execution begins

c  Control-parameter IDEQ, and time-arguments

	IDEQ = DEQ
	DT = (DJE - DCT0)/DCJUL
	T = DT
	DTSQ = DT * DT
	TSQ = DTSQ

c  Values of all elements for the instant DJE

	Do 100 K = 1 , 8
	    DLOCAL = Dmod (DCFEL(1,K)+DT*DCFEL(2,K)+DTSQ*DCFEL(3,K), DC2PI)
	    If(K.EQ.1) DML = DLOCAL
	    If(K .ne. 1) FORBEL(K-1) = DLOCAL
100	    Continue

	DEPS = Dmod (DCEPS(1)+DT*DCEPS(2)+DTSQ*DCEPS(3), DC2PI)

	Do 200 K = 1 , 17
	    SORBEL(K) = Dmod (CCSEL(1,K)+T*CCSEL(2,K)+TSQ*CCSEL(3,K), DC2PI)
200	    Continue

c  Secular perturbations in longitude

	Do 300 K = 1 , 4
	    A = Dmod (CCSEC(2,K)+T*CCSEC(3,K), DC2PI)
	    SN(K) = SIN(A)
300	    Continue

c  PERIODIC PERTURBATIONS OF THE EMB (EARTH-MOON BARYCENTER)

	PERTL = CCSEC(1,1)	    *SN(1) +CCSEC(1,2)*SN(2)
     1		+(CCSEC(1,3)+T*CCSEC3)*SN(3) +CCSEC(1,4)*SN(4)

	PERTLD = 0.0
	PERTR =  0.0
	PERTRD = 0.0

	DO 400 K = 1 , 15
	    A = Dmod (DCARGS(1,K)+DT*DCARGS(2,K), DC2PI)
	    COSA = COS(A)
	    SINA = SIN(A)
	    PERTL = PERTL + CCAMPS(1,K)*COSA + CCAMPS(2,K)*SINA
	    PERTR = PERTR + CCAMPS(3,K)*COSA + CCAMPS(4,K)*SINA
	    If(K.GE.11) Go to 400
	    PERTLD = PERTLD + (CCAMPS(2,K)*COSA-CCAMPS(1,K)*SINA)*CCAMPS(5,K)
	    PERTRD = PERTRD + (CCAMPS(4,K)*COSA-CCAMPS(3,K)*SINA)*CCAMPS(5,K)
400	    Continue

c  ELLIPTIC PART OF THE MOTION OF THE EMB

	ESQ = E * E
	DPARAM = DC1 - ESQ
	PARAM = DPARAM
	TWOE = E + E
	TWOG = G + G
	PHI = TWOE*((1.d0 - ESQ*(1.d0/8.d0))*SIN(G) + E*(5.d0/8.d0)*SIN(TWOG)
     1	      + ESQ*0.5416667d0*SIN(G+TWOG) )
	F = G + PHI
	SINF = SIN(F)
	COSF = COS(F)
	DPSI = DPARAM/(DC1 + E*COSF)
	PHID = TWOE*CCSGD*((1.d0+ESQ*1.5d0)*COSF+E*(1.25d0-SINF*SINF*0.5d0))
	PSID = CCSGD*E*SINF/SQRT(PARAM)

c  PERTURBED HELIOCENTRIC MOTION OF THE EMB

	D1PDRO = (DC1 + PERTR)
	DRD = D1PDRO*(PSID+DPSI*PERTRD)
	DRLD = D1PDRO*DPSI*(DCSLD+PHID+PERTLD)
	DTL = Dmod (DML+PHI+PERTL, DC2PI)
	DSINLS = Dsin(DTL)
	DCOSLS = Dcos(DTL)
	DXHD = DRD*DCOSLS - DRLD*DSINLS
	DYHD = DRD*DSINLS + DRLD*DCOSLS

c  INFLUENCE OF ECCENTRICITY, EVECTION AND VARIATION OF THE GEOCENTRIC
c  MOTION OF THE MOON

	PERTL =  0.d0
	PERTLD = 0.d0
	PERTP =  0.d0
	PERTPD = 0.d0

	Do 500 K = 1, 3
	    A = Dmod (DCARGM(1,K) + DT*DCARGM(2,K), DC2PI)
	    SINA = SIN(A)
	    COSA = COS(A)
	    PERTL   = PERTL  + CCAMPM(1,K)*SINA
	    PERTLD  = PERTLD + CCAMPM(2,K)*COSA
	    PERTP   = PERTP  + CCAMPM(3,K)*COSA
	    PERTPD  = PERTPD - CCAMPM(4,K)*SINA
500	    Continue

c   HELIOCENTRIC MOTION OF THE EARTH

	TL =  FORBEL(2) + PERTL
	SINLM = SIN(TL)
	COSLM = COS(TL)
	SIGMA = CCKM/(1.0 + PERTP)
	A = SIGMA*(CCMLD+PERTLD)
	B = SIGMA*PERTPD
	DXHD = DXHD + A*SINLM + B*COSLM
	DYHD = DYHD - A*COSLM + B*SINLM
	DZHD =      - SIGMA*CCFDI*COS(FORBEL(3))

c   BARYCENTRIC MOTION OF THE EARTH

	DXBD = DXHD*DC1MME
	DYBD = DYHD*DC1MME
	DZBD = DZHD*DC1MME

	Do 600 K = 1 , 4
	    PLON = FORBEL(K+3)
	    POMG = SORBEL(K+1)
	    PECC = SORBEL(K+9)
	    TL = Dmod (PLON+2.0*PECC*SIN(PLON-POMG), CC2PI)
	    SINLP(K) = SIN(TL)
	    COSLP(K) = COS(TL)
	    DXBD = DXBD + CCPAMV(K)*(SINLP(K) + PECC*SIN(POMG))
	    DYBD = DYBD - CCPAMV(K)*(COSLP(K) + PECC*COS(POMG))
	    DZBD = DZBD - CCPAMV(K)*SORBEL(K+13)*COS(PLON-SORBEL(K+5))
600	    Continue

c   TRANSITION TO MEAN EQUATOR OF DATE

	  DCOSEP = Dcos(DEPS)
	  DSINEP = Dsin(DEPS)
	  DYAHD = DCOSEP*DYHD - DSINEP*DZHD
	  DZAHD = DSINEP*DYHD + DCOSEP*DZHD
	  DYABD = DCOSEP*DYBD - DSINEP*DZBD
	  DZABD = DSINEP*DYBD + DCOSEP*DZBD
c
	  If (IDEQ .ne. 0)  Go to 700
	  DVELH(1) = DXHD
	  DVELH(2) = DYAHD
	  DVELH(3) = DZAHD
	  DVELB(1) = DXBD
	  DVELB(2) = DYABD
	  DVELB(3) = DZABD

	  Return

c   GENERAL PRECESSION FROM EPOCH DJE TO DEQ

700	DEQDAT = (DJE - DCT0 - DCBES)/DCTROP + DC1900

	Call PRE(DEQDAT,DEQ,DPREMA)

	Do 800 N = 1 , 3
	    DVELH(N) = DXHD*DPREMA(N,1) + DYAHD*DPREMA(N,2) + 
     1		       DZAHD*DPREMA(N,3)
	    DVELB(N) = DXBD*DPREMA(N,1) + DYABD*DPREMA(N,2) + 
     1		       DZABD*DPREMA(N,3)
800	    Continue

	Return

	End

c	1990	Return HCV AND BCV

c Apr  1 1992	Use altitude when computing velocity correction
c Aug 14 1992	Return 0 from JULDAY if year is <= 0

c Apr 13 1994	Declare more variables
c Jun 15 1994	Fix parameter descriptions
c Nov 16 1994	Remove ^L's between subroutines
c Nov 16 1994	Shorten three 80-character lines

c Jul 13 1995	Move Julian Date to SPP subroutine

c Jan  3 1997	Note that Julian Date is heliocentric
c Jan 24 1997	Make all constants explicitly double precision
c Jan 30 1997	Fix to match RVCORRECT
