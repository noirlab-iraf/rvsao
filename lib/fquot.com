# Parameters used by FQUOT only

common /fixq/ cpf,chi0,unc,niter

real	cpf		# Detector units per detected photon
real	chi0		# Initialization of chi-square
real	unc		# Convergence criterion for chi-square minimization
int	niter		# Max iterations for chi-square minimization

