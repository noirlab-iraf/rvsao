# File rvsao/Util/velcon.x
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# April 30, 1997

# Subroutines to convert between velocity in km/sec and z = dwl/wl

# VEL2Z		Convert velocity in km/sec to z = dwl/wl
# VEL2Z1	Convert velocity in km/sec to 1+z
# Z2VEL		Convert velocity in z to km/sec

double procedure vel2z1 (velocity)

double	velocity	# Velocity in km/sec

double	c		# Velocity of light in km/sec
double	z1		# 1+z (unitless)

begin
	c = 299792.5
	z1 = sqrt ((1.d0 + (velocity / c)) / (1.d0 - velocity / c))
	return (z1)
end

double procedure vel2z (velocity)

double	velocity	# Velocity in km/sec

double	c		# Velocity of light in km/sec
double	z		# dwl / wl (unitless)

begin
	c = 299792.5
	z = vel2z1 (velocity) - 1.d0
	return (z)
end

double procedure z2vel (z)

double	z		# dwl / wl (unitless)

double	c		# Velocity of light in km/sec
double	velocity	# Velocity in km/sec

begin
	velocity = z12vel (1.d0 + z)
	return (velocity)
end


double procedure z12vel (z)

double	z		# dwl / wl (unitless)

double	c		# Velocity of light in km/sec
double	velocity	# Velocity in km/sec
double	z12		# (1+z)^2

begin
	z12 = z1 * z1
	velocity = c * ((z12 - 1.d0) / (z12 + 1.d0))
	return (velocity)
end
