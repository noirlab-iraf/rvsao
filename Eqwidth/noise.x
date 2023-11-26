#		if (norm) {
#		    noise = (((flux1*norm1) + sflux1) / dsqrt (norm1 -1.d0)) +
#			(((flux2*norm2) + sflux2) / dsqrt (norm2 - 1.d0)) +
#			(((flux3*norm3) + sflux3) / dsqrt (norm3 - 1.d0))
#		    sigma = index * (noise / (flux1*norm1))
#		    }
#		else {
#		    noise = ((flux1 + sflux1) / dsqrt (norm1 - 1.d0)) +
#			((flux2 + sflux2) / dsqrt (norm2 - 1.d0)) +
#			((flux3 + sflux3) / dsqrt (norm3 - 1.d0))
#		    sigma = index * (noise / flux1)
#		    }
