1161,1162c1161,1162
< int	i, i1, i2
< double	a, b, r1, r2, w1, w2, di1, di2, wt
---
> int	i, i1, i2, fluxi, flux2
> double	a, b, r1, r2, w1, w2, di1, di2, wt, sigma
1227c1227,1229
< 	    flux = wt * spec[i1]
---
> 	    fluxi = wt * spec[i1]
> 	    flux = fluxi
> 	    flux2 = fluxi * fluxi
1233c1235,1237
< 		flux = flux + wt * spec[i]
---
> 		fluxi = wt * spec[i]
> 		flux = flux + fluxi
> 		flux2 = flux2 + (fluxi * fluxi)
1239c1243,1245
< 	    flux = flux + wt * spec[i2]
---
> 	    fluxi = wt * spec[i2]
> 	    flux = flux + fluxi
> 	    flux2 = flux2 + (fluxi * fluxi)
1240a1247,1249
> 
> 	    mean = flux / normb
> 	    sigma = (flux2 - (flux * flux)) / normb
