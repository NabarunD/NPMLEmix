# from new carina dataset let's extract ascension and declination
cnew = read.csv('carina_published.csv', header = TRUE)

# ascension is represented by RA_hour, RA_minute, RA_second
# needs to be converted to degrees
ra = (cnew[,'RA_hour'] + cnew[,'RA_minute']/60 + cnew[,'RA_second']/3600)/24 * 2 * pi
# how do we get that?
# declination represented in Declin_deg, Declin_arcmin, Declin_arcsec
# need to convert this to degrees
de = (cnew[,'Declin_deg'] + cnew[,'Declin_arcmin']/60 + cnew[,'Declin_arcsec']/3600)/360 * 2 * pi

# we need to know RA of the center of galaxy
# otherwise we can not proceed
# ra0 = ?
# de0 = ?

# given those we need to transform as follows
xi = cos(de) * sin(ra - ra0)/(sin(de0) * cos(de) + cos(de0) * cos(de) * cos(ra - ra0))
eta = (cos(de0) * sind(de) - sin(de0) * cos(de) * cos(ra - ra0))/
  (sin(de0) * sin(de) + cos(de0) * cos(de) * cos(ra - ra0))

R = sqrt(xi^2 + eta^2)

# from matt
# RA = 100.40292
# Dec = 50.96611
# do this....