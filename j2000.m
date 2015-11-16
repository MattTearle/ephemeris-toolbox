function j2k = j2000(t)
if ~isdatetime(t)
    error('Ephemeris:JulianDate','Date must be a datetime variable')
end
j2k = juliandate(t) - 2451545;
