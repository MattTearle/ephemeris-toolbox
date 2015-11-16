function x = eclipticpos(planet,t)

p = getorbitalelements(t,{planet,'Earth'});
% a1 = p(1).semimajor;
% e1 = p(1).eccentricity;
% I1 = p(1).inclination;
% L1 = p(1).meanlong;
% w1 = p(1).perihelion;
% Q1 = p(1).ascendingnode;
% a2 = p(2).semimajor;
% e2 = p(2).eccentricity;
% I2 = p(2).inclination;
% L2 = p(2).meanlong;
% w2 = p(2).perihelion;
% Q2 = p(2).ascendingnode;
dx = heliocentricposition(p(1)) - ...
    heliocentricposition(p(2));
x = cart2latlon(dx);
