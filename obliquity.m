function ob = obliquity(t)
% DE200: e = 23� 26' 21.45" - 46.815" T - 0.0006" T^2 + 0.00181" T^3
% T = Julian centuries from J2000
t = j2000(t)/36525;
c = [5.027777777777778e-07;
    -1.666666666666667e-07;
    -0.013004166666667;
    23.43929166666667;];
ob = angl(polyval(c,t));
