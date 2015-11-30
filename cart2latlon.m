function y = cart2latlon(x)
% Convert Cartesian coordinates (x,y,z) to latitude & longitude

lat = normalizeinrange(angl.atan(x(:,2,:),x(:,1,:)),1);
lon = normalizeinrange(angl.asin(x(:,3,:)./sqrt(sum(x.^2,2))),0);
y = [lat,lon];
