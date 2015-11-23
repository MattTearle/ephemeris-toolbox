function y = cart2latlon(x)

lat = angl.atan(x(:,2,:),x(:,1,:));
lon = angl.asin(x(:,3,:)./sqrt(sum(x.^2,2)));
y = [lat,lon];
