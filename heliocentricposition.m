function x = heliocentricposition(orbelements)

[numplanet,~,numtime] = size(orbelements);

% a = semimajor;
% e = eccentricity;
% I = inclination;
% L = meanlong;
% w = perihelion;
% Q = ascendingnode;
a = reshape(orbelements(:,1,:),numplanet,numtime);
e = reshape(orbelements(:,2,:),numplanet,numtime);
I = reshape(orbelements(:,3,:),numplanet,numtime);
L = reshape(orbelements(:,4,:),numplanet,numtime);
w = reshape(orbelements(:,5,:),numplanet,numtime);
Q = reshape(orbelements(:,6,:),numplanet,numtime);

anom = pi*mod(L - w,360)/180;
Kep = @(E) E - e.*sin(E) - anom;
dKep = @(E) 1 - e.*cos(E);
E = anom;
dE = Kep(E);
while (max(abs(dE(:))) > 1e-6)
    E = E - dE./dKep(E);
    dE = Kep(E);
end

x = zeros(numplanet,3,numtime);
for k = 1:numtime
    for j = 1:numplanet
        
        xprime = [a(j,k)*(cos(E(j,k))-e(j,k));...
            a(j,k)*sqrt(1-e(j,k)^2)*sin(E(j,k))];
        
        sQ = sind(Q(j,k));
        cQ = cosd(Q(j,k));
        cw = cosd(w(j,k)-Q(j,k));
        sw = sind(w(j,k)-Q(j,k));
        cI = cosd(I(j,k));
        sI = sind(I(j,k));
        x(j,:,k) = [cw*cQ-sw*sQ*cI,-sw*cQ-cw*sQ*cI;...
            cw*sQ+sw*cQ*cI,-sw*sQ+cw*cQ*cI;sw*sI,cw*sI]*xprime;
    end
end
