function x = heliocentricposition(orbelements)

% a = semimajor;
% e = eccentricity;
% I = inclination;
% L = meanlong;
% w = perihelion;
% Q = ascendingnode;

[numtime,~,numplanet] = size(orbelements);

a = reshape(orbelements(:,1,:),numtime,numplanet);
e = reshape(orbelements(:,2,:),numtime,numplanet);
I = reshape(orbelements(:,3,:),numtime,numplanet);
L = reshape(orbelements(:,4,:),numtime,numplanet);
w = reshape(orbelements(:,5,:),numtime,numplanet);
Q = reshape(orbelements(:,6,:),numtime,numplanet);

% a = orbelements(:,1,:);
% e = orbelements(:,2,:);
% I = orbelements(:,3,:);
% L = orbelements(:,4,:);
% w = orbelements(:,5,:);
% Q = orbelements(:,6,:);

% if nargin==1
%         e = a.eccentricity;
%         I = a.inclination;
%         L = a.meanlong;
%         w = a.perihelion;
%         Q = a.ascendingnode;
%         a = a.semimajor;
% elseif nargin~=6
%     error('Wrong number of inputs')
% end

anom = pi*mod(L - w,360)/180;
Kep = @(E) E - e.*sin(E) - anom;
dKep = @(E) 1 - e.*cos(E);
E = anom;
dE = Kep(E);
while (max(abs(dE(:))) > 1e-6)
    E = E - dE./dKep(E);
    dE = Kep(E);
end
% E = fzero(Kep,anom);

x = zeros(numtime,3,numplanet);
for k = 1:numplanet
    for j = 1:numtime
        
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
