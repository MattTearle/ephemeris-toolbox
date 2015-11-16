function planet = getorbitalelements(t,pnames)

allnames = {'Sun';'Mercury';'Venus';'Mars';'Jupiter';'Saturn';...
    'Uranus';'Neptune';'Pluto';'Earth'};
if nargin<2
    pnames = allnames;
    if nargin<1
        t = datetime('now');
    end
elseif ischar(pnames) && size(pnames,1)==1
    pnames = {pnames};
end
if ~iscellstr(pnames)
    error('Ephemeris:OrbitalElements:PlanetNameInput',...
        'Planet names must be given as a string or cell array of strings')
end
t = datetime(t);

[~,~,idx] = intersect(lower(pnames),lower(allnames),'stable');
numplanet = length(idx);
if (numplanet < length(pnames)) || (length(pnames) > length(allnames))
    error('Ephemeris:OrbitalElements:PlanetNames',...
        ['Invalid planet name.',char(10),'Allowable planets are: ',...
        sprintf('%s, ',allnames{1:end-1}),allnames{end},'.'])
end

OEAll = ...
[0 0 0 0 0 0;
 0.38709927      0.20563593      7.00497902      252.25032350     77.45779628     48.33076593
 0.72333566      0.00677672      3.39467605      181.97909950    131.60246718     76.67984255
 1.52371034      0.09339410      1.84969142       -4.55343205    -23.94362959     49.55953891
 5.20288700      0.04838624      1.30439695       34.39644051     14.72847983    100.47390909
 9.53667594      0.05386179      2.48599187       49.95424423     92.59887831    113.66242448
 19.18916464      0.04725744      0.77263783      313.23810451    170.95427630     74.01692503
 30.06992276      0.00859048      1.77004347      -55.12002969     44.96476227    131.78422574
 39.48211675      0.24882730     17.14001206      238.92903833    224.06891629    110.30393684
 1.00000261      0.01671123     -0.00001531      100.46457166    102.93768193      0.0];

DeltaAll = ...
[0 0 0 0 0 0;
 0.00000037      0.00001906     -0.00594749   149472.67411175      0.16047689     -0.12534081
 0.00000390     -0.00004107     -0.00078890    58517.81538729      0.00268329     -0.27769418
 0.00001847      0.00007882     -0.00813131    19140.30268499      0.44441088     -0.29257343
-0.00011607     -0.00013253     -0.00183714     3034.74612775      0.21252668      0.20469106
-0.00125060     -0.00050991      0.00193609     1222.49362201     -0.41897216     -0.28867794
-0.00196176     -0.00004397     -0.00242939      428.48202785      0.40805281      0.04240589
 0.00026291      0.00005105      0.00035372      218.45945325     -0.32241464     -0.00508664
-0.00031596      0.00005170      0.00004818      145.20780515     -0.04062942     -0.01183482
 0.00000562     -0.00004392     -0.01294668    35999.37244981      0.32327364      0.0];

OE = OEAll(idx,:);
Delta = DeltaAll(idx,:);
dt = j2000(t(:))/36525;

planet = zeros(length(dt),6,numplanet);
for k = 1:numplanet
    planet(:,:,k) = bsxfun(@plus,OE(k,:),bsxfun(@times,dt,Delta(k,:)));
end