function varargout = ephemeris(varargin)
[t,planets,outformat,outtype] = parseinputs(varargin);
oe = getorbitalelements(t,[planets;{'Earth'}]);
dx = oe(:,:,1:end-1) - oe(:,:,end)
cart2latlon(dx)

function [t,planets,outformat,outtype] = parseinputs(c)
% Start with empty arrays
t = [];
planets = {};
outformat = '';
outtype = '';
% Initialize flag
badinput = false;

% Create list of all possible objects to ask for
allnames = {'Sun';'Mercury';'Venus';'Mars';'Jupiter';...
    'Saturn';'Uranus';'Neptune';'Pluto'};

% Loop over all inputs
for k = 1:length(c)
    % Copy of current input
    this = c{k};
    % Initialize variable to test with input can be interpreted as a time
    thist = [];
    if ischar(this)
        lowthis = lower(this);
        if ismember(lowthis,{'ecliptic','equatorial','both'})
            if isempty(outformat)
                outformat = lowthis;
            else
                error('Multiple specifications of output format')
            end
        elseif ismember(lowthis,{'table','matrix'})
            if isempty(outtype)
                outtype = lowthis;
            else
                error('Multiple specifications of output type')
            end
        else
            [~,~,idx] = intersect(lower(this),lower(allnames),'stable');
            if isempty(idx)
                thist = this;
            else
                planets = allnames(idx);
            end
        end
    elseif iscellstr(this)
        [~,~,idx] = intersect(lower(this),lower(allnames),'stable');
        if (length(idx) < length(this))  || (length(this) > length(allnames))
            thist = this;
        else
            planets = allnames(idx);
        end
    elseif isempty(t)
        thist = this;
    else
        badinput = true;
    end
    if ~isempty(thist)
        if isempty(t)
            % Time hasn't been defined yet, so assume this is it
            % Pass the buck to DATETIME
            try
                t = datetime(thist);
            catch
                % Wasn't a time => don't know what it was
                badinput = true;
            end
        else
            % Two things that could be a time => bad
            badinput = true;
        end
    end
    if badinput
        error('Unrecognized input')
    end
end

% Set defaults for anything not yet defined
if isempty(t)
    t = datetime('now');
end
if isempty(planets)
    planets = allnames;
end
% TODO Check number of outputs
if isempty(outformat)
    outformat = 'both';
end
if isempty(outtype)
    outtype = 'table';
end
