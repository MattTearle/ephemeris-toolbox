function varargout = ephemeris(varargin)

% Get required inputs out of variable input list
[t,planets,outformat,outtype] = parseinputs(varargin,nargout);
% Get orbital element information and calculate heliocentric position of
% all planets, plus Earth (for reference)
oe = getorbitalelements(t,[planets;{'Earth'}]);
x = heliocentricposition(oe);
% Calculate cartesian position relative to Earth
dx = bsxfun(@minus,x(1:end-1,:,:),x(end,:,:));
% Convert geocentric cartesian position to desired coordinate frame
if ~strcmp(outformat,'equatorial')
    % ecliptic position
    position = cart2latlon(dx);
end
if ~strcmp(outformat,'ecliptic')
    % transform into equatorial
    x = dx;
    e = obliquity(reshape(t(:),1,1,[]));
    c = cos(e);
    s = sin(e);
    x(:,2,:) = bsxfun(@times,c,dx(:,2,:)) - bsxfun(@times,s,dx(:,3,:));
    x(:,3,:) = bsxfun(@times,s,dx(:,2,:)) + bsxfun(@times,c,dx(:,3,:));
    RAdec = cart2latlon(x);
    % RA typically given in time units
    RAdec(:,1,:) = convert(RAdec(:,1,:),'hours');
    if strcmp(outformat,'both')
        position = cat(2,position,RAdec);
    else
        position = RAdec;
    end
end

% Get output into the desired form
% Get dimensions
[np,no,nt] = size(position);
switch outtype
    case 'matrix'
        if nargout==1
            % Single array output
            varargout{1} = position;
        elseif (nargout==no)
            % Deal out components
            for k = 1:nargout
                varargout{k} = reshape(position(:,k,:),[np,nt]); %#ok<AGROW>
            end
        else
            error('Ephemeris:OutputNumberMismatch',...
                'Number of outputs doesn''t match output format')
        end
    case 'table'
        if np==1
            % One planet -> change order so rows = times
            position = permute(position,[3 2 1]);
        end
        % Now separate table from displayed output
        
        % Rearrange 3-D array into 2-D table
        position = reshape(permute(position,[1 3 2]),np*nt,no);
        % Expand planet names and times to match table rows
        time = repmat(t(:),np,1);
        planets = repelem(planets(:),nt,1);
        % Make table row names (for rearranged array) by combining planet
        % names with date strings
        rname = strcat(planets,{' '},cellstr(time));
        % Get variable names for table
        vname = {'Planet','Time'};
        switch outformat
            case 'ecliptic'
                vname = [vname,{'Longitude','Latitude'}];
            case 'equatorial'
                vname = [vname,{'RA','Declination'}];
            case 'both'
                vname = [vname,{'Longitude','Latitude','RA','Declination'}];
        end
        if nargout==1
            % Make output table
            varargout{1} = [table(planets,time),array2table(position)];
            varargout{1}.Properties.VariableNames = vname;
            varargout{1}.Properties.RowNames = rname;
        else
            error('Ephemeris:OutputNumberMismatch',...
                'Only one output can be specified when output format is ''table''')
        end
    case 'none'
        % Set labels for header lines and rows
        if np == 1;
            % One planet -> change order so rows = times
            position = permute(position,[3 2 1]);
            rname = cellstr(t(:));
            loopstr = planets(:);
            headerstr = 'Time';
        else
            loopstr = cellstr(t(:));
            rname = planets(:);
            headerstr = 'Planet';
        end
        % Add spacing to header based on length of elements in first col
        n = length(headerstr);
        longnm = max(n,size(char(rname),2));
        headerstr = [headerstr,repmat(' ',1,longnm-n+6)];
        % Loop over times (or single planet)
        for k = 1:length(loopstr)
            disp(loopstr{k})
            disp([headerstr,'Longitude      Latitude        RA      Declination'])
            % Get formatted display and split into cell array
            txt = regexp(fdisp(' %z    %1d    %t    %2d',position(:,:,k)),'\n','split');
            % Add row names and print
            txt = [rname,txt(:)]';
            fprintf(1,['%-',num2str(longnm),'s   %s\n'],txt{:});
        end
end


function [t,planets,outformat,outtype] = parseinputs(c,n)
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
        % Single string
        lowthis = lower(this);
        % See if input is valid output format
        if ismember(lowthis,{'ecliptic','equatorial','both'})
            % Check that there isn't already a format specified
            if isempty(outformat)
                outformat = lowthis;
            else
                error('Ephemeris:MultipleFormatSpecifications',...
                    'Multiple specifications of output format')
            end
        % See if input is valid output type
        elseif ismember(lowthis,{'table','matrix'})
            % Check that there isn't already a type specified
            if isempty(outtype)
                outtype = lowthis;
            else
                error('Ephemeris:MultipleTypeSpecifications',...
                    'Multiple specifications of output type')
            end
        % See if input is valid planet name
        else
            [~,~,idx] = intersect(lower(this),lower(allnames),'stable');
            if isempty(idx)
                thist = this;
            else
                planets = allnames(idx);
            end
        end
    elseif iscellstr(this)
        % Multiple strings => planet or datetimes
        [~,~,idx] = intersect(lower(this),lower(allnames),'stable');
        % If bad match to planet names, assume time (otherwise, planets)
        if (length(idx) < length(this))  || (length(this) > length(allnames))
            thist = this;
        else
            planets = allnames(idx);
        end
    elseif isempty(t)
        % Nothing else fits, so assume time (unless already specified)
        thist = this;
    else
        % Couldn't get a valid interpretation => bad input
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
        error('Ephemeris:BadEphemerisInput','Unrecognized input')
    end
end

% Set defaults for anything not yet defined
if isempty(t)
    t = datetime('now');
end
if isempty(planets)
    planets = allnames;
end
% Check that options make sense with the number of outputs
if isempty(outformat)
    if (n==0) || (n==4)
        outformat = 'both';
    else
        outformat = 'equatorial';
    end
elseif n==0
    warning('Output format ignored when there are no output arguments') %#ok<WNTAG>
    outformat = 'both';
end
if isempty(outtype)
    if n==0
        outtype = 'none';
    elseif n==1
        outtype = 'table';
    else
        outtype = 'matrix';
    end
elseif n==0
    warning('Output type ignored when there are no output arguments') %#ok<WNTAG>
    outtype = 'none';
end
