classdef angl
    % Angle variable: stores angles measured in degrees, radians, or hours
    %
    %  A = ANGL(X,U) creates an angle array A with values specified by the
    %  array X and units specified by the string or cell array of strings
    %  U. If U is a single string, all values in X will have the same
    %  units. If U is a cell array of strings, it must have the same
    %  dimensions as X. Units can be specified as 'degrees', 'radians', or
    %  'hours'.
    %
    %  A = ANGL(X) creates an angle array, using default units of
    %  'degrees'.
    %
    %  Examples:
    %  a = angl(45,'degrees')
    %  a = angl([45,90,200],'degrees')
    %  a = angl([45,6,pi/2],{'degrees','hours','radians'})
    %  a = angl(pi*rand(3,4),'radians')
    
    properties (Access = private)
        value
        units
    end
    
    methods
        % Constructor
        function c = angl(x,u)
            if nargin>0
                % Set default unit
                if nargin<2
                    u = 1;
                end
                % Validate input of values (x)
                if isempty(x)
                    error('Angle value cannot be empty')
                elseif ~isnumeric(x)
                    error('Angle values must be numeric')
                else
                    xsize = size(x);
                    xlen  = numel(x);
                end
                % Validate input of units (u)
                if ischar(u)
                    % For simpler programming further on, make u a cell array
                    u = {u};
                end
                if iscellstr(u)
                    u = lower(u);
                    unum = zeros(size(u));
                    unum(strcmp(u,'degrees')) = 1;
                    unum(strcmp(u,'radians')) = 2;
                    unum(strcmp(u,'hours')) = 3;
                elseif isnumeric(u)
                    unum = u;
                else
                    error('Ephemeris:Angle:BadUnits',...
                        'Units must be "degrees", "radians", or "hours"')
                end
                % Check units
                if ~all(ismember(unum(:),[1,2,3]))
                    error('Ephemeris:Angle:BadUnits',...
                        'Units must be "degrees", "radians", or "hours"')
                end
                % How are the units specified (individually or a single unit)?
                issameunits = (numel(u)==1);
                if issameunits
                    % Extract the one unit
                    thisu = unum(1);
                    % Multiple units => check dimensions against x
                elseif ~isequal(xsize,size(u))
                    error('Ephemeris:Angle:MismatchedDimensions',...
                        'Dimension mismatch of values and units')
                end
                % Build the angle!
                % Scalar
                if (xlen == 1)
                    c.value = x;
                    c.units = thisu;
                    % Array
                else
                    for k = xlen:-1:1
                        % Get units (if defined individually)
                        if ~issameunits
                            thisu = unum(k);
                        end
                        % Call scalar constructor on this single element
                        c(k).value = x(k);
                        c(k).units = thisu;
                    end % loop over array of values (x)
                    c = reshape(c,xsize);
                end % if (scalar vs array)
                % Ensure that everything is put into the base interval
                % (-360,360), (-2pi,2pi), (-24,24)
                c = normalize(c);
            end
        end
        
        % TODO: Comments!
        function varargout = fdisp(fstr,A)
            % FDISP(FORMAT,A) displays the array of angles A to the screen,
            % according to the FORMAT string.
            % X = FDISP(FORMAT,A) returns the formatted display to a string X.
            %
            % The FORMAT specification string is SPRINTF-style, but uses
            % angle-specific formats: %d, %t, and %z
            % %d = degrees/minutes/seconds: 42.314 -> +42°18'50"
            % %t = time (hours/minutes/seconds): 42.314 -> 02:49:15
            % %z = zodiac (degrees-sign/minutes/seconds): 42.314 -> 12 Tau 18'50"
            % Other SPRINTF markups such as \t and are valid
            %
            % A scalar format specifier is applied to all elements.
            % Multiple specifiers are applied to the columns of A; the
            % number of specifiers must match the number of columns of A.
            %
            % Example:
            % fdisp('%z\t%2d -- %t',angl([42,3.14,666]))
            %   12 Tau 00'00"	+03°08'24" -- 20:24:00
            
            % Find the angle-specific formats
            [idx1,idx2] = regexp(fstr,'%\d*[dtz]');
            % Number of formats & columns of A
            n = length(idx1);
            ncol = size(A,2);
            if n==1
                % Repeat single specifier ncol times
                fstr = deblank(repmat([fstr,'  '],1,ncol));
                [idx1,idx2] = regexp(fstr,'%\d*[dtz]');
            elseif ncol ~= n
                error('Ephemeris:FDisp:MismatchedDimensions',...
                    'Need to specify one format for each column')
            end
            % Make space for output text (flip size for transpose)
            outtext = cell(fliplr(size(A)));
            % Loop over columns of A
            for k = 1:ncol
                % Which format? Look at the last letter
                % Use the *form methods to generate the strings for the kth
                % column of A; store as the kth row of the output text
                switch fstr(idx2(k))
                    case 't'
                        outtext(k,:) = timeform(convert(A(:,k),'hours'))';
                    case 'z'
                        outtext(k,:) = zodiacform(convert(A(:,k),'degrees'))';
                    case 'd'
                        % Is there a number in front of the d?
                        if idx2(k) == idx1(k)+1 % no
                            outtext(k,:) = sexiform(convert(A(:,k),'degrees'))';
                        else % yes
                            n = str2double(fstr((idx1(k)+1):(idx2(k)-1)));
                            outtext(k,:) = sexiform(convert(A(:,k),'degrees'),n)';
                        end
                    otherwise
                        error('Unknown format specification')
                end
            end % loop over columns
            % Use SPRINTF to make the output string from the cell array
            % Replace all angle-specific formats with %s (meaning that
            % everything else in the format string will be preserved)
            outstr = sprintf([regexprep(fstr,'%\d*[dtz]','%s'),'\n'],outtext{:});
            % Strip off the last "\n"
            outstr(end) = [];
            
            % Return the string or display it to the screen
            if nargout
                varargout{1} = outstr;
            else
                disp(outstr)
            end
        end
        
        function disp(A)
            % DISP(A) displays the angle array A, without printing the
            % array name.
            
            % If 2-D array, just display
            if ismatrix(A)
                matdisp(A)
            else
                % Otherwise, have to loop over other dimensions
                d = size(A);
                % Make an array of all indices (in display order), ignoring
                % first two dimensions b/c each display is a matrix
                dimvec = angl.allindices(d(3:end));
                for k = 1:size(dimvec,1)
                    % Make a string of the indices
                    idx = sprintf('%d,',dimvec(k,:));
                    idx = ['(:,:,',idx(1:end-1),')'];
                    % Display "x(:,:,d1,d2,...) = "
                    disp([inputname(1),idx,' ='])
                    % Display contents by evaluating index statement
                    disp(eval(['A',idx]))
                end
            end
        end
        
        function display(A) %#ok<DISPLAY>
            % Display style = "name = ..." unless n-D array
            if ismatrix(A)
                disp([inputname(1),' ='])
            end
            % Evil eval, but it works...
            % Keeps name of input when passing on to DISP
            eval([inputname(1),' = A;']);
            eval(['disp(',inputname(1),')']);
        end
        
        function s = zodiacform(a)
            % S = ZODIACFORM(A) returns a cell array of strings
            % representing the angle variable A in Zodiac form
            % (deg SIGN min'sec")
            
            % Put into [0 360] range
            a = normalizeinrange(a,1);
            % Get degree/min/sec components
            [~,deg,min,sec] = dec2sex(a);
            % Determine zodiac sign (30-degree blocks)
            s = floor(deg/30);
            % Remainder (degrees within given sign)
            deg = deg - 30*s;
            zsign = {' Ari ';' Tau ';' Gem ';' Cnc ';' Leo ';' Vir ';' Lib ';' Sco ';' Sgr ';' Cap ';' Aqr ';' Psc '};
            % Build formatted string
            str = strcat(cellstr(num2str(deg,'%02d')),zsign(s+1),num2str(min,'%02d'),char(39),num2str(sec,'%02d'),'"');
            % Ensure output is same shape as input (dec2sex changes shape)
            s = reshape(str,size(a));
        end
        
        function s = timeform(a)
            % S = TIMEFORM(A) returns a cell array of strings representing
            % the angle variable A in "time form" (hh:mm:ss)
            
            % Put into [0 24] range
            a = normalizeinrange(a,1);
            % Get hour/min/sec components
            dms = dec2sex(a);
            % Build formatted string
            str = cellstr(strcat(num2str(dms(:,2),'%02d'),':',num2str(dms(:,3),'%02d'),':',num2str(dms(:,4),'%02d')));
            % Ensure output is same shape as input (dec2sex changes shape)
            s = reshape(str,size(a));
        end
        
        function s = sexiform(a,n)
            % S = SEXIFORM(A) returns a cell array of strings representing
            % the angle variable A in "sexigesimal form" (+/- deg°mm'ss")
            
            % Get degree/min/sec components
            [s,deg,min,sec] = dec2sex(a);
            % Build the +/- sign part of the output string
            sgn = repmat('+',size(s));
            sgn(s<0) = '-';
            % Create format string "%0nd", where n is given or determined
            if nargin<2
                % n not specified -- calculate it
                % find number of digits needed for each value of deg
                % (with a minimum of 1 digit); take the max
                fmt = ['%0',num2str(max(max(floor(log10(abs(deg)))+1),1),'%1d'),'d'];
            else
                fmt = ['%0',num2str(n,'%1d'),'d'];
            end
            % Build formatted string (char 176 = degree symbol)
            str = cellstr(strcat(sgn,num2str(deg,fmt),char(176),num2str(min,'%02d'),char(39),num2str(sec,'%02d'),'"'));
            % Ensure output is same shape as input (dec2sex changes shape)
            s = reshape(str,size(a));
        end
        
        % Mathematical operators
        function z = plus(x,y)
            
            % numeric + angle -> convert numeric to angle
            if isnumeric(x)
                % Assume number has same units as angle
                x = angl(x,y.units);
            end
            nx = numel(x);
            ny = numel(y);
            if (nx==1)
                % x is scalar
                z = y;
                if (ny==1)
                    % y is also scalar => scalar + scalar
                    if isnumeric(y)
                        % angle + numeric -> convert numeric to angle
                        % Assume number has same units as angle
                        y = angl(y,x.units);
                    elseif (x.units ~= y.units)
                        % angle + angle, but with different units
                        % convert to same units
                        y = convert(y,x.units);
                    end
                    % Do the addition
                    zval = x.value + y.value;
                    z = angl(zval,x.units);
                else
                    % y is array => scalar + array
                    % Use loop to do scalar expansion
                    for k = 1:ny
                        z(k) = x + y(k);
                    end
                end
            elseif (ny==1)
                % x is array, y is scalar => array + scalar
                z = x;
                % Use loop to do scalar expansion
                for k = 1:nx
                    z(k) = x(k) + y;
                end
            elseif isequal(size(x),size(y))
                % array + array => check for same dimensions
                z = x;
                % Use loop to perform elementwise addition
                for k = 1:nx
                    z(k) = x(k) + y(k);
                end
            else
                % array + array but not same dimensions
                error('Ephemeris:AngleMath:MismatchedDimensions',...
                    'Dimensions do not agree')
            end
        end
        
        % TODO: check minus math (range of negative output)
        function z = minus(x,y)
            % A - B means the angle between two lines from the origin with
            % polar angles A and B. The smaller of the two complementary
            % angles is taken. The result is therefore in the range
            % [-180 180] degrees ([-pi pi] rad, [-12 12] hrs).
            
            % numeric + angle -> convert numeric to angle
            if isnumeric(x)
                % Assume number has same units as angle
                x = angl(x,y.units);
            end
            nx = numel(x);
            ny = numel(y);
            % x is scalar
            if (nx==1)
                z = y;
                if (ny==1)
                    % y is also scalar => scalar - scalar
                    if isnumeric(y)
                        % angle - numeric -> convert numeric to angle
                        % Assume number has same units as angle
                        y = angl(y,x.units);
                    elseif (x.units ~= y.units)
                        % angle - angle, but with different units
                        % convert to same units
                        y = convert(y,x.units);
                    end
                    % Do the subtraction
                    zval = x.value - y.value;
                    z = angl(zval,x.units);
                    % Check range of result to determine which way the
                    % subtraction should go
                    zc = convert(z,'degrees');
                    if (zc.value>180)
                        zc.value = 360 - zc.value;
                        z = convert(zc,z.units);
                    end
                else
                    % y is array => scalar - array
                    % Use loop to do scalar expansion
                    for k = 1:ny
                        z(k) = x - y(k);
                    end
                end
            elseif (ny==1)
                % x is array, y is scalar => array - scalar
                z = x;
                % Use loop to do scalar expansion
                for k = 1:nx
                    z(k) = x(k) - y;
                end
            elseif isequal(size(x),size(y))
                % array - array => check for same dimensions
                z = x;
                % Use loop to perform elementwise subtraction
                for k = 1:nx
                    z(k) = x(k) - y(k);
                end
            else
                % array - array but not same dimensions
                error('Ephemeris:AngleMath:MismatchedDimensions',...
                    'Dimensions do not agree')
            end
        end
        
        function z = mtimes(x,y)
            % no matrix multiplication for angles
            % so mtimes (*) == times (.*)
            z = x.*y;
        end
        
        function z = times(x,y)
            
            % angle * numeric = numeric * angle
            if isnumeric(y)
                z = y.*x;
            else
                % now x must be numeric
                if ~isnumeric(x)
                    error('Ephemeris:AngleMath:NonScalarMultiply',...
                        'Only scalar multiplication is defined for angles')
                end
                % scalar multiply
                if isscalar(x)
                    z = y;
                    for k = 1:numel(z)
                        z(k).value = x*y(k).value;
                    end
                elseif isequal(size(x),size(y))
                    % Array multiply (numeric .* angle)
                    z = y;
                    for k = 1:numel(y)
                        z(k).value = x(k).value*y(k).value;
                    end
                else
                    error('Ephemeris:AngleMath:MismatchedDimensions',...
                        'Dimension mismatch for scalar multiplication')
                end
            end
            z = normalize(z);
        end
        
        function z = mrdivide(x,y)
            
            % Allow "x/y" if y is scalar (= x./y)
            if isscalar(y)
                z = x./y;
            else
                error('Ephemeris:AngleMath:NonScalarDivide',...
                    'Only scalar division is defined for angles')
            end
        end
        
        function z = rdivide(x,y)
            
            % only angle/numeric makes sense
            if isnumeric(y)
                if isscalar(y)
                    % angle/scalar
                    z = x;
                    for k = 1:numel(z)
                        z(k).value = z(k).value/y;
                    end
                elseif isequal(size(x),size(y))
                    % array divide (angle./numeric)
                    z = x;
                    for k = 1:numel(z)
                        z(k).value = z(k).value/y(k);
                    end
                else
                    error('Ephemeris:AngleMath:MismatchedDimensions',...
                        'Dimension mismatch for scalar division')
                end
            else
                error('Ephemeris:AngleMath:NonScalarDivide',...
                    'Only scalar division is defined for angles')
            end
            z = normalize(z);
        end
        
        function x = uminus(x)
            for k = 1:numel(x)
                x(k).value = -(x(k).value);
            end
        end
        
        function y = convert(x,newunits)
            % Y = CONVERT(X,NEWUNITS) converts the angle array X into an
            % angle array Y that has units of NEWUNITS.
            
            y = x;
            % check newunits
            if ischar(newunits)
                switch newunits
                    case 'degrees'
                        newunits = 1;
                    case 'radians'
                        newunits = 2;
                    case 'hours'
                        newunits = 3;
                    otherwise
                        newunits = 0;
                end
            end
            % String has been converted to numeric (1, 2, or 3, with 0
            % being an error). Anything else is disallowed.
            if ~isnumeric(newunits) || ~ismember(newunits,[1,2,3])
                error('Ephemeris:Angle:BadUnits',...
                    'Units must be "degrees", "radians", or "hours"')
            end
            % Convert each element of x
            for k = 1:numel(x)
                if (x(k).units == newunits)
                    % No conversion necessary
                    y(k) = x(k);
                else
                    % Convert...
                    switch x(k).units
                        case 1
                            % ...from degrees...
                            if (newunits == 2)
                                % ...into radians
                                yval = x(k).value*pi/180;
                            elseif (newunits == 3)
                                % ...into hours
                                yval = x(k).value/15;
                            end
                        case 2
                            % ...from radians...
                            if (newunits == 1)
                                % ...into degrees
                                yval = x(k).value*180/pi;
                            elseif (newunits == 3)
                                % ...into hours
                                yval = x(k).value*12/pi;
                            end
                        case 3
                            % ...from hours...
                            if (newunits == 1)
                                % ...into degrees
                                yval = x(k).value*15;
                            elseif (newunits == 2)
                                % ...into radians
                                yval = x(k).value*pi/12;
                            end
                    end
                    y(k).value = yval;
                    y(k).units = newunits; %= angl(yval,newunits);
                end
            end
        end
        
        function x = normalize(x)
            % Y = NORMALIZE(X) returns an angle Y that has the same
            % physical meaning as the angle X but with the value of the
            % angle restricted to [0 360) degrees, [0 2*pi) radians, or
            % [0 24) hours.
            
            % Use more general normalization method
            x = normalizeinrange(x,2);
        end
        
        function x = normalizeinrange(x,rng)
            % Y = NORMALIZEINRANGE(X) returns an angle Y that has the same
            % physical meaning as the angle X but with the value of the
            % angle restricted to lie in the range RNG. The range can be
            % specified as a two-element vector (e.g. [0 pi] or [-90 90])
            % or as the scalar values -1, 0, 1, 2, which correspond to:
            % RNG |  degrees     radians       hours
            % ------------------------------------------
            % -1  |  [-360,0)    [-2*pi 0)     [-24 0)
            %  0  |  [-180,180)  [-pi pi)      [-12 12)
            %  1  |  [0 360)     [0 2*pi)      [0 24)
            %  2  |  [-360 360)  [-2*pi 2*pi)  [-24 24)
            
            % Check input
            if ~isnumeric(rng) || ~isvector(rng) || (isscalar(rng) && ~ismember(rng,-1:2))
                error('Ephemeris:Angle:BadRange',...
                    'Range must be specified as a numeric vector [a,b] or scalar (-1, 0, 1, or 2)')
            end
            % Normalize each element
            for k = 1:numel(x)
                if isscalar(rng)
                    % Interpret scalar RNG (-1, 0, 1, or 2)
                    % Cunning trick: units are stored internally as 1, 2 or
                    % 3. Take 4*units (4, 8, or 12) and add RNG. Result is
                    % a unique value from 3 to 14:
                    switch 4*x(k).units + rng
                        % degrees
                        case 3
                            r = [-360,0];
                        case 4
                            r = [-180,180];
                        case 5
                            r = [0,360];
                        case 6
                            r = [-360,360];
                        % radians
                        case 7
                            r = [-2*pi,0];
                        case 8
                            r = [-pi,pi];
                        case 9
                            r = [0,2*pi];
                        case 10
                            r = [-2*pi,2*pi];
                        % hours
                        case 11
                            r = [-24,0];
                        case 12
                            r = [-12,12];
                        case 13
                            r = [0,24];
                        case 14
                            r = [-24,24];
                    end
                else
                    % Otherwise, just use given numeric interval
                    r = rng;
                end
                % Use MOD to wrap value into an interval of length
                % r(2)-r(1) (= diff(r)), with minimum value r(1)
                x(k).value = r(1) + mod(x(k).value+r(1),diff(r));
            end
        end
        
        function x = double(a,u)
            % X = DOUBLE(A) converts the angle array A to the double array
            % X containing the values of the angles in degrees.
            % X = DOUBLE(A,UNITS) uses the values of the angles in the
            % given units ("degrees", "radians", or "hours").
            % X = DOUBLE(A,'none') does not perform any conversion. The
            % values of the angles are left in their original units.
            
            % Default = degrees
            if nargin<2
                u = 1;
            end
            % Check for 'none' flag (convert unless U is "no*")
            if ~ischar(u) || ~strncmpi(u,'no',2)
                a = convert(a,u);
            end
            % Extract values and reshape to match input
            x = reshape([a.value],size(a));
        end
        
        function varargout = dec2sex(x,frac)
            % SDMS = DEC2SEX(X) converts the vector of decimal angles X to
            % sexigesimal form. The result is a matrix SDMS with four
            % columns: sign (+/- 1), degrees, minutes, seconds.
            if nargin<2
                frac = false;
            end
            if ~isvector(x)
                error('x cannot be a matrix.')
            end
            y = [x.value];
            y = y(:);
            sgn = sign(y);
            y = abs(y);
            d = floor(y);
            m = floor((y-d)*60);
            s = (y-d-m/60)*3600;
            if ~frac
                s = round(s);
            end
            idx = (abs(s-60) < 1e-10);
            s(idx) = 0;
            m(idx) = m(idx)+1;
            idx = (m == 60);
            m(idx) = 0;
            d(idx) = d(idx)+1;
%             d = sgn.*d;
            if nargout<2
                varargout{1} = [sgn,d,m,s];
            elseif nargout==4
                varargout = {sgn,d,m,s};
            else
                error('Wrong number of outputs')
            end
        end
        
        % Basic rounding math operations (all elementwise)
        function x = floor(x)
            for k = 1:numel(x)
                x(k).value = floor(x(k).value);
            end
        end
        function x = round(x)
            for k = 1:numel(x)
                x(k).value = round(x(k).value);
            end
        end
        function x = ceil(x)
            for k = 1:numel(x)
                x(k).value = ceil(x(k).value);
            end
        end
        
        % Trig functions
        % Convert to radians and apply elementwise
        function y = sin(x)
            z = convert(x,2);
            y = reshape(sin([z.value]),size(z));
        end
        function y = cos(x)
            z = convert(x,2);
            y = reshape(cos([z.value]),size(z));
        end
        function y = tan(x)
            z = convert(x,2);
            y = reshape(tan([z.value]),size(z));
        end

%         function plot(lon,lat,mag)
%             if size(lon,2)==2
%                 x = lon(:,1);
%                 y = lon(:,2);
%                 mag = lat;
%             else
%                 x = lon;
%                 y = lat;
%             end
%             x = normalizeinrange(convert(x,1),0);
%             y = normalizeinrange(convert(y,1),0);
%             scatter([x.value],[y.value],10.^(-mag),'w','filled')
%             axis([-200,200,-90,90])
%             set(gca,'color','k')
%         end
        
        function scatter(lon,lat,varargin)
            if size(lon,2)==2
                x = lon(:,1);
                y = lon(:,2);
                if nargin > 1
                    varargin = [{lat},varargin];
                end
            else
                x = lon;
                y = lat;
            end
            x = normalizeinrange(convert(x,1),0);
            y = normalizeinrange(convert(y,1),0);
            scatter([x.value],[y.value],varargin{:})
            axis([-200,200,-90,90])
        end
        
    end
    
    methods(Access=private)
        function matdisp(A)
            
            % Matrix (ndims = 2) display
            sizeA = size(A);    % Get dimensions
            vals = reshape([A.value],sizeA);    % Extract values
            % Logical flag to show whether any of the values in each column
            % have some fractional parts (T) or they are all integers (F)
            frac = ~all(round(vals)==vals,1);
            % Number of characters needed for the integer part
            iwdth = max(max(floor(log10(max(abs(vals),1)))+1+(vals<0),1));
            % Format string should be either "%nd" or "%n.4f"
            % (The first when FRAC is T, the second when it's F)
            % n is IWDTH is FRAC is T and (IWDTH+5) when FRAC is F (5 = 1
            % decimal point + 4 decimal values)
            fstr = cellstr(num2str(4*frac'));
            tstr = cellstr(num2str((iwdth+5*frac)'));
            % "d" = char(100), "f" = char(102)
            fmt = strcat('%',tstr,'.',fstr,cellstr(char(100+2*frac')));
            % This approach makes "%n.0d" instead of "%nd"; fix that
            fmt = regexprep(fmt,'.0d','d');
            
            % Now print
            % Loop over rows (each line)
            for j = 1:sizeA(1)
                fprintf(1,'    ');
                % Loop over columns
                for k = 1:sizeA(2)
                    % Get the current value to print
                    x = vals(j,k);
                    % Print the value + unit abbreviation
                    switch A(j,k).units
                        case 1
                            fprintf(1,[fmt{k},' deg'],x);
                        case 2
                            fprintf(1,[fmt{k},' rad'],x);
                        case 3
                            fprintf(1,[fmt{k},' hrs'],x);
                    end
                    % Add space, except at the end of the line
                    if k<sizeA(2)
                        fprintf(1,'    ');
                    end
                end % columns
                % New line
                fprintf(1,'\n');
            end % rows
        end
        

    end
    
    methods(Static)
        % Inverse trig functions
        function y = asin(x)
            y = angl(asind(x));
        end
        function y = acos(x)
            y = angl(acosd(x));
        end
        function z = atan(y,x)
            if nargin<2
                x = ones(size(y));
            end
            idx = (sign(x)<=0);
            z = angl(atand(y./x));
            z(idx) = z(idx)+180;
        end
        function z = acot(x,y)
            if nargin<2
                y =ones(size(x));
            end
            idx = (sign(y)<=0);
            z = angl(acotd(x./y));
            z(idx) = z(idx)+180;
        end
        
        function a = sex2dec(d,m,s,units)
            % Make angle 
            switch nargin
                case 1
                    units = 1;
                    x = d(:,1);
                    y = d(:,2);
                    z = d(:,3);
                case 2
                    if ischar(m)
                        switch m
                            case 'degrees'
                                units = 1;
                            case 'radians'
                                units = 2;
                            case 'hours'
                                units = 3;
                        end
                    else
                        units = m;
                    end
                    x = d(:,1);
                    y = d(:,2);
                    z = d(:,3);
                case 3
                    units = 1;
                    x = d;
                    y = m;
                    z = s;
                case 4
                    x = d;
                    y = m;
                    z = s;
            end
            a = angl((2*(x>=0)-1).*(abs(x) + y/60 + z/3600),units);
        end
    end
    
    methods(Static,Access=private)
        function idx = allindices(d)
            
            if isscalar(d)
                idx = (1:d)';
            else
                idx2 = angl.allindices(d(2:end));
                s = size(idx2,1);
                idx2 = repelem(idx2,d(1),1);
                idx = [repmat((1:d(1))',s,1),idx2];
            end
            
        end
    end
    
end

