function a = arrow(ax,startpt,endpt,varargin)
% ARROW(START,END) adds an arrow to the current axes from the (x,y) point
% defined by the 2-element vector START to the (x,y) point defined by the
% 2-element vector END.
%
% ARROW(AX,START,END) adds an arrow to the axes AX.
%
% ARROW(...,STR) uses the line style specification string STR to set the
% color and line style of the arrow. The string can be any line style
% specification string that is accepted by PLOT (e.g. 'k:', 'm-*', etc.).
% However markers are ignored.
%
% ARROW(...,NAME,VALUE) formats the arrow with additional options specified
% by one or more Name,Value pair arguments. Available options are:
% 'Color'        Color of the arrow (line and arrowhead)
%                String or 1-by-3 RGB vector
%                Default: 'k' ([0 0 0])
% 'HeadLength'   Length of the arrowhead in pts (1/72 inch)
%                Positive numeric scalar
%                Default: 10
% 'HeadWidth'    Width of the arrowhead in pts (1/72 inch)
%                Positive numeric scalar
%                Default: 10
%     Note that HeadLength is scaled if necessary to be no more than
%     half the length of the arrow. HeadWidth is scaled equivalently,
%     to maintain the aspect ratio.
% 'HeadIndent'   Amount of indentation of the back of the arrowhead
%                Numeric scalar between -1 and 1
%                Default: 0.2
% 'LineStyle'    Style of arrow line
%                String (any line style string accepted by PLOT)
%                Default: '-'
% 'LineWidth'    Width of arrow line, in pts (1/72 inch)
%                Positive numeric scalar
%                Default: 0.5
%
% A = ARROW(...) returns a graphics object group A containing the graphics
% objects of the arrow.
%
% Examples
% rng(123)
% x = rand(1,6);
% plot(x)
% arrow([4 x(4)],[1 x(1)],'Linewidth',2)
% arrow([2 x(2)],[6 x(6)],'m--','HeadIndent',-0.2,'HeadWidth',18)

if ishghandle(ax)
    if ~strcmp(ax.Type,'axes')
        s = inputname(1);
        if isempty(s)
            s = 'First input';
        end
        error([s,' is not a valid axes object'])
    end
else
    % Shift inputs if first input isn't an axes object
    if nargin > 2
        varargin = [{endpt},varargin];
    end
    endpt = startpt;
    startpt = ax;
    % Use current axes
    ax = gca;
end

% Get all display options by parsing the extra inputs
opts = parseoptions(varargin,ax);

% Calculate points that define the arrow
[linex,liney,headx,heady] = getpoints(ax,startpt,endpt,opts);

% Draw arrow
l = line(linex,liney,'Parent',ax,...
    'Color',opts.clr,'LineWidth',opts.lwidth,'LineStyle',opts.lstyle);
p = patch(headx,heady,opts.clr,'EdgeColor',opts.clr,'Parent',ax);
% Group line and head together
arr = hggroup('Parent',ax);
l.Parent = arr;
p.Parent = arr;

% Define callbacks to redraw arrow when view changes
% UPDATEPOINTS calls GETPOINTS and updates objects (L and P)
chngfn = @(~,~) updatepoints([],[],l,p,ax,startpt,endpt,opts);

% Invoke UPDATEPOINTS if there is any change to the axis limits...
zoomlistener = addlistener(ax,{'XLim','YLim'},'PostSet',chngfn);
% ...or figure size.
% To avoid stomping on any potential existing window resize callback, keep
% a copy of the current SizeChangedFcn value. Then build a callback with
% multiple actions.
fig = ax.Parent;
existingfn = fig.SizeChangedFcn;
data = addcallbacks(existingfn,arr,ax,startpt,endpt,opts);
fig.SizeChangedFcn = @(obj,evt) multicallbacks(obj,evt,data);
% Define function to remove listeners and callbacks on deletion
arr.DeleteFcn = @(arr,evt) arrowdelete(arr,evt,zoomlistener,fig);

% Force one redraw (in case arrow is put into blank axes)
drawnow
chngfn([],[])

% Return output, if needed
if nargout
    a = arr;
end

function [linex,liney,headx,heady] = getpoints(ax,startpt,endpt,opts)
% Aspect ratio of actual axis lengths
r = ax.PlotBoxAspectRatio(2)/ax.PlotBoxAspectRatio(1);

% Get dx and dy of arrow in relative units
dx = (endpt(1)-startpt(1))/diff(ax.XLim);
dy = (endpt(2)-startpt(2))/diff(ax.YLim);
% Get length of arrow in actual physical units
units = ax.Units;
ax.Units = 'inches';
L = norm([dx*ax.Position(3) dy*ax.Position(4)]);
ax.Units = units;

% For remaining calculations, need to factor in aspect ratio
dy = r*dy;

% Get angle to rotate arrow
theta = atan2(dy,dx);
ct = cos(theta);
st = sin(theta);

% Get the arrowhead dimensions in normalized units (arrow of length 1)
opts.hlength = opts.hlength/(72*L);
opts.hwidth = opts.hwidth/(144*L);
% Scale arrowhead dimensions if necessary
scl = min(0.5/opts.hlength,1);
opts.hlength = opts.hlength*scl;
opts.hwidth = opts.hwidth*scl;
% Define (x,y) locations of 3 points on back of arrowhead (outer corners
% and center of back)
xy = [(1 - opts.hlength) (1 - (1-opts.hindent)*opts.hlength) (1 - opts.hlength);
    opts.hwidth 0 -opts.hwidth];
% Scale and rotate
scl = norm([dx dy]);
M = [ct -st;st ct];
xy = M*xy*scl;
% Translate
x = xy(1,:);
y = xy(2,:);
x = x*diff(ax.XLim) + startpt(1);
y = y*diff(ax.YLim)/r + startpt(2);

% Create arrays for x and y points for line and head of arrow
% Line from startpt to back-center point of arrowhead
linex = [startpt(1) x(2)];
liney = [startpt(2) y(2)];
% Polygon around back of arrowhead to tip (= endpt)
headx = [x endpt(1)];
heady = [y endpt(2)];

function updatepoints(~,~,ln,hd,ax,startpt,endpt,opts)
[linex,liney,headx,heady] = getpoints(ax,startpt,endpt,opts);
ln.XData = linex;
ln.YData = liney;
hd.XData = headx;
hd.YData = heady;

function arrowdelete(arr,~,zoomlistener,fig)
% Delete axis listener
delete(zoomlistener)
% Remove arrow from figure resize function
fn = functions(fig.SizeChangedFcn);
data = fn.workspace{1}.data;
% Find arrow to be deleted
idx = (arr == data.arrow);
% Remove it from callback data
data.arrow(idx) = [];
data.ax(idx) = [];
data.startpt(idx,:) = [];
data.endpt(idx,:) = [];
data.opts(idx) = [];
% Update or reset callback
if numel(data.arrow) == 0
    fig.SizeChangedFcn = data.prev;
else
    fig.SizeChangedFcn = @(obj,evt) multicallbacks(obj,evt,data);
end

function optionstruct = parseoptions(opts,ax)

% Set defaults for properties that can be set by linespec string
% Useful for later checking of doubling up
clr = false;
lstyle = false;
% Set defaults for remaining properties
hlength = 10;
hindent = 0.2;
hwidth = 10;
lwidth = 1;

n = numel(opts);
if mod(n,2)
    % odd number of options => first should be linespec string
    try
        ax2 = axes('Parent',ax.Parent,'Visible','off');
        p = plot(ax2,[NaN NaN],[NaN NaN],opts{1});
        clr = p.Color;
        % If color is all 0/1 values then it's one of the ones with a
        % shorthand string => it was set by linespec string (so leave it).
        % If not, color was not actually set in linespec string
        if ~all(floor(clr)==clr)
            clr = false;
        end
        lstyle = p.LineStyle;
        % linespec string can set marker style, which is irrelevant
        if ~strcmp(p.Marker,'none')
            warning('Marker specification ignored (no markers on arrows)')
        end
        % Get rid of line
        delete(ax2)
    catch
        error('Unrecognized color/linestyle option')
    end
    % Throw away first option
    opts(1) = [];
    n = n-1;
end
% Remaining options are now name-value pairs
n = n/2;
opts = reshape(opts,2,n);
% Names
names = opts(1,:);
if iscellstr(names)
    names = lower(names);
else
    error('Property names must be strings')
end
% Values
values = opts(2,:);

% Check names
allowednames = {'color','headlength','headwidth','headindent','linewidth','linestyle'};
namecheck = find(~ismember(names,allowednames),1,'first');
if ~isempty(namecheck)
    error([names{namecheck},' is not a valid property name'])
end

% Assign values to appropriate properties
for k = 1:n
    switch names{k}
        case allowednames{1}
            if ~islogical(clr)
                % Was previously set (by linespec string)
                warning('Color set by line specification string overwritten by ''Color'' option')
            end
            clr = values{k};
        case allowednames{2}
            hlength = values{k};
        case allowednames{3}
            hwidth = values{k};
        case allowednames{4}
            hindent = values{k};
        case allowednames{5}
            lwidth = values{k};
        case allowednames{6}
            if ~islogical(lstyle)
                % Was previously set (by linespec string)
                warning('Line style set by line specification string overwritten by ''LineStyle'' option')
            end
            lstyle = values{k};
    end
end

% Set defaults for color and linestyle (if not yet changed)
if islogical(clr)
    clr = 'k';
end
if islogical(lstyle)
    lstyle = '-';
end

% If marker is given, but no linestyle (in linespec), we end up with no
% linestyle. This is bad.
if strcmp(lstyle,'none')
    lstyle = '-';
    warning('No linestyle given. Reverting to solid line')
end

% Check values before returning
if ~(ischar(clr) || (isnumeric(clr) && (numel(clr)==3)))
    error('Color must be specified as a string or a 3-element numeric vector')
end
if ~(isnumeric(hlength) && (hlength > 0))
    error('HeadLength must be a positive numeric scalar')
end
if ~(isnumeric(hwidth) && (hwidth > 0))
    error('HeadWidth must be a positive numeric scalar')
end
if ~(isnumeric(hindent) && (hindent >= -1) && (hindent <= 1))
    error('HeadIndent must be a numeric scalar from -1 to 1')
end
if ~(isnumeric(lwidth) && (lwidth > 0))
    error('LineWidth must be a positive numeric scalar')
end
if ~ischar(lstyle)
    error('LineStyle must be specified as a string')
end

optionstruct = struct('clr',clr,'hlength',hlength,'hwidth',hwidth,...
    'hindent',hindent,'lwidth',lwidth,'lstyle',lstyle);

function multicallbacks(obj,evt,data)
% Make a callback out of multiple callbacks, using data structure
% Run initial callback (if any -- prevtype = 0 if none)
f = data.prev;
switch data.prevtype
    case 1
        % function handle
        f(obj,evt);
    case 2
        % cell array
        f{1}(obj,evt,f{2:end});
    case 3
        % string (nonempty)
        eval(f);
end
% Run UPDATEPOINTS on each arrow in turn
for k = 1:length(data.arrow)
    arr = data.arrow(k);
    updatepoints(obj,evt,arr.Children(2),arr.Children(1),data.ax(k),...
        data.startpt(k,:),data.endpt(k,:),data.opts(k))
end

function data = addcallbacks(existingfn,arrow,ax,startpt,endpt,opts)

% Determine if there is already a figure resize callback
if isa(existingfn,'function_handle')
    fn = functions(existingfn);
    if isempty(strfind(fn.function,'multicallbacks')) || isempty(strfind(fn.file,'arrow'))
        % Existing fun handle, but not from here
        data = [];
    else
        % Already have arrow callbacks
        data = fn.workspace{1}.data;
    end
else
    data = [];
end
if isempty(data)
    % First time adding arrow resize callbacks
    data = struct('prev',existingfn,'prevtype',[],'arrow',[],...
        'ax',[],'startpt',[],'endpt',[],'opts',[]);
    % Save initial resize callback information
    if isa(existingfn,'function_handle')
        data.prevtype = 1;
    elseif iscell(existingfn)
        data.prevtype = 2;
    elseif isempty(existingfn)
        data.prevtype = 0;
    else
        data.prevtype = 3;
    end
end
% Add current arrow information to callback data
data.arrow = [data.arrow;arrow];
data.ax = [data.ax;ax];
data.startpt = [data.startpt;startpt(:)'];
data.endpt = [data.endpt;endpt(:)'];
data.opts = [data.opts;opts];
