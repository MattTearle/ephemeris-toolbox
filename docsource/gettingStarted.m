%% Getting Started with Ephemeris Toolbox
%
%% What is an Ephemeris?
% From the Greek for "journal/diary", an ephemeris is a table of the
% positions of celestial objects at a given time.
%
% Ephemeris Toolbox is a collection of MATLAB functions and applications
% related to earth-based observational astronomy
%
%% Applications
% There are three GUI-based applications available in Ephemeris Toolbox
% 
% # Ephemeris Table Calculator
% # Planisphere
% # Star Chart Creator
% 
% The Ephemeris Table Calculator generates a table of the locations of the
% Sun and planets (and Pluto!) for a given time

%% Working with Angles
% Because celestial locations are given in terms of angles, Ephemeris
% Toolbox provides a custom |angl| data type.  You do not need to use
% |angl| variables for most of the functions in Ephemeris Toolbox.
% However, because angles are commonly specified in various different unit
% systems in astronomical applications, using |angl| variables removes
% ambiguity and ensures consistency.  There are also conversion functions
% built into the |angl| data type, including utilities for switching
% between decimal (3.14) and sexigesimal (3°08'24") notation.
%
% Angles can be specified in degrees (default), radians, or hours:

%%
x = angl(3.14)
y = angl(12.34,'hours')

%%
% Angles can be arrays

z = angl(magic(4))

%%
% Mathematical operations are performed naturally, in standard MATLAB
% style, with unit conversions performed automatically as needed:

x+2*y
y+z

%%
% Subtraction is interpreted as the geometric angle between two polar
% angles

x = angl(350);
y = angl(10);
z = x-y
z = y-x

%%
% Numeric subtraction can be achieved using addition and unary minus:

z = x + (-y)
z = y + (-x)

%%
% Angles are normalized to [-360,360] degrees (or equivalent in rad/hrs).
% You can normalize to other standard ranges:

normalizeinrange(z,0) % 0 => [-180,180]
normalizeinrange(z,1) % 1 => [0,360]

%%
% You can create |angl| variables from sexigesimal values:

x = angl.sex2dec(-42,12,34)

%% Specifying Time
