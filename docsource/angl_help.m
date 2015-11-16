%% angl
% The |angl| data type
%
%% Description
% The |angl| data type is designed to store angles, including units which
% can be degrees, radians or time (1 hour of angle corresponds to 15
% degrees).
%
% The |angl| data type allows you to perform calculations without having to
% worry about units or modular arithmetic.  You can perform array
% calculations with angles using standard MATLAB notation.
%
%% Creation
% Use the |angl| function to create angles from numeric values or
% variables.
%
% |A = angl(X,U)| creates an angle array |A| with values specified by
% the array |X| and units specified by the string or cell array of
% strings |U|.  If |U| is a single string, all values in |X| will have the
% same units.  If |U| is a cell array of strings, it must have the same
% dimensions as |X|.  Units can be specified as 'degrees', 'radians', or
% 'hours'.
%
% |A = angl(X)| creates an angle array, using default units of 'degrees'.
%
%  Examples:
%%
x = angl(3.14)
y = angl([45,90,200],'degrees')
z = angl([45,6,pi/2],{'degrees','hours','radians'})

%% Methods
% You can perform the following arithmetic operations on |angl| variables:
%
% * Addition
% * Subtraction, which is interpreted as angle difference (ie angle between two polar angles)
% * Negation
% * Scalar multiplication and division
%
% The following mathematical functions are defined for |angl| variables,
% with the same interpretation as for any numeric variable:
%
% * |ceil|
% * |floor|
% * |round|
% * |cos|
% * |sin|
% * <tan_help.html |tan|>
%
% You can use the following functions with the |angl| data type:
%
% * <convert_help.html |convert|> - Change angle units
% * <dec2sex_help.html |dec2sex|> - Return angle as degrees, minutes, and seconds
% * <double_help.html |double|> - Convert |angl| to |double|
% * <fdisp_help.html |fdisp|> - Display angle with formatting
% * <normalize_help.html |normalize|> - Normalize angle to [-360,360] degrees (or equivalent)
% * <normalizeinrange_help.html |normalizeinrage|> - Normalize angle to given range
% * <plot_help.html |plot|> - 
% * <scatter_help.html |scatter|> - 
% * <sexiform_help.html |sexiform|> - 
% * <timeform_help.html |timeform|> - 
% * <zodiacform_help.html |zodiacform|> - 
%
% The following functions return |angl| variables:
%
% * acos
% * acot
% * asin
% * atan
% * <sex2dec_help.html |sex2dec|> - Convert sexigesimal (degrees, minutes, seconds) form to decimal |angl|

