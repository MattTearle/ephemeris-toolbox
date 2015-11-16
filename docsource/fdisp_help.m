%% FDISP  
% Formatted display of angles
%
%% Description
% |fdisp(format, A)| displays the array of angles |A| to the screen,
% according to the |format| string.
%
% |X = fdisp(format, A)| returns the formatted display to a string |X|.
%
% The |format| specification string is |sprintf|-style, but uses
% angle-specific formats: |%d|, |%t|, and |%z|
%
% |%d| = degrees/minutes/seconds: 42.314 -> +42°18'50"
%
% |%t| = time (hours/minutes/seconds): 42.314 -> 02:49:15
%
% |%z| = zodiac (degrees-sign/minutes/seconds): 42.314 -> 12 Tau 18'50"
%
% Other |sprintf| markups such as \t and are valid
%
% A scalar format specifier is applied to all elements. Multiple specifiers
% are applied to the columns of |A|; the number of specifiers must match
% the number of columns of |A|.
%
%% Example:
fdisp('%z\t%2d -- %t',angl([42,3.14,666]))

%% See also
% <angl_help.html |angl|>