function a = flipud(a)
% b = flipud(s) : Flip object in up/down direction
%
%   @iFunc/flipud function to flip object in up/down direction
%     With 2D data sets, the Y axis (vertical) is inverted.
%
% input:  s: object or array (iFunc)
% output: b: object or array (iFunc)
% ex:     b=flipud(a);
%
% Version: oct.. 23, 2018
% See also iFunc, iFunc/fliplr, fliplr, iFunc/flipud, flipud

a = iFunc_private_unary(a, 'flipud');
