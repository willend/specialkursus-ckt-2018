function a = atanh(a)
% b = atanh(s) : computes the inverse hyperbolic tangent of iFunc object
%
%   @iFunc/atanh function to compute the inverse hyperbolic tangent of data sets.
%
% input:  s: object or array (iFunc)
% output: b: object or array (iFunc)
% ex:     b=atanh(a);
%
% Version: oct.. 23, 2018
% See also iFunc, iFunc/cos, iFunc/acos, iFunc/sin, iFunc/asin, iFunc/tan, iFunc/atan

a = iFunc_private_unary(a, 'atanh');

