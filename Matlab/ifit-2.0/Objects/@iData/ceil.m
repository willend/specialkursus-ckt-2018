function a = ceil(a)
% b = ceil(s) : upper integer round of iData object
%
%   @iData/ceil function to round the elements of 's' to the nearest integers
%   towards plus infinity.
%
% input:  s: object or array (iData)
% output: b: object or array (iData)
% ex:     b=ceil(a);
%
% Version: oct.. 23, 2018
% See also iData, iData/floor, iData/ceil, iData/round

a = iData_private_unary(a, 'ceil');

