function y=size(s, dim)
% size(s) : get iData object size (number of elements)
%
%   @iData/size function to get iData object size
%
% input:  s: object or array (iData)
%         dim: optional dimension/rank to inquire
% output: v: size of the iData Signal (double)
% ex:     size(iData), size(iData,1)
%
% Version: oct.. 23, 2018
% See also iData, iData/disp, iData/get, length

% EF 23/09/07 iData implementation

if numel(s) > 1  % this is an array of iData
  if nargin > 1, y = builtin('size', s, dim);
  else           y = builtin('size', s); end
  return
end

if numel(s) == 0
  y=0;
else
  y = size(subsref(s,struct('type','.','subs','Signal')));
  if nargin > 1, 
    if dim > length(y), y=0; 
    else
      y = y(dim); 
    end
  end
end

