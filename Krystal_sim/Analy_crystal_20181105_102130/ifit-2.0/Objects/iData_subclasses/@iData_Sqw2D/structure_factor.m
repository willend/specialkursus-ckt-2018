function s = structure_factor(s)
% iData_Sqw2D: structure_factor: compute the structure factor
%  The structure factor is the integral of the dynamic structure factor along
%  the energy axis. It is representative of the material structure.
%  Its Fourier transform is the pair distribution function g(r).
%
%  This function is basically a call to trapz(s,1)
%
% References: Fischer, Barnes and Salmon, Rev Prog Phys 69 (2006) 233
%
% syntax:
%   sq = structure_factor(s)
%
% input:
%   s:  Sqw data set e.g. 2D data set with w as 1st axis (rows, meV), q as 2nd axis (Angs-1).
% output:
%   sq: S(q) data set
%
% Example: sq = structure_factor(iData_Sqw2D('SQW_coh_lGe.nc'));
%
% See also: trapz, iData/trapz
% (c) E.Farhi, ILL. License: EUPL.

  if nargin == 0, return; end

  % handle array of objects
  if numel(s) > 1
    for index=1:numel(s)
      s(index) = feval(mfilename, s(index));
    end
    return
  end
  
  s = Sqw_check(s);
  if isempty(s), return; end

  % compute integral
  s = iData(trapz(s,1)); % on q
  
  % reset axes
  title(s,'S(q)');
  if isempty(s.Label), s.Label='S(q)'; end
  
  if nargout == 0
    fig=figure; 
    plot(s);
    set(fig, 'NextPlot','new');
  end
