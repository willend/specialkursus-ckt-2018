function d = ifitpath
% ifitpath iFit library location
%
% Version: oct.. 23, 2018
% (c) E.Farhi, ILL. License: EUPL.

d = [ fileparts(which('iData/version')) filesep '..' filesep '..' filesep ];


