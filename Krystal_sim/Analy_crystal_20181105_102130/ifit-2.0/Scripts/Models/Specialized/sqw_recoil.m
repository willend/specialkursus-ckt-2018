function signal=sqw_recoil(varargin)
% model = sqw_recoil(p, q ,w, {signal}) : Recoil dispersion(Q) for a single harmonic oscillator
%
%   iFunc/sqw_recoil: a 2D S(q,w) with a recoil dispersion, isotropic harmonic oscillator 
%
%   This is a 2D recoil model for a single harmonic oscillator particle with given mass.
%   The mass defines the dispersion energy, while the harmonic oscillator energy
%   defines the Debye-Waller function, i.e. intensity vs Q.
%   This model is also known as 'impulse approximation' or 'short time' Gaussian. 
%
%   The dispersion has the form:
%      S(q,w) = 1/sqrt(2*pi*delta2).*exp( -(w - Er).^2/2./delta2 )
%   where
%      Er     = h2 q2/2/m
%      delta2 = Er.*coth(w0/2/kT)
%
%   You can build a recoil model for a given mass and oscillator energy with syntax:
%      sqw = sqw_recoil(mass)
%   or
%      sqw = sqw_recoil([ mass energy ])
%
%   You can of course tune other parameters once the model object has been created.
%   To model more than one oscillator/particle recoil, just add these models together.
%
%   Evaluate the model with syntax:
%     sqw(p, q, w)
%
% input:  p: sqw_recoil model parameters (double)
%             p(1)= Mass       Mass of the scattering unit [g/mol]
%             p(2)= Ex_Energy  Harmonic Excitation Energy [meV]
%             p(3)= Amplitude 
%             p(4)= Temperature  Temperature [K]  
%         q:  axis along wavevector/momentum in Angs-1 (row,double)
%         w:  axis along energy in meV (column,double)
% output: signal: model value
%
%
% Example:
%   s=sqw_recoil(2); % Deuterium mass
%   plot(iData(s, [], 0:.1:20, (0:100)))  % q=0:20, energy=0:100
%
% Reference: Schober, JNR 17 (2014) 109–357 - DOI 10.3233/JNR-140016
%
% Version: oct.. 23, 2018
% See also iData, iFunc
%   <a href="matlab:doc(iFunc,'Models')">iFunc:Models</a>
% (c) E.Farhi, ILL. License: EUPL.

signal.Name           = [ 'sqw_recoil dispersion(Q) for given mass and excitation energy [' mfilename ']' ];
signal.Description    = 'A 2D S(q,w) with a recoil dispersion.';

signal.Parameters     = {  ...
'Mass             Mass of the scattering unit [g/mol]' ...
'Ex_Energy        Harmonic Excitation Energy [meV]' ...
'Amplitude' ...
'Temperature  Temperature [K]' };
  
signal.Dimension      = 2;         % dimensionality of input space (axes) and result
signal.Guess = [ 1 10 1 10 ];

signal.Expression     = { ...
 'm = p(1); w0= p(2); T = p(4);' ...
 'q = x; w = y;' ...
 'if isvector(q) && isvector(w) && numel(q) ~= numel(w), [q,w] = meshgrid(q,w); end' ...
 'kb      = 1.38064852E-023;  % Boltzmann [J/K]' ...
 'mn      = 1.674927471E-027; % neutron mass [kg]' ...
 'e       = 1.60217662E-019;  % [C]' ...
 'HBAR    = 1.05457168e-34;   % Plank/2PI [m2 kg / s]' ...
 'meVtoK  = e/1000/kb;        % = 11.6045 = [meV] to [K]' ...
 'q2toE   = HBAR*HBAR/2/mn/e*1000*1e20; % = 2.0721 = [Angs^-2] to [meV] ' ...
 'kT      = T/meVtoK;                        % T [K] -> [meV] = 11.6045' ...
 'Er      = q2toE*q.^2/2/m;' ...
 'delta2  = Er.*coth(w0/2/kT);' ...
 'signal  = p(3)./sqrt(2*pi*delta2).*exp(- (w - Er).^2/2./delta2);' ...
 };

signal= iFunc(signal);
signal= iFunc_Sqw2D(signal); % overload Sqw2D flavour

if nargin == 1 && isnumeric(varargin)
  p = varargin{1};
  if numel(p) >= 1, signal.ParameterValues(1) = p(1); end
  if numel(p) >= 2, signal.ParameterValues(2) = p(2); end
end
