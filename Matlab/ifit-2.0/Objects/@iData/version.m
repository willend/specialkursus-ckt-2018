function [b, vers] = version(a,long_request)
% v = version(iData): iData class version
%
%   @iData/version: version of the iData class library
%
% Version: oct.. 23, 2018

% EF 23/09/07 iData impementation

vers = '2.0';
date = 'oct.. 23, 2018';
auth = 'E.Farhi, P. Willendrup and Y.Debab, (c) ILL/DS/CS <farhi@ill.eu> EUPL';
contrib = 'Eric Ludlam, Felix Morsdorf, Joe Conti, Douglas M. Schwarz, Alexandros Leontitsis, F. Sigworth, Argimiro R. Secchi, Sheela V. Belur, Javad Ivakpour, Nikolaus Hansen, Alexei Kuntsevich and Franz Kappel, C.T. Kelley, Brecht Donckels, Miroslav Balda, Paul Spencer, Juerg Schwizer, Petr Mikulik, David Gingras, Joachim Vandekerckhove, Yi Cao, Oliver Bunk, R. G. Abraham, Bruno Luong, J. Ollivier, D. Riley, E. Trautman, F. Esmonde-White, Y. Altman, Dirk-Jan Kroon, G. Toombes, A. Zheludev, A. Tennant and D. Mc Morrow, Hargreave and Hullah, A. Bouvet/A. Filhol, K. Yamaguchi, W. Falkena, J. Kohlbrecher, J. Rodriguez-Carvajal, W. Baumeister, A. Schmolck, V. Rathod, D. Valevski, J. Almeida, M. Nilsson, J. van Beek, D. Garcia, A. Grinsted, C. Kothe, J. Bialek, D-J Kroon, G. Flandin, M. A. Hopcroft, J. Hokanson, A. Silakov, M. Radin, J. Kantor, C. Pelizzari, J. Yeh, J. Dillon, L. Harriger, G. Romano, Y. Lengwile, W. Falkena, T. Perring, R. Ewings, A. Buts, J. van Duijn, I. Bustinduy, D. Whittaker, S. Toth, D. Alfe, C. Rossant, D. Schwarz, P. Maher, J. Dinale, S. Petit';

b = [ vers ' iFit/iData (' date ') by ' auth '. oct.. 23, 2018' ];
if nargin > 1
  b = [ b '** Licensed under the EUPL V.1.1 ** Contributions from ' contrib '. Send email to <ifit-users@mccode.org> to report bugs and requests. More on <ifit.mccode.org>.' ];
else
  
end

