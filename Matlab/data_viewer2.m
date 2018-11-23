% hygge mcstas
close all
clear all

%addpath(genpath('/Users/TummasN/Documents/GitHub/specialkursus-ckt-2018/Matlab'))       %replace with own path 
addpath(genpath('C:\Users\kedde\Documents\GitHub\specialkursus-ckt-2018\Matlab'))
%path = '/Users/TummasN/Documents/GitHub/specialkursus-ckt-2018/Krystal_sim/Analy_crystal_20181105_102130';      %path of simulation
path = 'C:\Users\kedde\Documents\Magnon_rerun_events_B';      %path of simulation

cd(path)
d = dir(fullfile(path,'**\analyser1_tmon.t_y'));

for idx = 1:45
    currentFile = fullfile(d(idx).folder);
end

%d=dir('*.t_y');
N=size(d);

Data = iData(d.folder);

y = Data.x;
tof = Data.y;

%5 meV short dimensions
Ef = 5e-3*1.602e-19; %
mn = 1.627e-27; %kg
hbar = 1.0546e-34; %J*s

kf = sqrt(2*mn*Ef/hbar^2);

Lsa = 1.544;
Lad = 1.623;

wa = 0.2976;
wd = 0.2976;
wad = wd*Lsa/(Lsa + Lad);

ya = y.*Lsa/(Lsa+Lad);

if abs(ya)>= -wd/2 & abs(ya) <= wd/2
    ya = ya;
else
    ya = [];
end

%Scattering angle
alpha = atan(ya./Lsa)*180/pi();

Lsaalpha = cos(alpha).*Lsa;
Ladalpha = cos(alpha).*Lad;

%Total length from sample to detector
L2 =abs(Lsaalpha + Ladalpha);
%The arrival time from sample to detector
t2 = L2./sqrt(2*Ef/mn);

%The initial arrival time from source to sample
t1 = tof(1:30)-t2;

%Incoming neutrons wavevector and energy
ki = mn*162./(hbar.*t1);
Ei = hbar^2.*ki.^2/(2*mn);

DeltaE = (Ei - Ef)/(1.602e-19)*1e3;

%The scattering vector size
Q = sqrt(ki.^2+kf^2-2.*ki.*kf.*cos(alpha.*pi()/180))./1e10;

figure(1)
plot(Data)
view(0,90)
xlim([tof(1) tof(end)])
ylim([y(1) y(end)])
%figure(2)
%plot(Q,DeltaE)




