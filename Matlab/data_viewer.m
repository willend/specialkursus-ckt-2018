% hygge mcstas
close all
clear all

%addpath(genpath('/Users/TummasN/Documents/GitHub/specialkursus-ckt-2018/Matlab'))       %replace with own path 
addpath(genpath('C:\Users\kedde\Documents\GitHub\specialkursus-ckt-2018\Matlab'))
%path = '/Users/TummasN/Documents/GitHub/specialkursus-ckt-2018/Krystal_sim/Analy_crystal_20181105_102130';      %path of simulation
path = 'C:\Users\kedde\Documents\Magnon_rerun_events_B\40';      %path of simulation

cd(path)
d=dir('*.t_y');
N=size(d);

Data = iData(d(1).name);

y = Data.x;
tof = Data.y;

plot(Data)
view(0,90)
xlim([tof(1) tof(end)])
ylim([y(1) y(end)])



