% hygge mcstas
close all
clear all

addpath(genpath('/Users/TummasN/Documents/GitHub/specialkursus-ckt-2018/Krystal_sim/Analy_crystal_20181105_102130'))
path = '/Users/TummasN/Documents/GitHub/specialkursus-ckt-2018/Krystal_sim/Analy_crystal_20181105_102130';
cd(path)
d=dir('*.dat');
N=size(d);

Data = iData(d(1).name);




