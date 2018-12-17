% hygge mcstas
close all
%clear all

%PW=pwd;

NuN=NaN;
matl_path='/Users/TummasN/Documents/GitHub/specialkursus-ckt-2018/Matlab';
fil = matfile('Datafile_B','Writable', true);

addpath(genpath('/Users/TummasN/Documents/GitHub/specialkursus-ckt-2018/Matlab'))       %replace with own path
addpath(genpath(matl_path))
%path = '/Users/TummasN/Documents/GitHub/specialkursus-ckt-2018/Krystal_sim/Analy_crystal_20181105_102130';      %path of simulation
cd(matl_path) 
path = '/Users/TummasN/Documents/GitHub/specialkursus-ckt-2018/BIFROST/Magnon_rerun_events_B_noEfocus_noSqq';      %path of simulation
%path = '/Users/TummasN/Documents/OneDrive - Danmarks Tekniske Universitet/DTU/5. semester/Mcstas/Magnon_rerun_A';      %path of simulation

% path = [PW '/../BIFROST/Magnon_rerun_A_gcc-4.9'];
% path2 = [PW '/../BIFROST/Magnon_rerun_B_gcc-4.9'];


y_matrix=cell(90,9,30);
tof_matrix=cell(90,9,256);
I_matrix=cell(90,9);
    
for iii=0:90
    cd(path)
    cd(num2str(iii))
    d=dir('*.t_y');
    N=size(d);
    
    for kkk = 0:8
        Data = iData(d(kkk+1).name);
        y = Data.x;
        tof = Data.y;
        I = Data.I;
        y_matrix{iii+1,kkk+1}=y;
        tof_matrix{iii+1,kkk+1}=tof;
        I_matrix{iii+1,kkk+1}=I;
    end
end
cd(matl_path)
fil.y=y_matrix;
fil.tof=tof_matrix;
fil.I=I_matrix;
