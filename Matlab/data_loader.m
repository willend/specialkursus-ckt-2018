close all
clear all

%load data
matl_path='C:\Users\kedde\Documents\GitHub\specialkursus-ckt-2018\Matlab';
cd(matl_path);
data=matfile('data_fil60.mat','Writable',false);

Kbarf = data.Kbar;
Hbarf = data.Hbar;
dEf = data.dEmatrix;
If = data.I_matrix;

figure
hold on
for ii=1:5
    ii
    Kbar=Kbarf{ii};
    Hbar=Hbarf{ii};
    dE=dEf{ii};
    %     Kbar=reshape(data.Kbar(ii,:,:),30,256);
    %     Hbar=reshape(data.Hbar(ii,:,:),30,256);
    %     dEmatrix=reshape(data.dEmatrix(ii,:,:),30,256);
    %     I=reshape(data.I_matrix(ii,:,:),30,256);
    surf(Hbar,Kbar,dE)
%     s.EdgeColor='none';
    shading interp
    xlabel('Q_H [�^{-1}]')
    ylabel('Q_K [�^{-1}]')
end