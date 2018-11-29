% hygge mcstas
close all
clear all

%addpath(genpath('/Users/TummasN/Documents/GitHub/specialkursus-ckt-2018/Matlab'))       %replace with own path
addpath(genpath('C:\Users\kedde\Documents\GitHub\specialkursus-ckt-2018\Matlab'))
%path = '/Users/TummasN/Documents/GitHub/specialkursus-ckt-2018/Krystal_sim/Analy_crystal_20181105_102130';      %path of simulation
path = 'C:\Users\kedde\Documents\Magnon_rerun_events_B\';      %path of simulation

%5 meV short dimensions in SI
Ef = 5e-3*1.602e-19; %
mn = 1.627e-27; %kg
hbar = 1.0546e-34; %J*s

kf = sqrt(2*mn*Ef/hbar^2);
%The first analyser crystal is positioned at a scattering angle of 40
%degrees relative to the sample
TwoTheta = 40*pi/180;

Lsa = 1.544; %distance from sample to analyser
Lad = 1.623; %distance from analyser to detector

wa = 0.1744; %width of analyser
wd = 0.2976; %width of detector
wad = wd*Lsa/(Lsa + Lad);


Qmatrix = zeros(46,30,256);
dEmatrix = zeros(46,30,256);

%The scattering plane coordinates for the magnon sample.
Hbar = zeros(46,30,256);  %The symmetric Q_H direction with the Miller index H.
Kbar = zeros(46,30,256);  %The symmetric Q_L direction with the Miller index L.

figure(1)
hold on
for nn = 40
    cd(path)
    cd(num2str(nn))
    %d = dir(fullfile(path,'**\analyser1_tmon.t_y'));
    
    %for idx = 1:45
    %    currentFile = fullfile(d(idx).folder);
    %end
    d=dir('**\analyser1_tmon.t_y');
    
    N=size(d);
    
    Data = iData(d.name);
    
    y = Data.x;
    tof = Data.y;
    I = Data.I;
    %I(I==0)=NaN;
    
    ya = y*wad/wd;
    
    %Scattering angle
    alpha = atan(ya./Lsa);
    
    Lsaalpha = Lsa./cos(alpha);
    Ladalpha = Lad./cos(alpha);
    
    %Total length from sample to detector
    L2 =abs(Lsaalpha + Ladalpha);
    %The arrival time from sample to detector
    t2 = L2./sqrt(2*Ef/mn);
    
    t = zeros(30,256);
    for ii = 1:30
       for kk = 1:256
            if I(ii,kk) == 0
                continue
            end
            t(ii,kk) = tof(kk);
       end
    %The initial arrival time from source to sample
        t1(ii,:) = t(ii,:)-t2(ii);
    end
    
    %Incoming neutrons wavevector and energy
    ki = mn*162./(hbar.*t1);
    Ei = hbar^2.*ki.^2/(2*mn);
    
    DeltaE = (Ei - Ef)./(1.602e-19)*1e3;
    
    %The scattering vector size
    Q = sqrt(ki.^2+kf^2-2*cos(TwoTheta+alpha)'.*ki.*kf)./1e10;
    
    %The scattering vectors 
    %Hbar_v = ((ki - kf.*cos(TwoTheta + alpha)').*cos(deg2rad(nn)) + kf.*sin(TwoTheta+alpha)'.*sin(deg2rad(nn)))./1e10; 
    %Kbar_v = ((ki - kf.*cos(TwoTheta + alpha)').*sin(deg2rad(nn)) - kf.*sin(TwoTheta+alpha)'.*cos(deg2rad(nn)))./1e10;
    
    Hbar_v = Q.*sin(TwoTheta + alpha)'.*cos(deg2rad(2*nn)) + Q.*cos(TwoTheta + alpha)'.*sin(deg2rad(2*nn));
    Kbar_v = Q.*sin(TwoTheta + alpha)'.*sin(deg2rad(2*nn)) - Q.*cos(TwoTheta + alpha)'.*cos(deg2rad(2*nn));
    
    DeltaE(find(t==0))=NaN;
    Q(find(t==0))=NaN;
    Hbar_v(find(t==0))=NaN;
    Kbar_v(find(t==0))=NaN;

    Qmatrix(nn+1,:,:) = Q;
    dEmatrix(nn+1,:,:) = DeltaE;
    Hbar(nn+1,:,:) = Hbar_v;
    Kbar(nn+1,:,:) = Kbar_v;
 
%     figure(1)
%     plot(Data)
%     view(0,90)
%     xlim([tof(1) tof(end)])
%     ylim([y(1) y(end)])
%     figure(1)
%     plot(Q,DeltaE,'.','color',[0,nn,nn]*0.02)
%     xlabel('Q [Å^{-1}]')
%     ylabel('Energy [meV]')
%     xlim([0 1])
    
    plot(Hbar_v,Kbar_v,'.','color',[0,nn,nn]*0.02)
    xlabel('Q_H [Å^{-1}]')
    ylabel('Q_K [Å^{-1}]')
end

%plot(Hbar(nn+1,:,:),Kbar(nn+1,:,:))
%     plot(Qmatrix(40,:,:),dEmatrix(40,:,:),'.','color',[0,0,0])
%     xlabel('Q [Å^{-1}]')
%     ylabel('Energy [meV]')


