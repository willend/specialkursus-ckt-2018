% hygge mcstas
close all
clear all

%addpath(genpath('/Users/TummasN/Documents/GitHub/specialkursus-ckt-2018/Matlab'))       %replace with own path
addpath(genpath('C:\Users\kedde\Documents\GitHub\specialkursus-ckt-2018\Matlab'))
%path = '/Users/TummasN/Documents/GitHub/specialkursus-ckt-2018/Krystal_sim/Analy_crystal_20181105_102130';      %path of simulation
path = 'C:\Users\kedde\Documents\Magnon_rerun_events_B\';      %path of simulation

Qmatrix = zeros(46,30);
dEmatrix = zeros(46,30);

for nn = 0:45
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
    
    %t = zeros(30,256);
    %for ii = 1:30
    %    for kk = 1:256
%             if I = 0
%                 continue
%             end
%             t(ii,kk) = tof(kk);
%         end
%     end
%     
%     
%     if max(I)==0
%         Qmatrix(nn+1,:) = Q;
%         dEmatrix(nn+1,:) = DeltaE;
%         continue
%     end
    
    Tmax = zeros(1,30);
    for ii=1:30
        if max(I(ii,:)) == 0
            continue
        end
        i_max= find(I(ii,:) == max(I(ii,:)));
        Tmax(ii)=tof(i_max);
    end
%     
%     
    
    %5 meV short dimensions in SI
    Ef = 5e-3*1.602e-19; %
    mn = 1.627e-27; %kg
    hbar = 1.0546e-34; %J*s
    
    kf = sqrt(2*mn*Ef/hbar^2);
    
    Lsa = 1.544; %distance from sample to analyser
    Lad = 1.623; %distance from analyser to detector
    
    wa = 0.1744; %width of analyser
    wd = 0.2976; %width of detector
    wad = wd*Lsa/(Lsa + Lad);
    
    ya = y*wad/wd;
    
    %Scattering angle
    alpha = atan(ya./Lsa);
    
    Lsaalpha = Lsa./cos(alpha);
    Ladalpha = Lad./cos(alpha);
    
    %Total length from sample to detector
    L2 =abs(Lsaalpha + Ladalpha);
    %The arrival time from sample to detector
    t2 = L2./sqrt(2*Ef/mn);
    
    %The initial arrival time from source to sample
    t1 = Tmax-t2;
    
    %Incoming neutrons wavevector and energy
    ki = mn*162./(hbar*t1);
    Ei = hbar^2*ki.^2/(2*mn);
    
    DeltaE = (Ei - Ef)/(1.602e-19)*1e3;
    
    %The scattering vector size
    Q = sqrt(ki.^2+kf^2-2.*ki*kf.*cos(alpha))./1e10;
    
    DeltaE(find(Tmax==0))=NaN;
    Q(find(Tmax==0))=NaN;
    
    Qmatrix(nn+1,:) = Q;
    dEmatrix(nn+1,:) = DeltaE;
    
 
%     figure(1)
%     plot(Data)
%     view(0,90)
%     xlim([tof(1) tof(end)])
%     ylim([y(1) y(end)])
%     figure(2)
%     plot(Q,DeltaE,'o-')
    %xlim([0 1])
end
plot(Qmatrix,dEmatrix,'o')


