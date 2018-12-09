% hygge mcstas
close all
%clear all

%PW=pwd;

NuN=NaN;
matl_path='/Users/TummasN/Documents/GitHub/specialkursus-ckt-2018/Matlab';
addpath(genpath('/Users/TummasN/Documents/GitHub/specialkursus-ckt-2018/Matlab'))       %replace with own path
addpath(genpath(matl_path))
%path = '/Users/TummasN/Documents/GitHub/specialkursus-ckt-2018/Krystal_sim/Analy_crystal_20181105_102130';      %path of simulation
cd(matl_path) 
data=matfile('data_fil60.mat','Writable',true);
path = '/Users/TummasN/Documents/GitHub/specialkursus-ckt-2018/BIFROST/Magnon_rerun_events_A_noEfocus_noSqq';      %path of simulation
%path = '/Users/TummasN/Documents/OneDrive - Danmarks Tekniske Universitet/DTU/5. semester/Mcstas/Magnon_rerun_A';      %path of simulation

% path = [PW '/../BIFROST/Magnon_rerun_A_gcc-4.9'];
% path2 = [PW '/../BIFROST/Magnon_rerun_B_gcc-4.9'];

%5 meV short dimensions in SI
Ef = 5e-3*1.602e-19; %
mn = 1.627e-27; %kg
hbar = 1.0546e-34; %J*s

kf = sqrt(2*mn*Ef/hbar^2);
%The first analyser crystal is positioned at a scattering angle of 40
%degrees relative to the sample
TwoTheta = 40*pi/180;

%Lsa = 1.544; %distance from sample to analyser
Lsa = [1.544 1.623 1.701];
Lad = 1.623; %distance from analyser to detector

%wa = [0.1744 0.1816 0.1886]; %width of analyser
%wd = 0.2976; %width of detector
wd = [0.2976 0.3048 0.3119];
wad = wd.*Lsa./(Lsa + Lad);

%Empty intensity matrix
In = zeros(30,30,256);

tic
figure(1)
hold on

Qqh=[];
Qqk=[];
Eee=[];
Iii=[];

for nn = 0:90
    
    for iii=0:0
        if mod(iii,2)==0
            path_i=path;
        else
            path_i=path;
        end
        cd(path_i)
        cd(num2str(nn))
        %d = dir(fullfile(path,'**\analyser1_tmon.t_y'));
        d=dir('*.t_y');
        N=size(d);
        
        
        for kkk = 0:8
            %         if mod(kkk,3)+1==3
            %             continue
            %         end
            Data = iData(d(kkk+1).name);
            y = Data.x;
            tof = Data.y;
            I = Data.I;
            
            Lsa_l=Lsa(mod(kkk,3)+1);
            wd_l=wd(mod(kkk,3)+1);
            wad_l=wad(mod(kkk,3)+1);
            %I(I==0)=NuN;
            
            ya = y.*wad_l./wd_l;
            
            %Scattering angle
            alpha = atan(ya./Lsa_l);
            
            Lsaalpha = Lsa_l./cos(alpha);
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
            
            %Intensity
            %     I(find(t==0))=NuN;
            %I(I == 0) = NuN;
            In(nn+1,:,:)=I;
            idxt1=find(I==0);
            t1(idxt1)=NuN;
            
            
            %Incoming neutrons wavevector and energy
            ki = mn*162./(hbar.*t1);
            Ei = hbar^2.*ki.^2/(2*mn);
            
            DeltaE = (Ei - Ef)./(1.602e-19)*1e3;
            
            %The scattering vector size
            Q = sqrt(ki.^2+kf^2-2*cos(TwoTheta+alpha+deg2rad(kkk*10)+deg2rad(iii*5))'.*ki.*kf)./1e10;
            Q(find(I==0))=NuN;
            
            beta = acos((-kf^2+ki.^2+(Q*1e10).^2)./(2*ki.*Q*1e10));
            two_theta=TwoTheta+alpha+deg2rad(kkk*10)+deg2rad(iii*5);
            phi=deg2rad(2*nn);
            
%             Hbar_v=cos(phi).*Q.*cos(beta)-Q.*sin(phi).*sin(beta);
%             Kbar_v=-Q.*sin(phi).*cos(beta)-Q.*cos(phi).*sin(beta);
            
            %The scattering vectors
%             Hbar_v = ((ki - kf.*cos(TwoTheta + deg2rad(kkk*10) + deg2rad(iii*5) + alpha)').*cos(deg2rad(nn*2)) - kf.*sin(TwoTheta + deg2rad(kkk*10)+deg2rad(iii*5) + alpha)'.*sin(deg2rad(nn*2)))./1e10;
%             Kbar_v = ((ki - kf.*cos(TwoTheta + deg2rad(kkk*10) +deg2rad(iii*5)+ alpha)').*sin(deg2rad(nn*2)) + kf.*sin(TwoTheta + deg2rad(kkk*10)+deg2rad(iii*5) + alpha)'.*cos(deg2rad(nn*2)))./1e10;
            Hbar_v = Q.*cos(beta-phi);
            Kbar_v = Q.*sin(beta-phi);

%             Hbar_s=zeros(30,256);Kbar_s=zeros(30,256);I_sis=zeros(30,256);E_lem=zeros(30,256);
% 
            %selecting E values         wtf? 
%             vals = [2.7 2.8];
%             Hbar_s(find(DeltaE>vals(1) & DeltaE<vals(2)))=Hbar_v(find(DeltaE>vals(1) & DeltaE<vals(2)));
%             Kbar_s(find(DeltaE>vals(1) & DeltaE<vals(2)))=Kbar_v(find(DeltaE>vals(1) & DeltaE<vals(2)));
%             %I_sis(find(DeltaE>vals(1) & DeltaE<vals(2)))=I(find(DeltaE>vals(1) & DeltaE<vals(2)));
%             E_lem(find(DeltaE>vals(1) & DeltaE<vals(2)))=DeltaE(find(DeltaE>vals(1) & DeltaE<vals(2)));
% 
%             Hbar_s(find(t==0))=NuN;
%             Kbar_s(find(t==0))=NuN;
            
            %         N=size(I_sis);
            %         if N(1)==0
            %             continue
            %         end
            
            
            Hbar_bob{nn+1}=Hbar_v;
            Kbar_bob{nn+1}=Kbar_v;
            I_bob{nn+1}=I;
            Q_bob{nn+1}=Q;
            deltaE_bob{nn+1}=DeltaE;

%             a = surf(Q,DeltaE,I);
%             shading interp
%             colormap jet
%             colorbar
%             a.EdgeColor='none';
%             xlabel('Q [�^{-1}]')
%             ylabel('\Delta E [meV]')
%             
            %     xlim([0 1])
            
            %s=surf(Hbar_s,Kbar_s,E_lem);
            
            s=surf(Hbar_v,Kbar_v,DeltaE,I);
            s.EdgeColor='none';
            shading interp
            xlabel('Q_H [Å^{-1}]')
            ylabel('Q_K [Å^{-1}]')
            
            Qhtmp=reshape(Hbar_v,30*256,1);
            Qktmp=reshape(Kbar_v,30*256,1);
            Etmp=reshape(DeltaE,30*256,1);
            Itmp=reshape(I,30*256,1);
                        
            idxNan=isnan(Qhtmp);
            Qhtmp(idxNan)=[];
            Qktmp(idxNan)=[];
            Etmp(idxNan)=[];
            Itmp(idxNan)=[];            
           
            Qqh=[Qqh;Qhtmp];
            Qqk=[Qqk;Qktmp];
            Eee=[Eee;Etmp];
            Iii=[Iii;Itmp];
        end
    end
end

toc
cd(PW)
data.Hbar = Hbar_bob;
data.Kbar = Kbar_bob;
data.dEmatrix = deltaE_bob;
data.I_matrix = I_bob;
data.Qmatrix = Q_bob;


% Define scattered interpolant in collected sparse dataset

F = scatteredInterpolant(Qqh,Qqk,Eee,Iii)
[Qa,Qb,eE]=meshgrid(linspace(0.5,3.5,100),linspace(-2,2,100),linspace(1,5,100));
iI=F(Qa,Qb,eE);

save workspace