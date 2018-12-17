% hygge mcstas
close all
%clear all

NuN=NaN;

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

%loading the datafile
%filA=matfile('Datafile_A','Writable', false);
fil=matfile('Datafile_B','Writable', true);
y_mat = fil.y;;
tof_mat = fil.tof;
I_mat = fil.I;


tic
figure(1)
hold on

Qqh=[];
Qqk=[];
Eee=[];
Iii=[];
QQq=[];

Hbar_m=zeros(90,30,256);Kbar_m=zeros(90,30,256);E_m=zeros(90,30,256);II_m=zeros(90,30,256);



for nn = 0:90
    
    for iii=1:1
        
        for kkk = 0:8
            y=y_mat(nn+1,kkk+1);y=y{1};
            tof=tof_mat(nn+1,kkk+1);tof=tof{1};
            I=I_mat(nn+1,kkk+1);I=I{1};
            
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
            
            %In(nn+1,:,:)=I;
            idxt1=find(I==0);
            t1(idxt1)=NuN;
            
            
            %Incoming neutrons wavevector and energy
            ki = mn*162./(hbar.*t1);
            Ei = hbar^2.*ki.^2/(2*mn);
            
            DeltaE = (Ei - Ef)./(1.602e-19)*1e3;
            
            %The scattering vector size
            Q = sqrt(ki.^2+kf^2-2*cos(TwoTheta+alpha+deg2rad(kkk*10)+deg2rad(iii*0.5))'.*ki.*kf)/1e10;
            Q(find(I==0))=NuN;
            
            beta = acos((-kf^2+ki.^2+(Q*1e10).^2)./(2*ki.*Q*1e10));
            two_theta=TwoTheta+alpha+deg2rad(kkk*10)+deg2rad(iii*0.5);
            phi=deg2rad(1*nn);

            %The scattering vectors
            Hbar_v = Q.*cos(beta-phi);
            Kbar_v = Q.*sin(beta-phi);
            
            Hbar_v(idxt1)=0;Kbar_v(idxt1)=0;DeltaE(idxt1)=0;I(idxt1)=0;
            
            Hbar_m(nn+1,:,:)=Hbar_v;Kbar_m(nn+1,:,:)=Kbar_v;E_m(nn+1,:,:)=DeltaE;II_m(nn+1,:,:)=I;

            Qhtmp=reshape(Hbar_v,30*256,1);
            Qktmp=reshape(Kbar_v,30*256,1);
            Etmp=reshape(DeltaE,30*256,1);
            Itmp=reshape(I,30*256,1);
                        
            idxNan=isnan(Qhtmp);
            Qhtmp(idxNan)=0;
            Qktmp(idxNan)=0;
            Etmp(idxNan)=0;
            Itmp(idxNan)=0;            
           
            Qqh=[Qqh;Qhtmp];
            Qqk=[Qqk;Qktmp];
            Eee=[Eee;Etmp];
            Iii=[Iii;Itmp];
            QQq=[QQq;Qhtmp.*Qktmp]; 
           
            
%             
%             a = surf(Q,DeltaE,I);
%             shading interp
%             colormap jet
%             colorbar
%             a.EdgeColor='none';
%             xlabel('Q [�^{-1}]')
%             ylabel('\Delta E [meV]')
%             
            nn
            kkk
        end
    end
end

toc

fil.I=Iii;fil.QH=Qqh;fil.dE=Eee;fil.QK=Qqk;


%cd(PW)

% Define scattered interpolant in collected sparse dataset



% F = scatteredInterpolant(Qqh,Qqk,Eee,Iii)
% [Qa,Qb,eE]=meshgrid(linspace(0.5,3.5,100),linspace(-2,2,100),linspace(1,5,100));
% iI=F(Qa,Qb,eE);
% 
% save workspace