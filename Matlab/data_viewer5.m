% hygge mcstas
close all
%clear all

NuN=NaN;

%addpath(genpath('/Users/TummasN/Documents/GitHub/specialkursus-ckt-2018/Matlab'))       %replace with own path
addpath(genpath('C:\Users\kedde\Documents\GitHub\specialkursus-ckt-2018\Matlab'))
%path = '/Users/TummasN/Documents/GitHub/specialkursus-ckt-2018/Krystal_sim/Analy_crystal_20181105_102130';      %path of simulation
path = 'C:\Users\kedde\Documents\Magnon_rerun_B\';      %path of simulation
path2 = 'C:\Users\kedde\Documents\Magnon_rerun_A\';      %path of simulation

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


Qmatrix = zeros(30,30,256);
dEmatrix = zeros(30,30,256);

%The scattering plane coordinates for the magnon sample.
Hbar = zeros(30,30,256);  %The symmetric Q_H direction with the Miller index H.
Kbar = zeros(30,30,256);  %The symmetric Q_L direction with the Miller index L.

%Empty intensity matrix
In = zeros(30,30,256);

tic
figure(1)
hold on


for nn = 0:29
    
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
            t1(find(I==0))=NuN;
            
            
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
            Hbar_v = Q.*cos(beta).*cos(phi) + Q.*sin(beta).*sin(phi);
            Kbar_v = Q.*cos(beta).*sin(phi) - Q.*sin(beta).*cos(phi);
            %Hbar_v = Q.*sin(TwoTheta + deg2rad(kkk*10) + alpha)'.*cos(deg2rad(2*nn)) + Q.*cos(TwoTheta + deg2rad(kkk*10) + alpha)'.*sin(deg2rad(2*nn));
            %Kbar_v = Q.*sin(TwoTheta + deg2rad(kkk*10) + alpha)'.*sin(deg2rad(2*nn)) - Q.*cos(TwoTheta + deg2rad(kkk*10) + alpha)'.*cos(deg2rad(2*nn));
            
            DeltaE(find(t==0))=NuN;
            Q(find(t==0))=NuN;
            Hbar_v(find(t==0))=NuN;
            Kbar_v(find(t==0))=NuN;
            
%             Hbar_s=zeros(30,256);Kbar_s=zeros(30,256);I_sis=zeros(30,256);E_lem=zeros(30,256);
% 
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
            
            Qmatrix(nn+1,:,:) = Q;
            dEmatrix(nn+1,:,:) = DeltaE;
            Hbar(nn+1,:,:) = Hbar_v;
            Kbar(nn+1,:,:) = Kbar_v;
            
            %     figure(1)
%             plot(Data)
%             view(0,90)
%             xlim([tof(1) tof(end)])
%             ylim([y(1) y(end)])
            %     figure(1)
%             plot(Q,DeltaE,'.','color',[0,nn,nn]*0.02)
%             
%             a = surf(Q,DeltaE,I);
%             shading interp
%             colormap jet
%             colorbar
%             a.EdgeColor='none';
%             xlabel('Q [Å^{-1}]')
%             ylabel('\Delta E [meV]')
%             
            %     xlim([0 1])
            
            %s=surf(Hbar_s,Kbar_s,E_lem);
            %s=surf(Hbar_v,Kbar_v,DeltaE);
            s.EdgeColor='none';
            shading interp
             %slice(In,Hbar_v,Kbar_v,DeltaE)

%             % %         colorMap = [linspace(0,1,256)',zeros(256,2)];
%             %         colormap(colorMap); colorbar;
%             %colormap(flipud(jet))
%             colorbar
%             s.EdgeColor='none';
            %plot(Hbar_v,Kbar_v,'.','color',[0,nn,nn]*0.02)
             xlabel('Q_H [Å^{-1}]')
             ylabel('Q_K [Å^{-1}]')
%             title('Allah u akhbar')
            %         zlabel('\Delta E [meV]')
        end
    end
    
end

%plot(Hbar(nn+1,:,:),Kbar(nn+1,:,:))
%     plot(Qmatrix(40,:,:),dEmatrix(40,:,:),'.','color',[0,0,0])
%     xlabel('Q [Å^{-1}]')
%     ylabel('Energy [meV]')

toc
