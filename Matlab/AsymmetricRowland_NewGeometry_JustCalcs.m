function [AnaArea, Accept, Acc] = AsymmetricRowland_NewGeometry_JustCalcs(En, width, L1, L2, delta, Nbla, Qflag)
close all;

%Program to calculate the analyzer geometries for Bifrost using the
%that vertical Q-acceptance should be constant
%As a first approximation, input is distance from sample and and energy

Acc = 2 %vertical angular acceptance
d = 3.355; %d-spacing of HOPG
%En is the analyzer energy
%width is angular width of the analyser blades
%L1 is the sample analyzer distance
%L2 is the analyzer detector distance
%delta is the spacing between blades (this might prove important via the mounting method)
%NBla is the

if Qflag == 1
    ki = sqrt(En/2.072);
    ki_norm = sqrt(2.7/2.072);
    renoko = (ki_norm/ki)*tand(2);
    Acc = atand(renoko);
end





Ana1 = [L1, En, width];

%i = 1;
%figure('position',[100 100 1300 500], 'color', 'w');

samp = [0 0];
Pa = [Ana1(1) 0]; % Set analyzer point
a5 = asind(9.0446/(sqrt(Ana1(2))*2*d))
Pd = [(L1 +L2*cosd(2*a5)) sind(-2*a5)*L2] % Calculate detector position in symmetric Rowland


[RowCen RowR] = calc_circle(samp, Pa, Pd); %Calculate the Rowland circle
th = 0:pi/7200:2*pi; %Mesh 
xunit = RowR * cos(th) + RowCen(1); %Coordinates for plot
yunit = RowR * sin(th) + RowCen(2); %Coordinates for plot

%subplot(1, 2,1)
%h = plot(xunit, yunit, 'color', 'k', 'linewidth', 2); % This is just plotting the Rowland circle


%hold on
x = [0 Pa(1) Pd(1)]; %x-coordinates of the three relevant points
y = [0 Pa(2) Pd(2)]; %y-coordinates of teh three relevant points   

%plot(x, y, 'o', 'markersize',9, 'color','k')

%line([0 Pa(1)], [0 0], 'color', 'r', 'linewidth', 2) %Line connecting sample and analyzer
%line([Pa(1) Pd(1)], [Pa(2) Pd(2)], 'color', 'r', 'linewidth', 2) %Line connecting analyzer and detector
PPG = [Ana1(1)*(1 - cosd(90-a5)) -sind(90-a5)*Ana1(1)]; %Calculate arbitrary point along this 
%line([Pa(1) PPG(1)], [Pa(2) PPG(2)], 'color', 'b', 'linewidth', 1.5, 'linestyle', '--')%Draw line perpendicular to PG scattering planes

%Below the plot is adjusted
%title('Rowland circle', 'fontname', 'times', 'fontsize', 15);
%set(gca, 'fontname', 'times', 'fontsize', 10);
%xlabel('Horizontal x-axis (along Q-channel) [mm]',  'fontname', 'times', 'fontsize', 12);
%ylabel('Vertical z-axis [mm]', 'fontname',  'times', 'fontsize', 12);
%set(gca, 'xlim', [-500 2500]);
%daspect([1 1 1])

%Old sloppy calculation of converage points on Rowland
%yR = [sind(Acc)*Ana1(1) -sind(Acc)*Ana1(1)];
%xR = [((sqrt(RowR^2 - (sind(Acc)*Ana1(1) - RowCen(2))^2)) + RowCen(1)) ((sqrt(RowR^2 - (-sind(Acc)*Ana1(1) - RowCen(2))^2)) + RowCen(1))];
%The calculation above was done merely assuming a symmetric geometry. We
%need to re-calculate a reference line following 

%Now we do it properly
mu = tand(Acc) %Upward slope
md = tand(-Acc) %Udownward slope

iko1 = RowR^2 * (1 + mu^2) - (RowCen(2) - mu*RowCen(1))^2;
iko2 = RowR^2 * (1 + md^2) - (RowCen(2) - md*RowCen(1))^2;

xR(1)= (RowCen(1) + RowCen(2)*mu + sqrt(iko1))/(1+mu^2)

xR(2)= (RowCen(1) + RowCen(2)*md + sqrt(iko2))/(1+md^2)

yR(1) = mu*xR(1);
yR(2) = md*xR(2);




%plot(xR, yR, 'color', 'b', 'linewidth', 3)

%For this, remember:
%delta is the spacing between blades (this might prove important via the mounting method)
%NBla is the

% nspace = Nbla-1;
% %First calculate the conversion between circle segments and y-value for center blade
% xcc = Ana1(1) - RowCen(1); 
% ycc = 0 - RowCen(2);
% Radc = atan(ycc/xcc); %angle of point on Rowland circle
% fraccc = (sin(Radc+0.01)-sin(Radc))/0.01;
% totspace = nspace*delta*fraccc;
% dy = (2*sind(Acc)*Ana1(1)-totspace)/Nbla; %The minus on is to leave space between blades
% 

%plot(RowCen(1), RowCen(2), 'o', 'linewidth', 2)
%line([RowCen(1) xR(1)], [RowCen(2) yR(1)], 'color', 'k')
%line([RowCen(1) xR(2)], [RowCen(2) yR(2)], 'color', 'k')

%Calculate angle span of analyser on Rowland circle
AngR1 = atand((xR(1) - RowCen(1))/(yR(1) - RowCen(2)));
AngR2 = atand((xR(2) - RowCen(1))/(yR(2) - RowCen(2)));

angdif = abs(AngR1-AngR2); %Ok

circum = 2*pi*RowR; %Circumference of Rowland circle - ok

dangle = (delta/circum)*360; %Angular equivalent of blade spacing on Rowland circle - ok


nspace = Nbla-1;

totspace = nspace*dangle; % ok

da = (angdif-totspace)/Nbla;

AngBla = [];

AngBla(1) = AngR1 + da/2;

xA(1) = sind(AngBla(1))*RowR + RowCen(1);
yA(1) = cosd(AngBla(1))*RowR + RowCen(2);

 %  plot(xA(i), yA(i), 's')
   
%Here we calculate the center positions - be smarter about this! Take into
%account the relationship |sin(ph1) - sin(phi2)| / |phi1-phi2|
%Old approach below
%totspace = nspace*delta;
%dy = (2*sind(Acc)*Ana1(1)-totspace)/Nbla; %The minus on is to leave space between blades

for i = 2:Nbla
    AngBla(i) = AngBla(i-1) + da + dangle;
    xA(i) = sind(AngBla(i))*RowR + RowCen(1);
    yA(i) = cosd(AngBla(i))*RowR + RowCen(2);
  %  plot(xA(i), yA(i), 's')
end

Segleng = (da/360)*circum


Wbla = Segleng; % Width of the blades
%Wbla = Segleng; % Width of the blades of no rounding is made

A_a5 = [];
da5 = [];

for j = 1:Nbla
    d1 = sqrt((xA(j)-Ana1(1))^2 + yA(j)^2);
    dd1 = acosd(1-d1^2/(2*RowR^2));
    da5(j) = dd1/2;
    if j > Nbla/2
        da5(j) = -dd1/2;
    end
    A_a5(j) =  a5-da5(j);
end
    

%for k = 1:Nbla
   % lina1 = line([0 xA(k)], [0 yA(k)], 'color', 'r', 'linestyle', '--') %Drawing sample-analyzer line for blade in question
 %   leng = sqrt(xA(k)^2 + yA(k)^2); %Defining the length of line
    %Below we draw the corresponding scattering line down to the middle of the
    %detector
   % line1b = line([xA(k) (xA(k)+cosd(-2*a5+da5(k))*leng*1.2)],[yA(k) (yA(k)+sind(-2*a5+da5(k))*leng*1.2)], 'color', 'r', 'linestyle', '--')
%end



%subplot(1, 2,2)
%h = plot(xunit, yunit, 'color', 'k', 'linewidth', 1.3, 'linestyle', '--'); 
%hold on%
%plot(x, y, 'o', 'markersize',9, 'color','k')
%line([0 Pa(1)], [0 0], 'color', 'r', 'linewidth', 2)
%line([Pa(1) Pd(1)], [Pa(2) Pd(2)], 'color', 'r', 'linewidth', 2)
PPG = [Ana1(1)*(1 - cosd(90-a5)) -sind(90-a5)*Ana1(1)];
%plot(xA, yA, 'o', 'markersize',6, 'color', 'b')
%hold on

%for i = 1:Nbla
%    anline1 = line([(xA(i)-sind(90-A_a5(i))*Wbla/2) (xA(i)+sind(90-A_a5(i))*Wbla/2)], ...
%    [(yA(i)+cosd(90-A_a5(i))*Wbla/2) (yA(i)-cosd(90-A_a5(i))*Wbla/2)], 'color', 'b', 'linewidth',3);
%end

%generalanline1 = line([(L1-sind(90-a5)*(Nbla*Wbla + (Nbla-1)*delta/sind(a5))/2) (L1+sind(90-a5)*(Nbla*Wbla + (Nbla-1)*delta/sind(a5))/2)], ...
%  [(0+cosd(90-a5)*(Nbla*Wbla + (Nbla-1)*delta/sind(a5))/2) (-cosd(90-a5)*(Nbla*Wbla + (Nbla-1)*delta/sind(a5))/2)], 'color', 'k', 'linewidth',0.7);


%for j = 1:Nbla
%lina1 = line([0 xA(j)], [0 yA(j)], 'color', 'r', 'linestyle', '--')
%leng = sqrt(xA(j)^2 + yA(j)^2);
%line1b = line([xA(j) (xA(j)+cosd(-2*a5+da5(j))*leng*1.2)],[yA(j) (yA(j)+sind(-2*a5+da5(j))*leng*1.2)], 'color', 'r', 'linestyle', '--')
%end

%plot(xR(1), yR(1), '*','color', 'm', 'linewidth', 1)
%plot(xR(2), yR(2), '*','color', 'm', 'linewidth', 1)

%set(gca, 'xlim', [Ana1(1)-80 Ana1(1)+80])
%set(gca, 'ylim', [-80 80])
%title('Analyzer @ E_f = 2.7 meV', 'fontname', 'times', 'fontsize', 15);
%set(gca, 'fontname', 'times', 'fontsize', 10);
%xlabel('Horizontal x-axis (along Q-channel) [mm]', 'fontname',  'times', 'fontsize', 12);
%ylabel('Vertical z-axis [mm]',  'fontname', 'times', 'fontsize', 12);

%daspect([1 1 1])

% 
% da5_1
% 
% a5
% 
% asind(sqrt((xA(1)-1000)^2+yA(1)^2)/RowR)
% yA(1)
% sqrt((xA(1)-1000)^2+yA(1)^2)/(2*pi*RowR)*360
% sqrt((xA(1)-1000)^2+yA(1)^2)/(4*pi*RowR)*360
% 
% (180/pi)*sqrt((xA(1)-1000)^2+yA(1)^2)*sind(a5)/1000
% 
% 
% da5_1
% da5_2
% da5_3
% da5_4
% da5_5
Wbla

widbla = round(Wbla*2)/2; %Calculating the width of each blade

Wana = 2*tand(Ana1(3)/2)*Ana1(1); %Calculating the width of the analyser

DetLength = Wana + 2*L2*tand(Ana1(3)/2) + 10; %Detectorlength using formula in note - extra 20 mm for margin - 1 cm on each side! 

Lbla = Wana + 2*tand(1)*L2*sind(a5);

area = widbla*Lbla;
AnaArea = area*Nbla/1000000;
%AnaArea = h*area;
dphi = []
Accept = 2*atand(0.5*Lbla/L1);
%save 'AnaData_2p7.mat' A_a5 xA yA width
	fid = fopen(['Bifrost_Analyzer_' num2str(Ana1(2)) 'meV_L1andL2is' num2str(L1) 'and' num2str(L2) 'mm_' num2str(Nbla) 'Blades.txt'], 'wt');
	h = length(xA);	
	fprintf(fid, 'Data in analyzer frame - central blade is at (0,0) \n');
    fprintf(fid, 'Y-axis is vertical - x-axis along the Q-channel center (analyzer axis) \n');
    fprintf(fid, ['Distance from sample: ' num2str(L1) ' mm \n']);
    fprintf(fid, ['Distance from detector: ' num2str(L2) ' mm \n']);
    fprintf(fid, ['Gap size: ' num2str(delta) ' mm \n']);
    fprintf(fid, ['Energy: ' num2str(En) ' meV \n']);
    fprintf(fid, ['Position of detector (in the global sample coordinate system - see figure) - X = ' num2str(Pd(1)) ' mm. Y = ' num2str(Pd(2)) ' mm.  \n']);
    fprintf(fid, ['2Theta coverage of analyzer: ' num2str(Ana1(3)) ' degrees \n']);
    fprintf(fid, ['Bragg angle: ' num2str(a5) ' degrees \n']);
    fprintf(fid, ['Length of each blade: ' num2str(Lbla, 4) ' mm \n']);
    fprintf(fid, ['Corresponding 2theta coverage: ' num2str(Accept, 4) ' degrees \n']);
    fprintf(fid, ['Proposed length of active detector tube: ' num2str(DetLength, 4) ' mm \n']);
    fprintf(fid, ['Width of each blade: ' num2str(widbla, 4) ' mm \n']);
    fprintf(fid, ['Total area of analyzer: ' num2str(h*area) ' mm \n']);
    fprintf(fid, ['Total area of all analyzers of this energy: ' num2str(9*h*area/1000000, 3) ' m^2 \n']);    
    fprintf(fid, ['Vertical angular coverage: ' num2str(Acc, 3) ' deg. \n']);    
    fprintf(fid, ['_____________________________________________________ \n']);
    fprintf(fid, ['Here follows a list of blade axis positions and analyzer take-off angle \n']); 
    fprintf(fid, ' Analyzer blade        x [mm]          y [mm]          angle [deg] \n');
	for n = 1:h
        n
	fprintf(fid, '%d \t \t %12.3f \t %12.3f \t %12.3f \n', [n xA(n)-Ana1(1) yA(n) A_a5(n)]);	
    if n > 1
        dphi(n-1) = abs(A_a5(n-1) - A_a5(n));
    end
    end
    fprintf(fid, ['_____________________________________________________ \n']);
    fprintf(fid, ['Average angle between adjacent blades: ' num2str(mean(dphi), 3) ' deg. \n']);    
    dphi
	fclose(fid)
    
