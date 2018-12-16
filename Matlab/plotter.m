close all

filA=matfile('Datafile_Q','Writable', false);
filB=matfile('Datafile_QB','Writable', false);

Qqh=filA.QH;
Qqk=filA.QK;
Eee=filA.dE;
Iii=filA.I;

Qqh=[Qqh filB.QH];
Qqk=[Qqk filB.QK];
Eee=[Eee filB.dE];
Iii=[Iii filB.I];

E_lim=[3.5 3.8];

E_ind1=find(Eee<E_lim(1)); E_ind2=reshape(find(Eee>E_lim(2)),length(find(Eee>E_lim(2))),1);
% E_in=[E_ind1 E_ind2];
Eee(E_ind1)=0;Eee(E_ind2)=0;
% Qqh(E_ind1)=NaN;Qqh(E_ind2)=NaN;
% Qqk(E_ind1)=NaN;Qqk(E_ind2)=NaN;
Iii(E_ind1)=0;Iii(E_ind2)=0;

sm=4;

Qqh=Qqh(1:sm:length(Qqh));Qqk=Qqk(1:sm:length(Qqk));Iii=Iii(1:sm:length(Iii));Eee=Eee(1:sm:length(Eee));

% F = scatteredInterpolant(Qqh,Qqk,Eee,Iii)
% 
% [Qa,Qb,Iii3]=meshgrid(linspace(0,3.5,100),linspace(-3,2,100),linspace(1,5,100));
% iI=F(Qa,Qb,Iii3);

figure(1)
x=0:0.1:3.5;y=-3:0.1:2;
z=zeros(length(y),length(x));
hold on
s=surf(x,y,z);
s.EdgeColor = 'none';
scatter3(Qqh,Qqk,Eee,1,Iii)
daspect([1 1 1])
xlabel('Q_H [Å^{-1}]')
ylabel('Q_K [Å^{-1}]')

% figure(2)
% plot3(Qa,Qb,iF)
% %mesh(Qa,Qb,)