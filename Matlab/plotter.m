close all

fil=matfile('Datafile_Q','Writable', false);

Qqh=fil.QH;
Qqk=fil.QK;
Eee=fil.dE;
Iii=fil.I;

E_lim=[3 3.2];

E_ind1=find(Eee<E_lim(1)); E_ind2=reshape(find(Eee>E_lim(2)),length(find(Eee>E_lim(2))),1);
% E_in=[E_ind1 E_ind2];
Eee(E_ind1)=0;Eee(E_ind2)=0;
% Qqh(E_ind1)=NaN;Qqh(E_ind2)=NaN;
% Qqk(E_ind1)=NaN;Qqk(E_ind2)=NaN;
Iii(E_ind1)=0;Iii(E_ind2)=0;

F = scatteredInterpolant(Qqh,Qqk,Eee,Iii)

[Qa,Qb,Iii]=meshgrid(linspace(0.5,3.5,100),linspace(-2,2,100),linspace(1,5,100));
iI=F(Qa,Qb,Iii);

plot3(Qqh,Qqk,Eee,1,Iii)
mesh(Qa,Qb,iI)
% rotate(h,zdir,45,center)