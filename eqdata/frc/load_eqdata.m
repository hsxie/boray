% 2021-04-26 13:10 Huasheng XIE, huashengxie@gmail.com, ENN
% load gseq FRC equilibrium data. Ma2021 NF paper#2 Fig.4
% Note that do not oppsite r and z index
% 21-05-15 16:37 update for new boray

close all;
clear; clc;

load('RR_C2_eqdata.mat');
% load('C2U_gseq.mat'); % 21-05-20 22:27
dr=(R(2,1)-R(1,1))*1; dz=(Z(1,2)-Z(1,1))*1;
% tmp=1/2;
%%
rg=R(1,1):dr:R(end,end);
zg=Z(1,1):dz:Z(end,end);

[rr,zz]=ndgrid(rg,zg);
[nr,nz]=size(rr);

tmp=abs(Bz(1,floor(nz/2)))/0.05;


psi=psi/tmp; Br=Br/tmp; Bz=Bz/tmp; P=P/tmp;
psi(isnan(psi))=0.0; Br(isnan(Br))=0.0; Bz(isnan(Bz))=0.0; P(isnan(P))=0.0;

fpsi=interp2(R.',Z.',psi.',rr,zz);
%%
fBr=interp2(R.',Z.',Br.',rr,zz);
fBz=interp2(R.',Z.',Bz.',rr,zz);
fBphi=0.*rr;
fP=interp2(R.',Z.',P.',rr,zz)/max(max(P)); fP=fP./max(max(fP));
fni=sqrt(fP)*0.24e20;
fne=fni;
% fti=400+0.*rr;
% fte=100+0.*rr;
fti=sqrt(fP)*800; % 21-05-17 23:08 
fte=sqrt(fP)*150;

fB=sqrt(fBr.^2+fBz.^2+fBphi.^2);

id0=find(fpsi==min(min(fpsi)));
B0=fB(id0);
Be=fB(floor(nr/2),nz);
R0=rr(id0);
Z0=zz(id0);
n0=fni(id0);

%%
qe=1.60217662e-19; % electron charge, coulombs
mp=1.6726219e-27; % proton mass, kg
me=9.1094e-31; % electron mass, kgs
qs=[1;-1]*qe;
ms=[2*mp;me];
fns0(1,:,:)=fni;
fns0(2,:,:)=fne;
fts0(1,:,:)=fti;
fts0(2,:,:)=fte;

if(1==1)
fdBrdr=0.*rr; fdBrdr(2:(nr-1),:)=(fBr(3:nr,:)-fBr(1:(nr-2),:))./(2*dr);
fdBrdz=0.*rr; fdBrdz(:,2:(nz-1))=(fBr(:,3:nz)-fBr(:,1:(nz-2)))./(2*dz);
fdBzdr=0.*rr; fdBzdr(2:(nr-1),:)=(fBz(3:nr,:)-fBz(1:(nr-2),:))./(2*dr);
fdBzdz=0.*rr; fdBzdz(:,2:(nz-1))=(fBz(:,3:nz)-fBz(:,1:(nz-2)))./(2*dz);
fdBphidr=0.*rr; fdBphidr(2:(nr-1),:)=(fBphi(3:nr,:)-fBphi(1:(nr-2),:))./(2*dr);
fdBphidz=0.*rr; fdBphidz(:,2:(nz-1))=(fBphi(:,3:nz)-fBphi(:,1:(nz-2)))./(2*dz);
fdBdr=0.*rr; fdBdr(2:(nr-1),:)=(fB(3:nr,:)-fB(1:(nr-2),:))./(2*dr);
fdBdz=0.*rr; fdBdz(:,2:(nz-1))=(fB(:,3:nz)-fB(:,1:(nz-2)))./(2*dz);

fdnedr=0.*rr; fdnedr(2:(nr-1),:)=(fne(3:nr,:)-fne(1:(nr-2),:))./(2*dr);
fdnedz=0.*rr; fdnedz(:,2:(nz-1))=(fne(:,3:nz)-fne(:,1:(nz-2)))./(2*dz);
fdnidr=0.*rr; fdnidr(2:(nr-1),:)=(fni(3:nr,:)-fni(1:(nr-2),:))./(2*dr);
fdnidz=0.*rr; fdnidz(:,2:(nz-1))=(fni(:,3:nz)-fni(:,1:(nz-2)))./(2*dz);

else  % 21-04-07 13:27 fixed r,z index, which seems should be opposite
fdBrdr=0.*rr; fdBrdr(:,2:(nr-1))=(fBr(:,3:nr)-fBr(:,1:(nr-2)))./(2*dr);
fdBrdz=0.*rr; fdBrdz(2:(nz-1),:)=(fBr(3:nz,:)-fBr(1:(nz-2),:))./(2*dz);
fdBzdr=0.*rr; fdBzdr(:,2:(nr-1))=(fBz(:,3:nr)-fBz(:,1:(nr-2)))./(2*dr);
fdBzdz=0.*rr; fdBzdz(2:(nz-1),:)=(fBz(3:nz,:)-fBz(1:(nz-2),:))./(2*dz);
fdBphidr=0.*rr; fdBphidr(:,2:(nr-1))=(fBphi(:,3:nr)-fBphi(:,1:(nr-2)))./(2*dr);
fdBphidz=0.*rr; fdBphidz(2:(nz-1),:)=(fBphi(3:nz,:)-fBphi(1:(nz-2),:))./(2*dz);
fdBdr=0.*rr; fdBdr(:,2:(nr-1))=(fB(:,3:nr)-fB(:,1:(nr-2)))./(2*dr);
fdBdz=0.*rr; fdBdz(2:(nz-1),:)=(fB(3:nz,:)-fB(1:(nz-2),:))./(2*dz);

fdnedr=0.*rr; fdnedr(:,2:(nr-1))=(fne(:,3:nr)-fne(:,1:(nr-2)))./(2*dr);
fdnedz=0.*rr; fdnedz(2:(nz-1),:)=(fne(3:nz,:)-fne(1:(nz-2),:))./(2*dz);
fdnidr=0.*rr; fdnidr(:,2:(nr-1))=(fni(:,3:nr)-fni(:,1:(nr-2)))./(2*dr);
fdnidz=0.*rr; fdnidz(2:(nz-1),:)=(fni(3:nz,:)-fni(1:(nz-2),:))./(2*dz);
end
fdns0dr(1,:,:)=fdnidr;
fdns0dz(1,:,:)=fdnidz;
fdns0dr(2,:,:)=fdnedr;
fdns0dz(2,:,:)=fdnedz;

S=length(qs);psilim=0.0;
% save('c2u_eqdata.mat','rg','zg','dr','dz','rr','zz','fB','fBr','fBz',...
%     'fBphi','fns0','fdBdr','fdBdz',...
%     'fdBrdr','fdBrdz','fdBzdr','fdBzdz','fdBphidr','fdBphidz',...
%     'fdns0dr','fdns0dz','fts0','qs','ms','R0','Z0','B0','n0','fpsi','S');

save('c2u_eqdata.mat','rg','zg','dr','dz','rr','zz',...
    'fB','fBr','fBz','fBphi','fns0','fts0','fdBdr','fdBdz',...
    'fdBrdr','fdBrdz','fdBzdr','fdBzdz','fdBphidr','fdBphidz',...
    'fdns0dr','fdns0dz','qs','ms','R0','Z0','B0','n0','fpsi','S','psilim');
%%
% figure;subplot(121);surf(rr,zz,fdBphidr);subplot(122);surf(rr,zz,fdBphidz);

%%
close all;
h=figure('unit','normalized','Position',[0.01 0.05 0.6 0.7],...
  'DefaultAxesFontSize',14);
subplot(421);
fpsiout=fpsi; fpsiout(fpsi<0)=NaN;
fpsiin=fpsi; fpsiin(fpsi>0)=NaN;
contour(zz,rr,fpsiout,50); hold on;
contour(zz,rr,fpsiin,10); hold on;
contour(zz,rr,fpsi,[0,0],'r--');
xlabel('Z');ylabel('R');title('\psi');
axis tight; colorbar;

subplot(422);
contour(zz,rr,fB,100);
xlabel('Z');ylabel('R');title('B');
axis tight; colorbar;

subplot(423);
contour(zz,rr,fBr,100);
xlabel('R');ylabel('Z');title('B_r');
axis tight; colorbar;

subplot(424);
contour(zz,rr,fBz,100);
xlabel('R');ylabel('Z');title('B_z');
axis tight; colorbar;

subplot(425);
contour(zz,rr,fni,100);
xlabel('R');ylabel('Z');title('n_{i}');
axis tight; colorbar;

subplot(426);
contour(zz,rr,fti,100);
xlabel('R');ylabel('Z');title('T_{i}');
axis tight; colorbar;

subplot(427);
contour(zz,rr,fdnedr,100);
xlabel('R');ylabel('Z');title('dn_e/dr');
axis tight; colorbar;

subplot(428);
contour(zz,rr,fdnedz,100);
xlabel('R');ylabel('Z');title('dn_e/dz');
axis tight; colorbar;


set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
  'PaperSize',[screenposition(3:4)]);

print(gcf,'-dpng',['eqdata_,B0=',num2str(B0,3),...
    'T,ni0=',num2str(n0,3),'.png']);

