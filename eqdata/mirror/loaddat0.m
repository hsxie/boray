% Hua-sheng XIE, huashengxie@gmail.com, ENN, 2021-04-28 16:50

close all; clear; clc;

icase=10;

load('B_rz.mat');

ir=2;iz=2;
Br=squeeze(Dat(1,1:ir:end,1:iz:end));
Bphi=squeeze(Dat(2,1:ir:end,1:iz:end));
Bz=squeeze(Dat(3,1:ir:end,1:iz:end));

[nr,nz]=size(Br);
rmin=0; dr=0.5e-3*ir;
zmin=0; dz=4e-3*iz;

rg=rmin+dr*(0:1:(nr-1));
zg=zmin+dz*(0:1:(nz-1));
zg=zg-mean(zg);
[rr,zz]=ndgrid(rg,zg);

Br(abs(Br)>3)=0;
Bz(abs(Bz)>3)=0;

if(icase==10)
    n0=1e18;
%     n0=3e18;
elseif(icase==11)
    n0=1e20;
    ctmp=10;
    Br=ctmp*Br;
    Bz=ctmp*Bz;
    Bphi=ctmp*Bphi;
end

sgm=4.5e-2;
% nrz=n0*exp(-rr.^2/sgm^2);
nrz=n0*exp(-rr.^2/sgm^2/2);

psi=0.*rr;
for jr=2:nr
    for jz=1:nz
        psi(jr,jz)=psi(jr-1,jz)+dr*0.5*(Bz(jr,jz)+...
            Bz(jr-1,jz))*0.5*(rr(jr,jz)+rr(jr-1,jz));
    end
end

%
% close all;
% subplot(221);
% contour(zz,rr,Br,100);
% subplot(222);
% contour(zz,rr,Bz,100);
% 
% subplot(223);
% contour(zz,rr,nrz,100);
% % plot(r,Bz(:,floor(nz/2)));
% 
% subplot(224);
% plot(z,Bz(1,:));

%


psi(isnan(psi))=0.0; Br(isnan(Br))=0.0; Bz(isnan(Bz))=0.0;
psi(isinf(psi))=0.0; Br(isinf(Br))=0.0; Bz(isinf(Bz))=0.0;

fpsi=psi;
fBr=Br;
fBz=Bz;
fBphi=0.*rr;
% if(icase==10)
fni=nrz;
% end
fne=fni;
fti=460+0.*rr;
fte=460+0.*rr;

fB=sqrt(fBr.^2+fBz.^2+fBphi.^2);


%%
qe=1.60217662e-19; % electron charge, coulombs
mp=1.6726219e-27; % proton mass, kg
me=9.1094e-31; % electron mass, kgs
qs=[1;-1]*qe;
ms=[mp;me];
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

B0=fB(1,floor(nz/2)); R0=sgm;
S=length(qs);
save('zjs_eqdata.mat','rg','zg','dr','dz','rr','zz','fB','fBr','fBz',...
    'fBphi','fns0','fdBdr','fdBdz',...
    'fdBrdr','fdBrdz','fdBzdr','fdBzdz','fdBphidr','fdBphidz',...
    'fdns0dr','fdns0dz','fts0','qs','ms','n0','R0','B0','fpsi','S');

% save('genray_eqdata.mat','rg','zg','dr','dz','rr','zz',...
%     'fB','fBr','fBz','fBphi','fns0','fts0','fdBdr','fdBdz',...
%     'fdBrdr','fdBrdz','fdBzdr','fdBzdz','fdBphidr','fdBphidz',...
%     'fdns0dr','fdns0dz','qs','ms','R0','Z0','B0','n0','fpsi','S');

%%
% figure;subplot(121);surf(rr,zz,fdBphidr);subplot(122);surf(rr,zz,fdBphidz);

%%
close all;
h=figure('unit','normalized','Position',[0.01 0.05 0.6 0.7],...
  'DefaultAxesFontSize',14);
subplot(421);
contour(zz,rr,fpsi,50);
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

print(gcf,'-dpng',['eqdata_icase=',num2str(icase),',B0=',num2str(B0,3),...
    'T,ni0=',num2str(n0,3),'.png']);

