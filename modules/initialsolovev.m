% Hua-sheng XIE, huashengxie@gmail.com, ENN-FTRC, 2021-05-16 10:20
% initialize Solovev equilibrium for ST, FRC and mirror
% clear;clc;

% iconfig=1;
% if(iconfig==1) % tokamak
%     R0=0.64; % m, major radius
%     B0=0.32; % T, B field at magnetic axis
%     q0=1.58; % safety factor at r=0
%     Rx=0.17; % r of X-point
%     E=1.5; % elongation
%     tau=0.8; %r 0.17,0.63,0.88; z~0.8
% elseif(iconfig==2) % frc
%     R0=0.5; % m, major radius
%     B0=0.5; % T, B field at (r=0,z=0)
%     E=2.0; % elongation
% elseif(iconfig==3) % mirror, to update
%     R0=0.5i; % m, major radius
%     B0=0.5; % T, B field at (r=0,z=0)
%     E=2.0; % elongation
% end
Lns=[0.9;0.9;0.9]*1.0; % density scale length
Lts=[0.8;0.8;0.8]; % temperature scale length

if(iconfig==1) % tokamak & ST
    kappa=2*E/sqrt(1-Rx^2/R0^2);
    psi0=B0*R0^2/(8*q0); % psi at O-point
    
    rmin=0.8*Rx; rmax=2.0*R0; % for plot grids
    zmax=R0*kappa*0.8; zmin=-zmax;
elseif(iconfig==2) % FRC
    Rx=0.0; tau=0.0;
    kappa=2*E/sqrt(1-Rx^2/R0^2);
    psi0=B0*R0^2/4; % normalized psi with B(r=0,z=0)=B0
    
    rmin=0.0*Rx; rmax=2.0*R0; % for plot grids
    zmax=R0*kappa; zmin=-zmax;
elseif(iconfig==3) % mirror
    Rx=0.0; tau=0.0;
%     kappa=2*E;
    kappa=2*E/sqrt(1-Rx^2/R0^2);
    psi0=B0*abs(R0^2)/4; % psi at O-point
    
    rmin=0.0*Rx; rmax=2.0*abs(R0); % for plot grids
    zmax=abs(R0)*kappa; zmin=-zmax;
end
Zx=abs(E*sqrt(tau*R0^2*log(Rx^2/R0^2+1e-10)+(2+tau)*(R0^2-Rx^2))); % Z of X-point

ffpsi=@(R,Z)psi0/R0^4*((R.^2-R0^2).^2+Z.^2/E^2.*(R.^2-Rx^2)-...
    tau*R0^2*(R.^2.*log((R+1e-10).^2/R0^2)-(R.^2-R0^2)-(R.^2-R0^2).^2/(2*R0^2)));
ffBr=@(R,Z)-2*psi0./(R+1e-10)./R0^4.*(Z./E^2.*(R.^2-Rx^2));
ffBz=@(R,Z)2*psi0/R0^4*(2*(R.^2-R0^2)+Z.^2/E^2-...
    tau*R0^2.*(log((R+1e-10).^2/R0^2)-(R.^2-R0^2)/R0^2));
if(iconfig==1)
    ffBphi=@(R,Z)B0*R0./(R+1e-10);
    ffdBphidr=@(R,Z)-B0*R0./(R+1e-10).^2;
elseif(iconfig==2 || iconfig==3)
    ffBphi=@(R,Z)0.*R;
    ffdBphidr=@(R,Z)0.*R;
end
ffdBphidz=@(R,Z)0.*R;
ffB=@(R,Z)sqrt(ffBr(R,Z).^2+ffBz(R,Z).^2+ffBphi(R,Z).^2);
ffdBrdr=@(R,Z)-2*psi0./R0^4.*(Z./E^2.*(1+Rx^2./R.^2));
ffdBrdz=@(R,Z)-2*psi0./(R+1e-10)./R0^4.*(1./E^2.*(R.^2-Rx^2));
ffdBzdr=@(R,Z)4*psi0/R0^4*R.*(2+tau-tau*R0^2./(R+1e-10).^2);
ffdBzdz=@(R,Z)4*psi0/R0^4*(Z./E^2);
ffdBdr=@(R,Z)(ffBr(R,Z).*ffdBrdr(R,Z)+ffBz(R,Z).*ffdBzdr(R,Z)+...
    ffBphi(R,Z).*ffdBphidr(R,Z))./ffB(R,Z);
ffdBdz=@(R,Z)(ffBr(R,Z).*ffdBrdz(R,Z)+ffBz(R,Z).*ffdBzdz(R,Z)+...
    ffBphi(R,Z).*ffdBphidz(R,Z))./ffB(R,Z);

Zx(isnan(Zx))=0.0;
psix=ffpsi(Rx,Zx); % psi at separatrix X-point

qe=1.60217662e-19; % electron charge, coulombs
mp=1.6726219e-27; % proton mass, kg
me=9.1094e-31; % electron mass, kgs
qs=qs0*qe;
ms=ms0*me;
S=length(qs);

% set density and temperature profiles
ffns0=@(s,R,Z)ns00(s)*exp(-ffpsi(R,Z)./(psix*Lns(s)^2));
ffts0=@(s,R,Z)ts00(s)*exp(-ffpsi(R,Z)./(psix*Lts(s)^2));
ffdns0dpsi=@(s,R,Z)-ffns0(s,R,Z)./(psix*Lns(s)^2);
ffdns0dr=@(s,R,Z)R.*ffBz(R,Z).*ffdns0dpsi(s,R,Z);
ffdns0dz=@(s,R,Z)-R.*ffBr(R,Z).*ffdns0dpsi(s,R,Z);

%%
dr=0.02*(rmax-rmin);
dz=0.02*(zmax-zmin);
rg=rmin:dr:rmax;
zg=zmin:dz:zmax;
[rr,zz]=ndgrid(rg,zg);
[nr,nz]=size(rr);

[fpsi,fB,fBr,fBz,fBphi,fns0,fts0,fdBdr,fdBdz,fdBrdr,fdBrdz,fdBzdr,fdBzdz,...
    fdBphidr,fdBphidz,fdns0dr,fdns0dz,iexit]=calpars_solovev(rr,zz,nr,nz,...
    ffpsi,ffB,ffBr,ffBz,ffBphi,ffns0,ffts0,ffdBdr,ffdBdz,ffdBrdr,...
    ffdBrdz,ffdBzdr,ffdBzdz,ffdBphidr,ffdBphidz,ffdns0dr,ffdns0dz,S);

n0=max(max(fns0(1,:,:)));
%%
if(1==0)
    close all;
    if(iconfig==1)
        xx=rr;yy=zz;
    elseif(iconfig==2 ||iconfig==3)
        xx=zz;yy=rr;
    end
    
    [ccpsi,hpsi]=contour(xx,yy,fpsi,[psix,psix],'visible','off');
    if(isempty(ccpsi))
        ccpsi=zeros(2,2);
    end
    
    subplot(231);
    contour(xx,yy,fpsi,10.^(-2:0.1:0.5)*psix);colorbar;hold on;
    plot(ccpsi(1,2:ccpsi(2,1)),ccpsi(2,2:ccpsi(2,1)),'r--');
    axis equal; title('\psi');
    subplot(232);
    contour(xx,yy,fB,100);colorbar;hold on;
    plot(ccpsi(1,2:ccpsi(2,1)),ccpsi(2,2:ccpsi(2,1)),'r--');
    axis equal; title('B');
    subplot(233);
    contour(xx,yy,squeeze(fns0(1,:,:)),100);colorbar;hold on;
    plot(ccpsi(1,2:ccpsi(2,1)),ccpsi(2,2:ccpsi(2,1)),'r--');
    axis equal; title('n_{s1}');
    subplot(234);
    contour(xx,yy,squeeze(fdns0dz(1,:,:)),100);colorbar;hold on;
    plot(ccpsi(1,2:ccpsi(2,1)),ccpsi(2,2:ccpsi(2,1)),'r--');
    axis equal; title('dn_{s1}/dz');
    subplot(235);
    plot(rg,squeeze(fns0(1,:,floor(nz/2))));hold on;
    xlabel('r');ylabel('n_{s0}');
    subplot(236);
    plot(rg,squeeze(fts0(1,:,floor(nz/2))));hold on;
    xlabel('r');ylabel('t_{s0}');
%     figure;plot(rg,fB(:,floor(nz/2)),rg,fBphi(:,floor(nz/2)));
end