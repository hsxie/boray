% Hua-sheng XIE, huashengxie@gmail.com, ENN, 2021-05-04 21:39
% Bi-linear, (r,z) uniform grid.
% Use the interp2d coefficients fc(1:4,:)

function [jr,jz,psi,B,Br,Bz,Bphi,ns0,ts0,dBdr,dBdz,dBrdr,dBrdz,dBzdr,dBzdz,...
    dBphidr,dBphidz,dns0dr,dns0dz,iexit]=calpars(ra,za,...
    rg,zg,dr,dz,fcpsi,fcB,fcBr,fcBz,fcBphi,fcns0,fcts0,fcdBdr,fcdBdz,fcdBrdr,...
    fcdBrdz,fcdBzdr,fcdBzdz,fcdBphidr,fcdBphidz,fcdns0dr,fcdns0dz,S,nr,nz)

% S=size(fcns0,1);
jr=floor((ra-rg(1))/dr)+1; jz=floor((za-zg(1))/dz)+1;
if ((jr>=nr) || (jr<1) || (jz>=nz) || (jz<1))
    
    iexit=1;
%     hr=0; hz=0;
    psi=0; B=0; Br=0; Bz=0; Bphi=0; ns0=0; ts0=0;
    dBdr=0; dBdz=0; dBrdr=0; dBrdz=0; dBzdr=0; dBzdz=0;
    dBphidr=0; dBphidz=0; dns0dr=0; dns0dz=0;
    
else
    hr=ra-rg(jr); hz=za-zg(jz);

    psi=fcpsi(1,jr,jz)+fcpsi(2,jr,jz)*hr+fcpsi(3,jr,jz)*hz+fcpsi(4,jr,jz)*hz*hr;
    B=fcB(1,jr,jz)+fcB(2,jr,jz)*hr+fcB(3,jr,jz)*hz+fcB(4,jr,jz)*hz*hr;
    dBdr=fcdBdr(1,jr,jz)+fcdBdr(2,jr,jz)*hr+fcdBdr(3,jr,jz)*hz+...
        fcdBdr(4,jr,jz)*hz*hr;
    dBdz=fcdBdz(1,jr,jz)+fcdBdz(2,jr,jz)*hr+fcdBdz(3,jr,jz)*hz+...
        fcdBdz(4,jr,jz)*hz*hr;
    
    Br=fcBr(1,jr,jz)+fcBr(2,jr,jz)*hr+fcBr(3,jr,jz)*hz+fcBr(4,jr,jz)*hz*hr;
    dBrdr=fcdBrdr(1,jr,jz)+fcdBrdr(2,jr,jz)*hr+fcdBrdr(3,jr,jz)*hz+...
        fcdBrdr(4,jr,jz)*hz*hr;
    dBrdz=fcdBrdz(1,jr,jz)+fcdBrdz(2,jr,jz)*hr+fcdBrdz(3,jr,jz)*hz+...
        fcdBrdz(4,jr,jz)*hz*hr;
    
    Bz=fcBz(1,jr,jz)+fcBz(2,jr,jz)*hr+fcBz(3,jr,jz)*hz+fcBz(4,jr,jz)*hz*hr;
    dBzdr=fcdBzdr(1,jr,jz)+fcdBzdr(2,jr,jz)*hr+fcdBzdr(3,jr,jz)*hz+...
        fcdBzdr(4,jr,jz)*hz*hr;
    dBzdz=fcdBzdz(1,jr,jz)+fcdBzdz(2,jr,jz)*hr+fcdBzdz(3,jr,jz)*hz+...
        fcdBzdz(4,jr,jz)*hz*hr;
    
    Bphi=fcBphi(1,jr,jz)+fcBphi(2,jr,jz)*hr+fcBphi(3,jr,jz)*hz+...
        fcBphi(4,jr,jz)*hz*hr;
    dBphidr=fcdBphidr(1,jr,jz)+fcdBphidr(2,jr,jz)*hr+fcdBphidr(3,jr,jz)*hz+...
        fcdBphidr(4,jr,jz)*hz*hr;
    dBphidz=fcdBphidz(1,jr,jz)+fcdBphidz(2,jr,jz)*hr+fcdBphidz(3,jr,jz)*hz+...
        fcdBphidz(4,jr,jz)*hz*hr;
    
    ns0=zeros(S,1);
    dns0dr=0.*ns0; dns0dz=0.*ns0; ts0=0.*ns0;
    for s=1:S
        ns0(s)=fcns0(s,1,jr,jz)+fcns0(s,2,jr,jz)*hr+fcns0(s,3,jr,jz)*hz+...
            fcns0(s,4,jr,jz)*hz*hr;
        dns0dr(s)=fcdns0dr(s,1,jr,jz)+fcdns0dr(s,2,jr,jz)*hr+...
            fcdns0dr(s,3,jr,jz)*hz+fcdns0dr(s,4,jr,jz)*hz*hr;
        dns0dz(s)=fcdns0dz(s,1,jr,jz)+fcdns0dz(s,2,jr,jz)*hr+...
            fcdns0dz(s,3,jr,jz)*hz+fcdns0dz(s,4,jr,jz)*hz*hr;
        
        ts0(s)=fcts0(s,1,jr,jz)+fcts0(s,2,jr,jz)*hr+fcts0(s,3,jr,jz)*hz+...
            fcts0(s,4,jr,jz)*hz*hr;
    end
    iexit=0;
end
end