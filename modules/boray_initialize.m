% 2021-05-05 13:25 Huasheng XIE, huashengxie@gmail.com, ENN
% initialize.m

if(jray==1)
    % default SI unit
    c2=(2.99792458E8)^2; % speed of ligth c^2
    epsilon0=8.854187817E-12;
    c=sqrt(c2);
    
    if(numeq==1)
        % initialize plasma density ns0(r,z) and magnetic field B(r,z)
        load(['../input/',eqfile]);
        
        % initial interp2d coefficients
        [fcB,fcBr,fcBz,fcBphi,fcpsi,fcns0,fcts0,fcdBdr,fcdBdz,fcdBrdr,...
            fcdBrdz,fcdBzdr,fcdBzdz,fcdBphidr,fcdBphidz,fcdns0dr,fcdns0dz]=...
            fcinterp(rg,zg,fB,fBr,fBz,fBphi,fpsi,fns0,fts0,fdBdr,fdBdz,...
            fdBrdr,fdBrdz,fdBzdr,fdBzdz,fdBphidr,fdBphidz,fdns0dr,fdns0dz,S);
        
        [nr,nz]=size(fB);
    elseif(numeq==0) % 21-05-16 17:09 for Solovev analytical equilibrium
        run initialsolovev;
    end
    
    % wave frequency
    w=2*pi*f; % rad/s
    w2=w^2;
    
    if(~exist(savepath,'dir')) % in case savepath not exist
        mkdir(savepath);
    end
    
end

%
if(numeq==1)
    [jr,jz,psi,B,Br,Bz,Bphi,ns0,ts0,dBdr,dBdz,dBrdr,dBrdz,dBzdr,dBzdz,...
      dBphidr,dBphidz,dns0dr,dns0dz,iexit]=calpars(r,z,rg,zg,dr,dz,...
      fcpsi,fcB,fcBr,fcBz,fcBphi,fcns0,fcts0,fcdBdr,fcdBdz,fcdBrdr,...
      fcdBrdz,fcdBzdr,fcdBzdz,fcdBphidr,fcdBphidz,fcdns0dr,fcdns0dz,S,nr,nz);
elseif(numeq==0)
    [psi,B,Br,Bz,Bphi,ns0,ts0,dBdr,dBdz,dBrdr,dBrdz,dBzdr,dBzdz,...
    dBphidr,dBphidz,dns0dr,dns0dz,iexit]=calpars_solovev(r,z,1,1,...
    ffpsi,ffB,ffBr,ffBz,ffBphi,ffns0,ffts0,ffdBdr,ffdBdz,ffdBrdr,...
    ffdBrdz,ffdBzdr,ffdBzdz,ffdBphidr,ffdBphidz,ffdns0dr,ffdns0dz,S);
end
if(iexit)
    disp('The ray is out of range! it=0');
    return;
end

%
wcs=qs.*B./ms; % gyro frequency
wps2=ns0.*qs.^2./(epsilon0*ms);

% eps1=S, eps2=D, eps3=P
eps1=1-sum(wps2./(w2-wcs.^2));
eps2=sum((wcs./w).*wps2./(w2-wcs.^2));
eps3=1-sum(wps2./w2);
fkpar2=@(kr)((kr*Br+kz*Bz+nphi/r*Bphi)/B).^2;
fkper2=@(kr) (kr.^2+kz^2+nphi^2/r^2)-fkpar2(kr);

% define the function to solve the initial kr
fDkr=@(kr)eps1*(fkper2(kr)*c2/w2).^2-...
    ((eps1+eps3)*(eps1-fkpar2(kr)*c2/w2)-eps2^2).*fkper2(kr)*c2/w2+...
    eps3*((eps1-fkpar2(kr)*c2/w2).^2-eps2^2);

