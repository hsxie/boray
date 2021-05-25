% Hua-sheng XIE, huashengxie@gmail.com, ENN, 2021-05-04 21:30 
% Bi-linear, (r,z) uniform grid
% Calculate the interp2d coefficients fc(1:4,:)

function [fcB,fcBr,fcBz,fcBphi,fcpsi,fcns0,fcts0,fcdBdr,fcdBdz,fcdBrdr,...
    fcdBrdz,fcdBzdr,fcdBzdz,fcdBphidr,fcdBphidz,fcdns0dr,fcdns0dz]=...
    fcinterp(rg,zg,fB,fBr,fBz,fBphi,fpsi,fns0,fts0,fdBdr,fdBdz,...
    fdBrdr,fdBrdz,fdBzdr,fdBzdz,fdBphidr,fdBphidz,fdns0dr,fdns0dz,S)
    
% load(fname); % including the eqdata
% S=length(qs);
[nr,nz]=size(fB);

fcB=cinterp2d(rg,zg,fB);
fcBr=cinterp2d(rg,zg,fBr);
fcBz=cinterp2d(rg,zg,fBz);
fcBphi=cinterp2d(rg,zg,fBphi);
fcpsi=cinterp2d(rg,zg,fpsi);

fcdBdr=cinterp2d(rg,zg,fdBdr);
fcdBdz=cinterp2d(rg,zg,fdBdz);
fcdBrdr=cinterp2d(rg,zg,fdBrdr);
fcdBrdz=cinterp2d(rg,zg,fdBrdz);
fcdBzdr=cinterp2d(rg,zg,fdBzdr);
fcdBzdz=cinterp2d(rg,zg,fdBzdz);
fcdBphidr=cinterp2d(rg,zg,fdBphidr);
fcdBphidz=cinterp2d(rg,zg,fdBphidz);

fcns0=zeros(S,4,nr,nz);
fcts0=0.*fcns0;
fcdns0dr=0.*fcns0;
fcdns0dz=0.*fcns0;
for s=1:S
    fcns0(s,:,:,:)=cinterp2d(rg,zg,squeeze(fns0(s,:,:)));
    fcts0(s,:,:,:)=cinterp2d(rg,zg,squeeze(fts0(s,:,:)));
    fcdns0dr(s,:,:,:)=cinterp2d(rg,zg,squeeze(fdns0dr(s,:,:)));
    fcdns0dz(s,:,:,:)=cinterp2d(rg,zg,squeeze(fdns0dz(s,:,:)));
end

end