% Hua-sheng XIE, huashengxie@gmail.com, ENN, 2021-05-16 17:15
% Solovev equilibrium for ST, FRC and mirror

function [psi,B,Br,Bz,Bphi,ns0,ts0,dBdr,dBdz,dBrdr,dBrdz,dBzdr,dBzdz,...
    dBphidr,dBphidz,dns0dr,dns0dz,iexit]=calpars_solovev(ra,za,nra,nza,...
    ffpsi,ffB,ffBr,ffBz,ffBphi,ffns0,ffts0,ffdBdr,ffdBdz,ffdBrdr,...
    ffdBrdz,ffdBzdr,ffdBzdz,ffdBphidr,ffdBphidz,ffdns0dr,ffdns0dz,S)

psi=ffpsi(ra,za);
B=ffB(ra,za);
Br=ffBr(ra,za);
Bz=ffBz(ra,za);
Bphi=ffBphi(ra,za);
dBdr=ffdBdr(ra,za);
dBdz=ffdBdz(ra,za);
dBrdr=ffdBrdr(ra,za);
dBrdz=ffdBrdz(ra,za);
dBzdr=ffdBzdr(ra,za);
dBzdz=ffdBzdz(ra,za);
dBphidr=ffdBphidr(ra,za);
dBphidz=ffdBphidz(ra,za);

ns0=zeros(S,nra,nza);
ts0=0.*ns0;
dns0dr=0.*ns0;
dns0dz=0.*ns0;
for s=1:S
    ns0(s,:,:)=ffns0(s,ra,za);
    ts0(s,:,:)=ffts0(s,ra,za);
    dns0dr(s,:,:)=ffdns0dr(s,ra,za);
    dns0dz(s,:,:)=ffdns0dz(s,ra,za);
end
iexit=0;

end