% Hua-sheng XIE, huashengxie@gmail.com, ENN, 2021-04-20 10:28
% Cold Plasma Dispersion Reletion
function [ww,detD]=colddr(qs,ms,fns0,fB0,fkx,fkz,fw)
% size:
%   =S. qs, ms
%   =npoint*S. ns0
%   =npoint. B0, kx, kz, w

% qs=[1,-1]*1.60217662e-19;
% ms=[1836,1]*9.1094e-31;
% ns0=[1e18,1e18];
% B0=1.0;
% kx=0.0;
% kz=1.0;
% w=1e9;

% 1. default SI unit
c2=(2.99792458E8)^2; % speed of ligth c^2
epsilon0=8.854187817E-12;
% mu0=1/(c2*epsilon0);
% kB=1.38064852e-23;
% qe=1.60217662e-19; % electron charge, coulombs
% mp=1.6726219e-27; % proton mass, kg
% me=9.1094e-31; % electron mass, kg

% S=length(qs);
% [S,npoint]=size(fns0);
[npoint,S]=size(fns0);

ww=zeros(npoint,1);
detD=zeros(npoint,1);
for jp=1:npoint
    
    ns0=squeeze(fns0(jp,:));
    B0=fB0(jp);
    kx=fkx(jp);
    kz=fkz(jp);
    w=fw(jp);
    
    wps=sqrt(ns0.*qs.^2./ms/epsilon0); % plasma frequency
    wcs=B0*qs./ms; % cyclotron frequency
    wps2=wps.^2;
    
%     k=sqrt(kz^2+kx^2);
    
    % % Part 1: use matrix method to solve omega with given kx,kz
    NN=3*S+6;
    SJ=3*S;
    % eigen matrix
    M=sparse(NN,NN);
    for s=1:S
        ind=(s-1)*3;
        % dv ~ v
        M=M+sparse(ind+2,ind+1,-1i*wcs(s),NN,NN)+...
            sparse(ind+1,ind+2,1i*wcs(s),NN,NN);
        % dv ~ E
        M=M+sparse(ind+1,SJ+1,1i*qs(s)/ms(s),NN,NN)+...
            sparse(ind+2,SJ+2,1i*qs(s)/ms(s),NN,NN)+...
            sparse(ind+3,SJ+3,1i*qs(s)/ms(s),NN,NN);
        % dE ~ v
        M=M+sparse(SJ+1,ind+1,-1i*qs(s)*ns0(s)/epsilon0,NN,NN)+...
            sparse(SJ+2,ind+2,-1i*qs(s)*ns0(s)/epsilon0,NN,NN)+...
            sparse(SJ+3,ind+3,-1i*qs(s)*ns0(s)/epsilon0,NN,NN);
    end
    % E(B)
    M=M+sparse(SJ+1,SJ+5,c2*kz,NN,NN)+...
        sparse(SJ+2,SJ+4,-c2*kz,NN,NN)+...
        sparse(SJ+2,SJ+6,c2*kx,NN,NN)+...
        sparse(SJ+3,SJ+5,-c2*kx,NN,NN);
    % B(E)
    M=M+sparse(SJ+4,SJ+2,-kz,NN,NN)+...
        sparse(SJ+5,SJ+1,kz,NN,NN)+...
        sparse(SJ+5,SJ+3,-kx,NN,NN)+...
        sparse(SJ+6,SJ+2,kx,NN,NN);  
    % ww=eigs(M,NN);
    ww(jp)=eigs(M,1,w);
    
    % % Part 2: calculate detD
    w2=w*w;
    % eps1=S, eps2=D, eps3=P
    eps1=1-sum(wps2./(w2-wcs.^2));
    eps2=sum((wcs./w).*wps2./(w2-wcs.^2));
    eps3=1-sum(wps2./w2);
    kpar2=kz*kz;
    kper2=kx*kx;
    detD(jp)=eps1*(kper2*c2/w2).^2-...
        ((eps1+eps3)*(eps1-kpar2*c2/w2)-eps2^2).*kper2*c2/w2+...
        eps3*((eps1-kpar2*c2/w2).^2-eps2^2);
end

end




