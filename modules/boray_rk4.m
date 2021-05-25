% 2021-05-05 13:34 Huasheng XIE, huashengxie@gmail.com, ENN
% rk4.m
% 21-05-17 10:42 to correction the Hamilitonian along ray

yy=zeros(nt,18+2*S)+NaN;
for it=1:nt
    %         it
    yy(it,1:6)=[r,phi,z,kr,nphi,kz];
    
    %     Nphi=nphi*c/w/r
    irk4=1;
    if(irk4==1)
        % % RK-4, 1st step
        if(numeq==1)
            [jr,jz,psi,B,Br,Bz,Bphi,ns0,ts0,dBdr,dBdz,dBrdr,dBrdz,dBzdr,dBzdz,...
                dBphidr,dBphidz,dns0dr,dns0dz,iexit]=calpars(r,z,...
                rg,zg,dr,dz,fcpsi,fcB,fcBr,fcBz,fcBphi,fcns0,fcts0,fcdBdr,fcdBdz,fcdBrdr,...
                fcdBrdz,fcdBzdr,fcdBzdz,fcdBphidr,fcdBphidz,fcdns0dr,fcdns0dz,S,nr,nz);
        elseif(numeq==0)
            [psi,B,Br,Bz,Bphi,ns0,ts0,dBdr,dBdz,dBrdr,dBrdz,dBzdr,dBzdz,...
                dBphidr,dBphidz,dns0dr,dns0dz,iexit]=calpars_solovev(r,z,1,1,...
                ffpsi,ffB,ffBr,ffBz,ffBphi,ffns0,ffts0,ffdBdr,ffdBdz,ffdBrdr,...
                ffdBrdz,ffdBzdr,ffdBzdz,ffdBphidr,ffdBphidz,ffdns0dr,ffdns0dz,S);
        end
        if(iexit)
            disp(['The ray is out of range! it=',num2str(it)]); return;
        end
        
        [kpar1,kper21,dFdr1,dFdphi1,dFdz1,dFdkr1,dFdnphi1,dFdkz1,dFdw1,D1]=dydt(...
            r,phi,z,kr,nphi,kz,...
            B,Br,Bz,Bphi,ns0,dBdr,dBdz,dBrdr,dBrdz,dBzdr,dBzdz,...
            dBphidr,dBphidz,dns0dr,dns0dz,qs,ms,w,c2);
        
        % % RK-4, 2nd step
        rtmp=yy(it,1)-dFdkr1/dFdw1*dt;
        phitmp=yy(it,2)-dFdnphi1/dFdw1*dt;
        ztmp=yy(it,3)-dFdkz1/dFdw1*dt;
        krtmp=yy(it,4)+dFdr1/dFdw1*dt;
        nphitmp=yy(it,5)+dFdphi1/dFdw1*dt;
        kztmp=yy(it,6)+dFdz1/dFdw1*dt;
        
        if(numeq==1)
            [jr,jz,psi,B,Br,Bz,Bphi,ns0,ts0,dBdr,dBdz,dBrdr,dBrdz,dBzdr,dBzdz,...
                dBphidr,dBphidz,dns0dr,dns0dz,iexit]=calpars(rtmp,ztmp,...
                rg,zg,dr,dz,fcpsi,fcB,fcBr,fcBz,fcBphi,fcns0,fcts0,fcdBdr,fcdBdz,fcdBrdr,...
                fcdBrdz,fcdBzdr,fcdBzdz,fcdBphidr,fcdBphidz,fcdns0dr,fcdns0dz,S,nr,nz);
        elseif(numeq==0)
            [psi,B,Br,Bz,Bphi,ns0,ts0,dBdr,dBdz,dBrdr,dBrdz,dBzdr,dBzdz,...
                dBphidr,dBphidz,dns0dr,dns0dz,iexit]=calpars_solovev(rtmp,ztmp,1,1,...
                ffpsi,ffB,ffBr,ffBz,ffBphi,ffns0,ffts0,ffdBdr,ffdBdz,ffdBrdr,...
                ffdBrdz,ffdBzdr,ffdBzdz,ffdBphidr,ffdBphidz,ffdns0dr,ffdns0dz,S);
        end
        if(iexit)
            disp(['The ray is out of range! it=',num2str(it)]); return;
        end
        
        [kpar2,kper22,dFdr2,dFdphi2,dFdz2,dFdkr2,dFdnphi2,dFdkz2,dFdw2,D2]=dydt(...
            rtmp,phitmp,ztmp,krtmp,nphitmp,kztmp,...
            B,Br,Bz,Bphi,ns0,dBdr,dBdz,dBrdr,dBrdz,dBzdr,dBzdz,...
            dBphidr,dBphidz,dns0dr,dns0dz,qs,ms,w,c2);
        
        
        % % RK-4, 3rd step
        rtmp=yy(it,1)-dFdkr2/dFdw2*dt;
        phitmp=yy(it,2)-dFdnphi2/dFdw2*dt;
        ztmp=yy(it,3)-dFdkz2/dFdw2*dt;
        krtmp=yy(it,4)+dFdr2/dFdw2*dt;
        nphitmp=yy(it,5)+dFdphi2/dFdw2*dt;
        kztmp=yy(it,6)+dFdz2/dFdw2*dt;
        
        if(numeq==1)
            [jr,jz,psi,B,Br,Bz,Bphi,ns0,ts0,dBdr,dBdz,dBrdr,dBrdz,dBzdr,dBzdz,...
                dBphidr,dBphidz,dns0dr,dns0dz,iexit]=calpars(rtmp,ztmp,...
                rg,zg,dr,dz,fcpsi,fcB,fcBr,fcBz,fcBphi,fcns0,fcts0,fcdBdr,fcdBdz,fcdBrdr,...
                fcdBrdz,fcdBzdr,fcdBzdz,fcdBphidr,fcdBphidz,fcdns0dr,fcdns0dz,S,nr,nz);
        elseif(numeq==0)
            [psi,B,Br,Bz,Bphi,ns0,ts0,dBdr,dBdz,dBrdr,dBrdz,dBzdr,dBzdz,...
                dBphidr,dBphidz,dns0dr,dns0dz,iexit]=calpars_solovev(rtmp,ztmp,1,1,...
                ffpsi,ffB,ffBr,ffBz,ffBphi,ffns0,ffts0,ffdBdr,ffdBdz,ffdBrdr,...
                ffdBrdz,ffdBzdr,ffdBzdz,ffdBphidr,ffdBphidz,ffdns0dr,ffdns0dz,S);
        end
        if(iexit)
            disp(['The ray is out of range! it=',num2str(it)]); return;
        end
        
        [kpar3,kper23,dFdr3,dFdphi3,dFdz3,dFdkr3,dFdnphi3,dFdkz3,dFdw3,D3]=dydt(...
            rtmp,phitmp,ztmp,krtmp,nphitmp,kztmp,...
            B,Br,Bz,Bphi,ns0,dBdr,dBdz,dBrdr,dBrdz,dBzdr,dBzdz,...
            dBphidr,dBphidz,dns0dr,dns0dz,qs,ms,w,c2);
        
        % % RK-4, 4th step
        rtmp=yy(it,1)-dFdkr3/dFdw3*dt;
        phitmp=yy(it,2)-dFdnphi3/dFdw3*dt;
        ztmp=yy(it,3)-dFdkz3/dFdw3*dt;
        krtmp=yy(it,4)+dFdr3/dFdw3*dt;
        nphitmp=yy(it,5)+dFdphi3/dFdw3*dt;
        kztmp=yy(it,6)+dFdz3/dFdw3*dt;
        
        if(numeq==1)
            [jr,jz,psi,B,Br,Bz,Bphi,ns0,ts0,dBdr,dBdz,dBrdr,dBrdz,dBzdr,dBzdz,...
                dBphidr,dBphidz,dns0dr,dns0dz,iexit]=calpars(rtmp,ztmp,...
                rg,zg,dr,dz,fcpsi,fcB,fcBr,fcBz,fcBphi,fcns0,fcts0,fcdBdr,fcdBdz,fcdBrdr,...
                fcdBrdz,fcdBzdr,fcdBzdz,fcdBphidr,fcdBphidz,fcdns0dr,fcdns0dz,S,nr,nz);
        elseif(numeq==0)
            [psi,B,Br,Bz,Bphi,ns0,ts0,dBdr,dBdz,dBrdr,dBrdz,dBzdr,dBzdz,...
                dBphidr,dBphidz,dns0dr,dns0dz,iexit]=calpars_solovev(rtmp,ztmp,1,1,...
                ffpsi,ffB,ffBr,ffBz,ffBphi,ffns0,ffts0,ffdBdr,ffdBdz,ffdBrdr,...
                ffdBrdz,ffdBzdr,ffdBzdz,ffdBphidr,ffdBphidz,ffdns0dr,ffdns0dz,S);
        end
        if(iexit)
            disp(['The ray is out of range! it=',num2str(it)]); return;
        end
        
        [kpar4,kper24,dFdr4,dFdphi4,dFdz4,dFdkr4,dFdnphi4,dFdkz4,dFdw4,D4]=dydt(...
            rtmp,phitmp,ztmp,krtmp,nphitmp,kztmp,...
            B,Br,Bz,Bphi,ns0,dBdr,dBdz,dBrdr,dBrdz,dBzdr,dBzdz,...
            dBphidr,dBphidz,dns0dr,dns0dz,qs,ms,w,c2);
        
        % % RK-4, push
        
        drdt=-(dFdkr1/dFdw1+2.0*dFdkr2/dFdw2+2.0*dFdkr3/dFdw3+dFdkr4/dFdw4)/6.0;
        dphidt=-(dFdnphi1/dFdw1+2.0*dFdnphi2/dFdw2+2.0*dFdnphi3/dFdw3+dFdnphi4/dFdw4)/6.0;
        dzdt=-(dFdkz1/dFdw1+2.0*dFdkz2/dFdw2+2.0*dFdkz3/dFdw3+dFdkz4/dFdw4)/6.0;
        dkrdt=(dFdr1/dFdw1+2.0*dFdr2/dFdw2+2.0*dFdr3/dFdw3+dFdr4/dFdw4)/6.0;
        dnphidt=(dFdphi1/dFdw1+2.0*dFdphi2/dFdw2+2.0*dFdphi3/dFdw3+dFdphi4/dFdw4)/6.0;
        dkzdt=(dFdz1/dFdw1+2.0*dFdz2/dFdw2+2.0*dFdz3/dFdw3+dFdz4/dFdw4)/6.0;
        
        r=yy(it,1)+drdt*dt;
        phi=yy(it,2)+dphidt*dt;
        z=yy(it,3)+dzdt*dt;
        kr=yy(it,4)+dkrdt*dt;
        nphi=yy(it,5)+dnphidt*dt;
        kz=yy(it,6)+dkzdt*dt;
    elseif(irk4==0) % 21-05-17 10:41 add new method to conserve Hamilitonian
        if(numeq==1)
            [jr,jz,psi,B,Br,Bz,Bphi,ns0,ts0,dBdr,dBdz,dBrdr,dBrdz,dBzdr,dBzdz,...
                dBphidr,dBphidz,dns0dr,dns0dz,iexit]=calpars(r,z,...
                rg,zg,dr,dz,fcpsi,fcB,fcBr,fcBz,fcBphi,fcns0,fcts0,fcdBdr,fcdBdz,fcdBrdr,...
                fcdBrdz,fcdBzdr,fcdBzdz,fcdBphidr,fcdBphidz,fcdns0dr,fcdns0dz,S,nr,nz);
        elseif(numeq==0)
            [psi,B,Br,Bz,Bphi,ns0,ts0,dBdr,dBdz,dBrdr,dBrdz,dBzdr,dBzdz,...
                dBphidr,dBphidz,dns0dr,dns0dz,iexit]=calpars_solovev(r,z,1,1,...
                ffpsi,ffB,ffBr,ffBz,ffBphi,ffns0,ffts0,ffdBdr,ffdBdz,ffdBrdr,...
                ffdBrdz,ffdBzdr,ffdBzdz,ffdBphidr,ffdBphidz,ffdns0dr,ffdns0dz,S);
        end
        if(iexit)
            disp(['The ray is out of range! it=',num2str(it)]); return;
        end
        
        [kpar1,kper21,dFdr1,dFdphi1,dFdz1,dFdkr1,dFdnphi1,dFdkz1,dFdw1,D1]=dydt(...
            r,phi,z,kr,nphi,kz,...
            B,Br,Bz,Bphi,ns0,dBdr,dBdz,dBrdr,dBrdz,dBzdr,dBzdz,...
            dBphidr,dBphidz,dns0dr,dns0dz,qs,ms,w,c2);
        
        % % 1st step, only r
        r=yy(it,1)-dFdkr1/dFdw1*dt;
        phi=yy(it,2)-dFdnphi1/dFdw1*dt;
        z=yy(it,3)-dFdkz1/dFdw1*dt;
        
        if(numeq==1)
            [jr,jz,psi,B,Br,Bz,Bphi,ns0,ts0,dBdr,dBdz,dBrdr,dBrdz,dBzdr,dBzdz,...
                dBphidr,dBphidz,dns0dr,dns0dz,iexit]=calpars(r,z,...
                rg,zg,dr,dz,fcpsi,fcB,fcBr,fcBz,fcBphi,fcns0,fcts0,fcdBdr,fcdBdz,fcdBrdr,...
                fcdBrdz,fcdBzdr,fcdBzdz,fcdBphidr,fcdBphidz,fcdns0dr,fcdns0dz,S,nr,nz);
        elseif(numeq==0)
            [psi,B,Br,Bz,Bphi,ns0,ts0,dBdr,dBdz,dBrdr,dBrdz,dBzdr,dBzdz,...
                dBphidr,dBphidz,dns0dr,dns0dz,iexit]=calpars_solovev(r,z,1,1,...
                ffpsi,ffB,ffBr,ffBz,ffBphi,ffns0,ffts0,ffdBdr,ffdBdz,ffdBrdr,...
                ffdBrdz,ffdBzdr,ffdBzdz,ffdBphidr,ffdBphidz,ffdns0dr,ffdns0dz,S);
        end
        if(iexit)
            disp(['The ray is out of range! it=',num2str(it)]); return;
        end
        
        [kpar2,kper22,dFdr2,dFdphi2,dFdz2,dFdkr2,dFdnphi2,dFdkz2,dFdw2,D2]=dydt(...
            r,phi,z,kr,nphi,kz,...
            B,Br,Bz,Bphi,ns0,dBdr,dBdz,dBrdr,dBrdz,dBzdr,dBzdz,...
            dBphidr,dBphidz,dns0dr,dns0dz,qs,ms,w,c2);
        
        % % 2nd step, only k
        kr=yy(it,4)+dFdr2/dFdw2*dt;
        nphi=yy(it,5)+dFdphi2/dFdw2*dt;
        kz=yy(it,6)+dFdz2/dFdw2*dt;
        
        drdt=-dFdkr1/dFdw1;
        dphidt=-dFdnphi1/dFdw1;
        dzdt=-dFdkz1/dFdw1;
        dkrdt=dFdr2/dFdw2;
        dnphidt=dFdphi2/dFdw2;
        dkzdt=dFdz2/dFdw2;
    end
    
    %     D=(D1+2*D2+2*D3+D4)/6.0;
    
    % %
    
    if(r<min(rg))
        kr=-abs(kr);
    end
    %     psilim =   -2.1978970401869986E-004;
    %     if(psilim-psi<0.005*abs(psilim)) % 21-05-13 20:11
    %     if(psilim-psi<-0.004)
    %         kr=-kr;
    %     end
    if(imag(r))
        r=NaN;
    end
    
    yy(it,7)=it*dt;
    yy(it,8)=D1;
    yy(it,9)=drdt;
    yy(it,10)=dphidt;
    yy(it,11)=dzdt;
    yy(it,12)=dkrdt;
    yy(it,13)=dnphidt;
    yy(it,14)=dkzdt;
    yy(it,15)=kpar1;
    yy(it,16)=sqrt(kper21);
    yy(it,17)=B;
    for s=1:S
        yy(it,16+2*s)=ns0(s);
        yy(it,17+2*s)=ts0(s);
    end
    yy(it,18+2*S)=psi;
end
% save([savepath,'ray.mat'],'yy','f','qs','ms');
