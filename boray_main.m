% 2021-03-26 14:15 Huasheng XIE, huashengxie@gmail.com, ENN
% boray.m, ray tracing for 2D axisymmetric equilibrium.
% Cold plasma dispersion relation for ray tracing, hot kinetic plasma for absorb.
% coordinate (r,phi,z,kr,nphi,kz), wave vector k=sqrt(kr^2+kz^2+nphi^2/r^2)
% 21-05-02 14:19 Rewrite to use new interp2d to speed up.
% 21-05-06 09:37 Test OK for single ray.
% 11:53 Test multi-ray OK.
% 21-05-14 08:02 with singularity remove at omega~omega_cs of cold DR
% 21-05-16 07:37 update, treat several species have the same singularity
% 21-05-16 17:05 add Solovev analytical equilibrium profile
% Todo: 1. Current Driven. 2. Singularity of hot DR at resonant point.
%  3. Better analytical profile. 4. Further speed up and improve interp2d.
% 
% Cite: [Xie2021] Hua-sheng XIE, Debabrata Banerjee, Yu-kun BAI, 
%   Han-yue ZHAO and Jing-chun LI, BORAY: An Axisymmetric Ray Tracing Code
%   Supports Both Closed and Open Field Lines Plasmas, 2021, 
%   https://arxiv.org/abs/2105.12014.

close all; clear; clc;

% # intialize the ray
% wave position
icase=1;
% nray=1;
if(icase==1) % load initial data from genray results
    numeq=1; % =0, analytical equilibrium; =1, numerical equilibrium
    
    % # modify here
    eqfile='../eqdata/genray/EAST/genray_eqdata.mat'; % put it in the './input/' directory
    savepath='../output/genray/EAST/'; % '..' means from 'modules' directory
    load('./eqdata/genray/EAST/ray_all.mat');
    
    yyallray0=yyallray;
%     yyallray=yyallray(1:120,:,:);
%     yyallray(:,5:6,:)=-yyallray(:,5:6,:);
    
    ist=20;
%     
%     rayid=[1,2,3];
%     yray0=[squeeze(yyallray(ist,1:6,rayid(1)));
%         squeeze(yyallray(ist,1:6,rayid(2)));
%         squeeze(yyallray(ist,1:6,rayid(3)));];
    rayid=[1];
    yray0=[squeeze(yyallray(ist,1:6,rayid(1)));];
%     yray0=[squeeze(yyallray(ist,1:6,rayid(1)));
%         squeeze(yyallray(ist,1:6,rayid(2)));
%         squeeze(yyallray(ist,1:6,rayid(3)));
%         ];

    dt0=0.0003; nt0=5000;
    
elseif(icase==2) % set the initial data data directly
    numeq=1; % =0, analytical equilibrium; =1, numerical equilibrium
    
    eqfile='../eqdata/mirror/zjs_eqdata.mat'; % put it in the './input/' directory
    savepath='../output/zjs_mirror/'; % '..' means from 'modules' directory

    f=160e6; % Hz
    % [r,phih,z,kr_guess,nphi,kz]
    % change kr_guess for possible multi-modes, i.e., O/X
    yray0=[0.125,0.0,-0.2,900,0.0,-16.3;
        0.125,0.0,0.0,900,0.0,-16.3;
        0.125,0.0,0.2,900,0.0,-16.3;
        0.125,0.0,0.4,900,0.0,-16.3;];
    dt0=0.01; nt0=4000;
    
%     eqfile='../eqdata/frc/c2u_eqdata.mat'; % put it in the './input/' directory
%     savepath='../output/frc/'; % '..' means from 'modules' directory
%     
%     f=7.0e6; % Hz
%     yray0=[0.6,0.0,0.0,-90,-70.0,0.5;
%         0.6,0.0,0.0,-1,-70.0,1.0;
%         0.6,0.0,0.0,-1,-70.0,2.0;];
% %     yray0=[0.6,0.0,0.0,-116,-50.0,3.0;];
%     dt0=0.2; nt0=5000;

elseif(icase==3) % 21-05-16 07:51 for analytical profile
    
    % set Solovev equilibrium for ST, FRC & mirror in initialsolovev.m
    numeq=0;
    
    savepath='../output/analytical/'; % '..' means from 'modules' directory
    iconfig=1; R0=0.64; B0=0.32; q0=1.58; Rx=0.17; E=1.5; tau=0.8;
    qs0=[-1;1;2]; ms0=[1;1836;4*1836];
    ns00=[5.5e18;4.95e18;2.75e17]; % m^-3
    ts00=[200;50;50]; % eV
    
    f=4.5e6; % Hz
%     yray0=[0.8504   -0.1125   -0.1896   -5.4910   -1.6517    0.7289];
    yray0=[0.85   0.0   -0.2   -5.4910   -1.6    0.8];
    
    dt0=0.02; nt0=4000;
    
end
nray=size(yray0,1);
disp('BO-Ray (v210525) seted the equilibrium data for ray tracing.');
%%

for jray=1:nray
%%
    r=yray0(jray,1); phi=yray0(jray,2); z=yray0(jray,3);
    kr_guess=yray0(jray,4); nphi=yray0(jray,5); kz=yray0(jray,6);

run ./modules/boray_initialize;
disp('Finish the initialize. Later to solve the initial kr.');
%%
if(1==0)
    close all;clc;
    % # plot to find the initial kr
    nphi=-1.5; kz=1.0;
    
    fkpar2=@(kr)((kr*Br+kz*Bz+nphi/r*Bphi)/B).^2;
fkper2=@(kr) (kr.^2+kz^2+nphi^2/r^2)-fkpar2(kr);

% define the function to solve the initial kr
fDkr=@(kr)eps1*(fkper2(kr)*c2/w2).^2-...
    ((eps1+eps3)*(eps1-fkpar2(kr)*c2/w2)-eps2^2).*fkper2(kr)*c2/w2+...
    eps3*((eps1-fkpar2(kr)*c2/w2).^2-eps2^2);
    
    krr=-1000:0.5:1000; krr=krr/5;
    subplot(121);plot(krr,fDkr(krr),krr,0.*krr,'--');
    xlabel('k_r');ylabel('D(\omega,k)');
    subplot(122);plot(krr,sqrt(fkpar2(krr)),krr,sqrt(fkper2(krr)),'--');
    xlabel('k_r');ylabel('k');legend('k_{||}','k_\perp');legend('boxoff');
    coeffa=eps1
    eps2
    eps3
    fkpar2(1)*c2/w2
    coeffb=-((eps1+eps3)*(eps1-fkpar2(1)*c2/w2)-eps2^2)
    coeffc=eps3*((eps1-fkpar2(1)*c2/w2).^2-eps2^2)
    roots([coeffa,coeffb,coeffc])
%     ylim([-1,1]*1e1);
    % min(fDkr(krr))
end
%%
% # multi-solutions may exists, choice the correct initial value to
% determine which mode to tracing
options=optimoptions('fsolve','Display','none','TolX',1e-8);
% kr=fsolve(fDkr,-500.0,options) % O-mode
% kr=fsolve(fDkr,-0.1,options) % X-mode
kr=fsolve(fDkr,kr_guess,options)

%%

% # calculate the ray tracing use RK4
tmp=2; % rescaling the dt
nt=nt0*tmp; dt=dt0/c/tmp; % time steps and total time

if(jray==1)
    yyray=zeros(nt,18+2*S,nray)+NaN; %  store for ray tracing
    yrayp=zeros(nt,6+(1+S),nray)+NaN; % store for power absorb
end
runtime1=cputime;
disp(['Begin to calculate the ray tracing use cold',...
    ' plasma dispersion relation.']);
run ./modules/boray_rk4;
yyray(:,:,jray)=yy;
runtime1=cputime-runtime1;
disp(['Finish the ray tracing calculation. Runtime=',...
    num2str(runtime1),'s']);

%% plot the ray
run ./modules/boray_plot;
disp(['The figure is stored in ',savepath]);

%% calculate power
jcalpower=1; % =1, calculate the power absorb; =0, only calculate tracing
if(jcalpower==1)
    jeach=0; % =1, calculate the damping from each species; =0, only all.
    joutw=1; % =1, use matrix method to calculate omega; =0, use dD_r
    N=3; % number of sum_N for kinetic DR. Increase N to make sure convergent
    J=8; % =4,8,12. number of sum_J for Z(zeta) function
    
    runtime2=cputime;
    disp(['Begin to calculate the power absorb along the ray use hot',...
        ' plasma dispersion relation.']);
    run ./modules/boray_calpower;
    
    tmp=[tp,rp,zp,wwh,Ptp];
    yrayp(1:size(tmp,1),1:size(tmp,2),jray)=tmp;
    runtime2=cputime-runtime2;
    disp(['Finish the power absorb calculation. Runtime=',...
        num2str(runtime2),'s']); 
end
%%
if(jcalpower==1)
    run ./modules/boray_plotpower; % 21-05-05 22:18 to update
end
%%
end
