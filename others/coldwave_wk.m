% 2021-02-19 14:19 Hua-sheng XIE, ENN, huashengxie@gmail.com
% Plot the cold plasma dispersion relation for fixed theta.
% Using Swanson 2-species notation, yields fifth order polynomial. 
% Otherwise, we need use BO/PDRF (Xie2021).
% Ref: [Xie14] H.S. Xie, J. Zhu and Z.W. Ma, Darwin model in plasma physics
%      revisited, Phys. Scr. 89 (2014)105602.
% 21-02-20 14:48 should be ok now

close all; clear; clc;

% tt=[0,10,20,30,50,60,80,85,89,90]*pi/180;
% tt=[0.0001,90]*pi/180;
% tt=[0.0001,30,89.999]*pi/180;
tt=30*pi/180;
for theta=tt
    
% theta=0*pi/180;
% 1. default SI unit
c2=(2.99792458E8)^2; % speed of ligth c^2
epsilon0=8.854187817E-12;
mu0=1/(c2*epsilon0);
kB=1.38064852e-23;
qe=1.60217662e-19; % electron charge, coulombs
mp=1.6726219e-27; % proton mass, kg
me=9.1094e-31; % electron mass, kg

% 2. set parameters, modify here for your own case
B0=1.0/2; % background B field, Tesla

m_i=mp; % ion mass
m_e=me; % electron mass
q_i=qe; % ion charge
q_e=-qe; % electron charge
n_e=10e18; % electron density, m^-3
n_i=abs(n_e*q_e/q_i); % ion density

% kc=(0.001:0.01:10)*1e11; % k*c
kc=10.^(-3:0.1:3)*1e10; % k*c
nkc=length(kc);

% 3. calculate other parameters
wci=q_i*B0/m_i; % the ion gyro frequency
wce=abs(q_e*B0/m_e); % the electron gyro frequency
wpi=sqrt(n_i*q_i^2/(epsilon0*m_i)); % the ion plasma frequency
wpe=sqrt(n_e*q_e^2/(epsilon0*m_e)); % the electron plasma frequency
wp=sqrt(wpi^2+wpe^2); % the total plasma frequency

kc2=kc.*kc;
kc4=kc2.*kc2;
wci2=wci*wci;
wce2=wce*wce;
wpe2=wpe*wpe;
wpi2=wpi*wpi;
wp2=wp*wp;
wextra=wp2+wci*wce;
wextra2=wextra*wextra;
cos2theta=cos(theta).*cos(theta);
wR=wce/2+sqrt(wpe2+wce2/4);
wL=-wce/2+sqrt(wpe2+wce2/4);

% the polynomial coefficients
polyc8=-(2*kc2+wce2+wci2+3*wp2);
polyc6=(kc4+(2*kc2+wp2)*(wce2+wci2+2*wp2)+wextra2);
polyc4=-(kc4*(wce2+wci2+wp2)+2*kc2*wextra2+...
    kc2*wp2*(wce2+wci2-wci*wce).*(1+cos2theta)+wp2*wextra2);
polyc2=(kc4.*(wp2*(wce2+wci2-wci*wce)*cos2theta+wci*wce*wextra)+kc2*wp2* ...
    wci*wce*wextra.*(1+cos2theta));
polyc0=-kc4*wci2*wce2*wp2.*cos2theta;

% 4. solve the dispersion relation
ww=zeros(10,nkc);
for jk=1:nkc
    disppolynomial=[1, 0, polyc8(jk), 0, polyc6(jk), 0, ...
        polyc4(jk), 0, polyc2(jk), 0, polyc0(jk)];
    wtemp=sort(real(roots(disppolynomial)),'descend');
    ww(:,jk)=wtemp;
end

% 5. find some other interesting parameters 
% the elements of the dielectric tensor, using Swansons notation
dielS=1-wpe2./(ww.^2-wce2)-wpi2./(ww.^2-wci2);
dielD=-wce*wpe2./(ww.*(ww.^2-wce2))+wci*wpi2./(ww.*(ww.^2-wci2));
dielP=1-(wpe2+wpi2)./(ww.^2);

kkc2=repmat(kc2,10,1);
% kkx=kkc2*sin(theta); % 21-04-20 11:05 wrong
% kkz=kkc2*cos(theta);
kkx=sqrt(kkc2)*sin(theta);
kkz=sqrt(kkc2)*cos(theta);

n2=kkc2./(ww.^2);

dielxx=dielS-n2.*cos(theta).*cos(theta);
dielxy=-1i*dielD;
dielxz=n2.*cos(theta).*sin(theta)+1e-10; % 21-02-20 13:28, to avoid =0
dielyy=dielS-n2;
dielzz=dielP-n2.*sin(theta).*sin(theta);

Ex=-dielzz./dielxz;
Ey=dielxy./dielyy.*Ex;
Ez=1;
Eperp=sqrt(Ex.*conj(Ex)+Ey.*conj(Ey));
Etot=sqrt(Ex.*conj(Ex)+Ey.*conj(Ey)+1);
EparK=(kkx.*Ex+kkz.*Ez)./sqrt(kkc2);
Epolar=-2*imag(Ex.*conj(Ey))./(Eperp.*Eperp+0e-10); % to avoid Inf

Bx=-kkz.*Ey./ww; % to update with correct unit c, epsilon0, mu0
By=(kkz.*Ex-kkx.*Ez)./ww;
Bz=kkx.*Ey./ww;
Btot=sqrt(Bx.*conj(Bx)+By.*conj(By)+Bz.*conj(Bz));

% temp=length(kc_x);
% dk_x=kc_x(2);dk_z=kc_z(2);
% dw_x=diff(ww,1,2); dw_z=diff(ww,1,3);
% dw_x(1,temp,temp)=0; dw_z(1,temp,temp)=0;
% v_x=dw_x/dk_x; v_z=dw_z/dk_z;
% pola(5,:,:,:)=sqrt(v_x.*v_x+v_z.*v_z); % value of the group vel.

pola(1,:,:)=kkx;
pola(2,:,:)=kkz;
pola(3,:,:)=ww;
pola(4,:,:)=Btot./Etot; % degree of electromagnetism
pola(5,:,:)=abs(EparK)./Etot; % degree of longitudinality
pola(6,:,:)=Ez./Etot; % degree of parallelity
pola(7,:,:)=Epolar; % ellipticity


%% plot
close all;
h=figure('unit','normalized','Position',[0.01 0.05 0.85 0.7],...
  'DefaultAxesFontSize',16);

fww=ww/(2*pi); % rad/s -> Hz
kk=kc/sqrt(c2); % k
fci=wci/(2*pi);
fce=wce/(2*pi);
fpi=wpi/(2*pi);
fpe=wpe/(2*pi);
fp=wp/(2*pi);
% flh=sqrt(fce*fci); % LHW
flh=1/sqrt(1/(fci^2+fpi^2)+1/(fce*fci)); % LHW
fuh=sqrt(fce^2+fpe^2); % UHW
fR=wR/(2*pi); % right cutoff
fL=wL/(2*pi); % left cutoff

% semilogy(kc,ww(1,:),kc,ww(2,:),kc,ww(3,:),kc,ww(4,:),kc,ww(5,:));
% loglog(kk,fww(1,:),kk,fww(2,:),kk,fww(3,:),kk,fww(4,:),kk,fww(5,:),...
%     'linewidth',3);hold on;

for jsub=1:2
    subplot(1,2,jsub);
loglog(kk,0.*kk+fci,'k--',kk,0.*kk+fce,'k--',kk,0.*kk+flh,'k--',...
    kk,kc/(2*pi),'k--',kk,0.*kk+fuh,'k--','linewidth',1);
hold on;

loglog(kk,0.*kk+fR,'m:',kk,0.*kk+fL,'m:',kk,0.*kk+fp,'m:','linewidth',1);
hold on;

ymin=1e5;ymax=1e14;
if(jsub==1)
    
patch([kk,kk(end)],[fww(1,:),NaN],[squeeze(pola(6,1,:));NaN],...
    'EdgeColor','interp','MarkerFaceColor',...
    'flat','FaceColor','none','LineWidth',3,'LineStyle','-'); hold on;
patch([kk,kk(end)],[fww(2,:),NaN],[squeeze(pola(6,2,:));NaN],...
    'EdgeColor','interp','MarkerFaceColor',...
    'flat','FaceColor','none','LineWidth',3,'LineStyle','-'); hold on;
patch([kk,kk(end)],[fww(3,:),NaN],[squeeze(pola(6,3,:));NaN],...
    'EdgeColor','interp','MarkerFaceColor',...
    'flat','FaceColor','none','LineWidth',3,'LineStyle','-'); hold on;
patch([kk,kk(end)],[fww(4,:),NaN],[squeeze(pola(6,4,:));NaN],...
    'EdgeColor','interp','MarkerFaceColor',...
    'flat','FaceColor','none','LineWidth',3,'LineStyle','-'); hold on;
patch([kk,kk(end)],[fww(5,:),NaN],[squeeze(pola(6,5,:));NaN],...
    'EdgeColor','interp','MarkerFaceColor',...
    'flat','FaceColor','none','LineWidth',3,'LineStyle','-'); hold on;
% title(['B_0=',num2str(B0),'T, n_e=',num2str(n_e,3),...
%     'm^{-3}, m_i/m_p=',num2str(m_i/mp),', q_i/e=',num2str(q_i/qe),...
%     ', \theta=',num2str(theta*180/pi),'\circ']);

text(10*min(kk),0.01*ymax,...
    ['B_{0}=',num2str(B0),'T',10,...
    'n_{e}=',num2str(n_e),'m^{-3}',10,...
    'm_i/m_p=',num2str(m_i/mp),10,...
    'q_i/e=',num2str(q_i/qe),10,...
    '\theta=',num2str(theta*180/pi),'\circ'],'fontsize',14);

text(2*min(kk),0.5*ymax,['(a) E_z/E_{tot}, degree of parallelity,',...
    ' X/O for 0/1'],'fontsize',14);
else
    
patch([kk,kk(end)],[fww(1,:),NaN],[squeeze(pola(7,1,:));NaN],...
    'EdgeColor','interp','MarkerFaceColor',...
    'flat','FaceColor','none','LineWidth',3,'LineStyle','-'); hold on;
patch([kk,kk(end)],[fww(2,:),NaN],[squeeze(pola(7,2,:));NaN],...
    'EdgeColor','interp','MarkerFaceColor',...
    'flat','FaceColor','none','LineWidth',3,'LineStyle','-'); hold on;
patch([kk,kk(end)],[fww(3,:),NaN],[squeeze(pola(7,3,:));NaN],...
    'EdgeColor','interp','MarkerFaceColor',...
    'flat','FaceColor','none','LineWidth',3,'LineStyle','-'); hold on;
patch([kk,kk(end)],[fww(4,:),NaN],[squeeze(pola(7,4,:));NaN],...
    'EdgeColor','interp','MarkerFaceColor',...
    'flat','FaceColor','none','LineWidth',3,'LineStyle','-'); hold on;
patch([kk,kk(end)],[fww(5,:),NaN],[squeeze(pola(7,5,:));NaN],...
    'EdgeColor','interp','MarkerFaceColor',...
    'flat','FaceColor','none','LineWidth',3,'LineStyle','-'); hold on;

text(10*min(kk),0.012*ymax,...
    ['f_{uh}=',num2str(fuh,3),'Hz',10,...
    'f_{R}=',num2str(fR,3),'Hz',10,...
    'f_{L}=',num2str(fL,3),'Hz',10,...
    'f_{p}=',num2str(fp,3),'Hz'],'fontsize',14);

text(0.01*max(kk),10*ymin,...
    ['f_{ci}=',num2str(fci,3),'Hz',10,...
    'f_{ce}=',num2str(fce,3),'Hz',10,...
    'f_{lh}=',num2str(flh,3),'Hz'],'fontsize',14);

text(2*min(kk),0.5*ymax,'(b) ellipticity, L/R for -1/1','fontsize',14);
end

colorbar;

grid on;grid minor; box on; xlim([min(kk),max(kk)]);ylim([ymin,ymax]);%axis tight;
xlabel('k=2\pi/\lambda[m^{-1}]'); ylabel('f=\omega/2\pi[Hz]');

text(0.6*max(kk),fci,'f_{ci}','fontsize',14);
text(0.6*max(kk),fce,'f_{ce}','fontsize',14);
text(0.6*max(kk),flh,'f_{lh}','fontsize',14);
text(0.6*max(kk),fuh,'f_{uh}','fontsize',14);
text(min(kk),fp,'f_{p}','fontsize',14);
text(min(kk),fR,'f_{R}','fontsize',14);
text(min(kk),fL,'f_{L}','fontsize',14);
text(0.2*max(kk),max(kc/2/pi),'kc/2\pi','fontsize',14);
end

% save figure
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
  'PaperSize',[screenposition(3:4)]);

print(gcf,'-dpng',['coldwave_','B0=',num2str(B0),',ne=',num2str(n_e,3),...
    ',mi=',num2str(m_i/mp),',qi=',num2str(q_i/qe),...
    ',theta=',num2str(theta*180/pi),'.png']);

end




