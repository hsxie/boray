% 2021-02-20 15:04 Hua-sheng XIE, ENN, huashengxie@gmail.com
% Plot the cold plasma dispersion relation (N^2,omega), for fixed theta.
% Using Gurnett2005 book Chap 4 notations.

close all; clear; clc;
    
theta=70*pi/180;
% 1. default SI unit
c2=(2.99792458E8)^2; % speed of ligth c^2
epsilon0=8.854187817E-12;
mu0=1/(c2*epsilon0);
kB=1.38064852e-23;
qe=1.60217662e-19; % electron charge, coulombs
mp=1.6726219e-27; % proton mass, kg
me=9.1094e-31; % electron mass, kg

% 2. set parameters, modify here for your own case
B0=0.9; % background B field, Tesla

m_i=mp; % ion mass
m_e=me; % electron mass
q_i=qe; % ion charge
q_e=-qe; % electron charge
n_e=8e17; % electron density, m^-3
n_i=abs(n_e*q_e/q_i); % ion density

ww=10.^(5:0.001:11)*1e1*(2*pi); % omega
nw=length(ww);

% 3. calculate other parameters
wci=q_i*B0/m_i; % the ion gyro frequency
wce=abs(q_e*B0/m_e); % the electron gyro frequency
wpi=sqrt(n_i*q_i^2/(epsilon0*m_i)); % the ion plasma frequency
wpe=sqrt(n_e*q_e^2/(epsilon0*m_e)); % the electron plasma frequency
wp=sqrt(wpi^2+wpe^2); % the total plasma frequency

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
nA2=1+wpi2/wci2+wpe2/wce2;

dielS=1-wpe2./(ww.^2-wce2)-wpi2./(ww.^2-wci2);
dielD=-wce*wpe2./(ww.*(ww.^2-wce2))+wci*wpi2./(ww.*(ww.^2-wci2));
dielP=1-(wpe2+wpi2)./(ww.^2);
dielR=dielS+dielD;
dielL=dielS-dielD;
dielA=dielS*sin(theta)^2+dielP*cos(theta)^2;
dielB=dielR.*dielL*sin(theta)^2+dielP.*dielS*(1+cos(theta)^2);

% the polynomial coefficients
polyc4=dielA;
polyc2=-dielB;
polyc0=dielR.*dielL.*dielP;

% 4. solve the dispersion relation
nn2=zeros(2,nw);
for jw=1:nw
    disppolynomial=[polyc4(jw), polyc2(jw), polyc0(jw)];
    n2temp=sort(real(roots(disppolynomial)),'descend');
    nn2(:,jw)=n2temp;
end

%% plot
close all;
h=figure('unit','normalized','Position',[0.01 0.05 0.6 0.6],...
  'DefaultAxesFontSize',16);

fww=ww/(2*pi); % rad/s -> Hz
% kk=kc/sqrt(c2); % k
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

ytmp=[1e-2,1e6];
loglog([fci,fci],ytmp,'k--',[fce,fce],ytmp,'k--',...
    [flh,flh],ytmp,'k--',[fuh,fuh],ytmp,'m--',[fR,fR],ytmp,'g--',...
    [fL,fL],ytmp,'g--',[fp,fp],ytmp,'m--','linewidth',1);
hold on;
loglog([min(fww),max(fww)],[1,1],'k:',[min(fww),max(fww)],[nA2,nA2],'k:','linewidth',1);
hold on;

nn2p=nn2; nn2p(nn2p<0)=NaN;
nn2m=nn2; nn2m(nn2m>0)=NaN;
nn2p(abs(diff(nn2p,2,2))>5e2)=NaN;
nn2m(abs(diff(nn2m,2,2))>5e2)=NaN;

loglog(fww,abs(nn2p(1,:)),'b-',fww,abs(nn2p(2,:)),'b-','linewidth',3); hold on;
loglog(fww,abs(nn2m(1,:)),'b:',fww,abs(nn2m(2,:)),'b:','linewidth',2); hold on;

xlim([min(fww),max(fww)]);
ylim([1e-2,1e5]);

title(['B_0=',num2str(B0),'T, n_e=',num2str(n_e,3),...
    'm^{-3}, m_i/m_p=',num2str(m_i/mp),', q_i/e=',num2str(q_i/qe),...
    ', \theta=',num2str(theta*180/pi),'\circ']);


text(0.05*max(fww),3e1,...
    ['f_{uh}=',num2str(fuh,3),'Hz',10,...
    'f_{R}=',num2str(fR,3),'Hz',10,...
    'f_{L}=',num2str(fL,3),'Hz',10,...
    'f_{p}=',num2str(fp,3),'Hz',10,...
    'f_{ci}=',num2str(fci,3),'Hz',10,...
    'f_{ce}=',num2str(fce,3),'Hz',10,...
    'f_{lh}=',num2str(flh,3),'Hz'],'fontsize',14);

text(0.05*max(fww),2e4,...
    ['- solid line: n^2>0',10,'. dot line: n^2<0'],'color','b','fontsize',14);

% 
% colorbar;
% 
grid on;grid minor; box on; 
ylabel('|n^2|=|k^2c^2/\omega^2|'); xlabel('f=\omega/2\pi[Hz]');
% 
text(fci,1,'f_{ci}','fontsize',14);
text(fce,1,'f_{ce}','fontsize',14);
text(flh,1,'f_{lh}','fontsize',14);

text(fuh,1e5,'f_{uh}','fontsize',14,'color','m');
text(fp,1e5,'f_{p}','fontsize',14,'color','m');
text(fR,1e-2,'f_{R}','fontsize',14,'color','g');
text(fL,1e-2,'f_{L}','fontsize',14,'color','g');


text(min(fww),nA2,'n_A^2','fontsize',14,'color','k');
% text(min(fww),1,'speed of light','fontsize',14,'color','k');

% save figure
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
  'PaperSize',[screenposition(3:4)]);

print(gcf,'-dpng',['coldwave_N2w_','B0=',num2str(B0),',ne=',num2str(n_e,3),...
    ',mi=',num2str(m_i/mp),',qi=',num2str(q_i/qe),...
    ',theta=',num2str(theta*180/pi),'.png']);

