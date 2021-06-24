% 2021-02-25 13:06 Hua-sheng XIE, ENN, huashengxie@gmail.com
% Plot the cold plasma dispersion relation for CMA diagram.

close all; clear; clc;

Mi=1836; % Mi=e*mi/(qi*me)

x=10.^(-3:0.00001:8);
y=10.^(-2:0.00001:4);

xR=1./(1./(1-y)+1./(Mi+y)); % R=0
xL=1./(1./(1+y)+1./(Mi-y)); % L=0
xP=1/(1+1/Mi)+0.*y; % P=0
xS=(1./(1./(1-y.^2)+Mi./(Mi^2-y.^2))); % S=0
yRinf=1+0.*x; % R=inf
yLinf=Mi+0.*x; % L=inf

xRL_PS=((1./(1-y)+1./(Mi+y))+(1./(1+y)+1./(Mi-y))-(1+1/Mi)-(1./(1-y.^2)+Mi./(Mi^2-y.^2))...
    )./((1./(1-y)+1./(Mi+y)).*(1./(1+y)+1./(Mi-y))-(1+1/Mi).*(1./(1-y.^2)+Mi./(Mi^2-y.^2)));

% fRL=(1-(1./(1-y)+1./(Mi+y)).*xRL_PS).*(1-(1./(1+y)+1./(Mi-y)).*xRL_PS);
% fPS=(1-(1+1/Mi).*xRL_PS).*(1-(1./(1-y.^2)+Mi./(Mi^2-y.^2)).*xRL_PS);

xR(xR<0)=NaN;
xL(xL<0)=NaN;
xS(xS<0)=NaN;
xRL_PS(xRL_PS<0)=NaN;

%%
close all;
h=figure('unit','normalized','Position',[0.01 0.05 0.5 0.7],...
  'DefaultAxesFontSize',16);

ax1=axes('Position',[0.1,0.1,0.8,0.8]);

loglog(xR,y,xL,y,xS,y,xRL_PS,y,'-.',xP,y,'--',x,yRinf,'--',x,yLinf,'--','linewidth',3);
legend('R=0','L=0','S=0 (hybrid resonance)','RL=PS',...
    'P=0','R=\infty (\omega_{ce} resonance)',...
    'L=\infty (\omega_{ci} resonance)','location','best');
legend('boxoff');
xlim([min(x),max(x)]);
ylim([min(y),max(y)]);
xlabel('x=\omega_{pe}^2/\omega^2');ylabel('y=|\omega_{ce}|/\omega');
text(min(x),Mi,['M_i=',num2str(Mi)],'fontsize',14);

% 
c2=(2.99792458E8)^2; % speed of ligth c^2
epsilon0=8.854187817E-12;
mu0=1/(c2*epsilon0);
kB=1.38064852e-23;
qe=1.60217662e-19; % electron charge, coulombs
mp=1.6726219e-27; % proton mass, kg
me=9.1094e-31; % electron mass, kg

% 2. set parameters, modify here for your own case
m_i=mp; % ion mass
m_e=me; % electron mass
q_i=qe; % ion charge
q_e=-qe; % electron charge

fw=28e9;  % Hz, wave frequency
% fw=13.56e6;
% fw=2.45e9;
w=fw*(2*pi); % wave frequency

Bmin=0.2; Bmax=1.2;
nemin=5e15; nemax=8e18;
xmin=(nemin*q_e^2/(epsilon0*m_e))/w^2;
xmax=(nemax*q_e^2/(epsilon0*m_e))/w^2;
ymin=abs(q_e*Bmin/m_e)/w;
ymax=abs(q_e*Bmax/m_e)/w;
hold on;
loglog([xmin,xmax],[ymin,ymax],':','linewidth',2);
% loglog(xmin,ymin,'.',xmax,ymax,'.',[xmin,xmax],[ymin,ymax],':','linewidth',2);
hold on;
patch([xmin,xmax,xmax,xmin,xmin],[ymin,ymin,ymax,ymax,ymin],'y-.'); alpha(0.05);
text(xmax,ymax,'max');
text(xmin,ymin,'min');
% annotation('textarrow',[xmin,xmax],[ymin,ymax],'String','parameter range');
title(['f=',num2str(fw,3),'Hz, n_{e}=[',num2str(nemin,3),...
    ',',num2str(nemax,3),...
    ']m^{-3}, B=[',num2str(Bmin),',',num2str(Bmax),']T']);

%%
tt=0:0.005*pi:2*pi;
jtmp=floor((length(tt)-1)/4)+1;
% xyjj=[xmin,xmax;
%     ymin,ymax];
xyjj=[xmin,xmax,1e-2,1e2,1e4,1e-2,1e5,2e1,2e-2,0.63,0.5,1.5,1e2;
    ymin,ymax,5e-2,1e-1,4e0,1e1,3e2,1e2,4e3,0.4,0.9,0.8,4e3];
for jj=1:size(xyjj,2)
%     xp=xmax; yp=ymax;
    xp=xyjj(1,jj);
    yp=xyjj(2,jj);
    
    Rp=1-xp/(1-yp)-xp/(Mi+yp);
    Lp=1-xp/(1+yp)-xp/(Mi-yp);
    Pp=1-xp*(1+1/Mi);
    Sp=1-xp/(1-yp^2)-Mi*xp/(Mi^2-yp^2);
    Dp=-yp*xp/(1-yp^2)+yp*xp/(Mi^2-yp^2);
    
    Apt=Sp*sin(tt).^2+Pp*cos(tt).^2;
    Bpt=Rp*Lp*sin(tt).^2+Pp*Sp*(1+cos(tt).^2);
    Fpt=sqrt(Bpt.^2-4*Apt*Rp*Lp*Pp);
    n2ta=(Bpt+Fpt)./(2*Apt);
    n2tb=(Bpt-Fpt)./(2*Apt);
    n2ta(real(n2ta)<=0)=NaN;
    n2tb(real(n2tb)<=0)=NaN;

    for jab=1:2
        
        if(jab==1)
            n2t=n2ta;
        else
            n2t=n2tb;
        end
        dielxx=Sp-n2t.*cos(tt).^2;
        dielxy=-1i*Dp;
        dielxz=n2t.*cos(tt).*sin(tt)+1e-10;
        dielyy=Sp-n2t;
        dielzz=Pp-n2t.*sin(tt).^2;
        
        Ez=1;
        if(jab==1)
            Exa=-dielzz./dielxz;
            Eya=dielxy./dielyy.*Exa;
            RLa=sign(-1i*Eya(1)/Exa(1));
            if(RLa>0)
                jRLa=1;
            elseif(RLa<0)
                jRLa=2;
            else
                jRLa=3;
            end
            tmp=sqrt(Exa.*conj(Exa)+Eya.*conj(Eya))./Ez;
%             OXa=(tmp(jtmp)<1);
            if(isnan(tmp(jtmp)))
                jOXa=3;
            elseif(tmp(jtmp)>1)
                jOXa=2;
            else
                jOXa=1;
            end
        else
            Exb=-dielzz./dielxz;
            Eyb=dielxy./dielyy.*Exb;
            RLb=sign(-1i*Eyb(1)/Exb(1));            
            if(RLb>0)
                jRLb=1;
            elseif(RLb<0)
                jRLb=2;
            else
                jRLb=3;
            end
            tmp=sqrt(Exb.*conj(Exb)+Eyb.*conj(Eyb))./Ez;
%             OXb=(tmp(jtmp)<1);
            if(isnan(tmp(jtmp)))
                jOXb=3;
            elseif(tmp(jtmp)>1)
                jOXb=2;
            else
                jOXb=1;
            end
        end
        strRL='RL ';
        strOX='OX ';
        
%         Eperp=sqrt(Ex.*conj(Ex)+Ey.*conj(Ey));
%         Etot=sqrt(Ex.*conj(Ex)+Ey.*conj(Ey)+1);
%         EparK=(kkx.*Ex+kkz.*Ez)./sqrt(kkc2);
%         Epolar=-2*imag(Ex.*conj(Ey))./(Eperp.*Eperp+0e-10); % to avoid Inf
    end


%     axes('Position',[0.3+0.12*jj,0.3+0.12*jj,0.1,0.1]);
%     axes('Position',[0,0+0.12*(jj-1),0.1,0.1]);
    
xtmp1=0.1+0.8*(log10(xp)-log10(min(x)))/(log10(max(x))-log10(min(x)));
ytmp1=0.1+0.8*(log10(yp)-log10(min(y)))/(log10(max(y))-log10(min(y)));
    axes('Position',[xtmp1-0.02,ytmp1-0.02,0.04,0.04]);
    
    vgxa=sqrt(1./n2ta).*cos(tt+0.5*pi); 
    vgya=sqrt(1./n2ta).*sin(tt+0.5*pi);
    vgxb=sqrt(1./n2tb).*cos(tt+0.5*pi); 
    vgyb=sqrt(1./n2tb).*sin(tt+0.5*pi);
%     plot(sqrt(n2ta).*cos(tt+0.5*pi),sqrt(n2ta).*sin(tt+0.5*pi),'r',...
%         sqrt(n2tb).*cos(tt+0.5*pi),sqrt(n2tb).*sin(tt+0.5*pi),'b','linewidth',2);
    plot(vgxa,vgya,'r',vgxb,vgyb,'b','linewidth',1);
    hold on;
    
%     figure;
%     patch([sqrt(n2ta).*cos(tt+0.5*pi),NaN],[sqrt(n2ta).*sin(tt+0.5*pi),NaN],[real(-1i*Eya./Exa),1],...
%         'EdgeColor','interp','MarkerFaceColor',...
%     'flat','FaceColor','none','LineWidth',1,'LineStyle','-');
%     hold on;
%     patch([sqrt(n2tb).*cos(tt+0.5*pi),NaN],[sqrt(n2tb).*sin(tt+0.5*pi),NaN],[real(-1i*Eyb./Exb),1],...
%         'EdgeColor','interp','MarkerFaceColor',...
%     'flat','FaceColor','none','LineWidth',1,'LineStyle','-');
    text(sqrt(1./n2ta(1)).*cos(tt(1)+0.5*pi),sqrt(1./n2ta(1)).*sin(tt(1)+0.5*pi),strRL(jRLa),'Color','r');
    text(sqrt(1./n2tb(1)).*cos(tt(1)+1.5*pi),sqrt(1./n2tb(1)).*sin(tt(1)+1.5*pi),strRL(jRLb),'Color','b');
    text(sqrt(1./n2ta(jtmp)).*cos(tt(jtmp)+0.5*pi),sqrt(1./n2ta(jtmp)).*sin(tt(jtmp)+0.5*pi),strOX(jOXa),'Color','r');
    text(sqrt(1./n2tb(jtmp)).*cos(tt(jtmp)+1.5*pi),sqrt(1./n2tb(jtmp)).*sin(tt(jtmp)+1.5*pi),strOX(jOXb),'Color','b');
    axis off; 
    axis equal;
    
end


% save figure
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
  'PaperSize',[screenposition(3:4)]);

print(gcf,'-dpng',['pltcma_','Mi=',num2str(Mi),...
    'fw=',num2str(fw,3),',B=',num2str(Bmin),...
    ',',num2str(Bmax),',ne=',num2str(nemin,3),',',num2str(nemin,3),'.png']);


