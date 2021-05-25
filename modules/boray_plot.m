% 2021-05-05 14:40 Huasheng XIE, huashengxie@gmail.com, ENN
% boray_plot.m

close all;
h=figure('unit','normalized','Position',[0.01 0.05 0.6 0.6],...
    'DefaultAxesFontSize',14);

subplot(231);
if(icase==1)
    contour(rr,zz,fpsi,100); hold on;
    plot(yy(:,1),yy(:,3),yy(1,1),yy(1,3),'rx','linewidth',2); hold on;
    yyid=yyallray(ist:end,:,rayid(jray)); % 21-05-06 11:13 to update
    plot(yyid(:,1),yyid(:,3),'r:','linewidth',1); hold on;
elseif(icase==2 || icase==3)
    contour(rr,zz,fpsi,100); hold on;
    plot(yy(:,1),yy(:,3),yy(1,1),yy(1,3),'rx','linewidth',2); hold on;
    
else
    surf(rr,zz,fB,'linewidth',2);
    zlabel('B[T]');
end

if(1==0)
    
    Bn1=f/(28e9/(max(abs(ms./qs))/min(abs(ms./qs))));Bn1=Bn1(1);
    [Cc1,hc1]=contour(rr,zz,fB,[1,1]*Bn1,'visible','off');
    [Cc2,hc2]=contour(rr,zz,fB,[1/2,1/2]*Bn1,'visible','off');
    %         [Cc3,hc3]=contour(rr,zz,fB,[1/3,1/3]*Bn1); hold on;
    % xlabel('z');ylabel('r');
    % title('(d) n_s(r,z)');
    %         contour(zz,rr,squeeze(fns0(1,:,:)),100); hold on;
    plot(Cc1(1,2:Cc1(2,1)),Cc1(2,2:Cc1(2,1)),'--','linewidth',1); hold on;
%     plot(Cc2(1,2:Cc2(2,1)),Cc2(2,2:Cc2(2,1)),'--','linewidth',1); hold on;
    %         plot(Cc3(2,2:end),Cc3(1,2:end),'.','linewidth',1); hold on;
end

xlabel('R');ylabel('Z');
title(['ray runtime=',num2str(runtime1),'s']);

subplot(232);
contour(rr,zz,squeeze(fns0(1,:,:)),50); colorbar;
xlabel('R');ylabel('Z'); title(['n_s[m^{-3}], S=',num2str(S)]);

subplot(233);
plot(yy(:,1).*cos(yy(:,2)),yy(:,1).*sin(yy(:,2)),...
    yy(1,1).*cos(yy(1,2)),yy(1,1).*sin(yy(1,2)),'rx',...
    'linewidth',2);hold on;
if(icase==1)
    plot(yyid(:,1).*cos(yyid(:,2)),yyid(:,1).*sin(yyid(:,2)),'r:',...
        'linewidth',1);hold on;
    legend('boray','genray','location','best');legend('boxoff');
end
tt=0:0.01*pi:2*pi;
plot(min(rg)*cos(tt),min(rg)*sin(tt),max(rg)*cos(tt),max(rg)*sin(tt));
axis equal;
xlabel('X');ylabel('Y');

title(['r=',num2str(yy(1,1),3),...
    ', \phi=',num2str(yy(1,2),3),', z=',num2str(yy(1,3),3)]);

subplot(234);
plot(yy(:,1),yy(:,8),'linewidth',2); % should fDR=0
xlabel('R');ylabel('D(\omega,k)');

subplot(235);
plot(yy(:,7),yy(:,4),yy(:,7),yy(:,5),yy(:,7),yy(:,6),...
    'linewidth',2); hold on;
xlabel('t');ylabel('k');
legend('k_r','n_\phi','k_z','location','best');
legend('boxoff');
title(['n0=',num2str(n0,3),...
    'm^{-3}, f=',num2str(f/1e6,4),'MHz, k_r=',num2str(yy(1,4),3),...
    ', n_\phi=',num2str(yy(1,5),3),', k_z=',num2str(yy(1,6),3),...
    ', dt=',num2str(dt,2),', nt=',num2str(nt)]);

subplot(236);
if(icase==1)
    yyidt=1:length(yyid(:,1));
    plot(yyidt,yyid(:,4),yyidt,yyid(:,5),yyidt,yyid(:,6),...
        'linewidth',2); hold on;
    legend('k_r','n_\phi','k_z','location','best');
    xlabel('t');ylabel('k'); title('genray');
    legend('boxoff');
else
    plot(yy(:,7),yy(:,9),yy(:,7),yy(:,10).*yy(:,1),...
        yy(:,7),yy(:,11),'linewidth',2); hold on;
    xlabel('t');ylabel('v_g');
    legend('vg_r','vg_\phi','vg_z','location','best');
    legend('boxoff');
end

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);

strtmp='';

print(gcf,'-dpng',[savepath,'boray_f=',num2str(f/1e9,3),...
    'GHz,y0=',num2str(yy(1,1),3),',',num2str(yy(1,2),3),',',...
    num2str(yy(1,3),3),',',num2str(yy(1,4),3),',',num2str(yy(1,5),3),...
    ',',num2str(yy(1,6),3),',dt=',num2str(dt,2),',nt=',num2str(nt),...
    ',S=',num2str(S),',jray=',num2str(jray),strtmp,'.png']);

%%
if(jray==nray)
h=figure('unit','normalized','Position',[0.01 0.05 0.3 0.6],...
    'DefaultAxesFontSize',14);

if(icase==1)
    contour(rr,zz,fpsi,100); hold on;
    
    for j=1:nray
        plot(yyray(:,1,j),yyray(:,3,j),...
            yyray(1,1,j),yyray(1,3,j),'rx','linewidth',2); hold on;
        
        yyid=yyallray(ist:end,:,rayid(j)); % 21-05-06 11:13 to update
        plot(yyid(:,1),yyid(:,3),'r:','linewidth',1); hold on;
    end
    xlabel('R'); ylabel('Z');
    
elseif(icase==2 || icase==3)
    
    contour(rr,zz,fpsi,100); hold on;
    for j=1:nray
        plot(yyray(:,1,j),yyray(:,3,j),...
            yyray(1,1,j),yyray(1,3,j),'rx','linewidth',2); hold on;
    end
    xlabel('R'); ylabel('Z');
    title(['f=',num2str(f/1e6,4),'MHz']);
    
end

if(1==0)
    
    Bn1=f/(28e9/(max(abs(ms./qs))/min(abs(ms./qs))));Bn1=Bn1(1);
    [Cc1,hc1]=contour(rr,zz,fB,[1,1]*Bn1,'visible','off');
    [Cc2,hc2]=contour(rr,zz,fB,[1/2,1/2]*Bn1,'visible','off');
    %         [Cc3,hc3]=contour(rr,zz,fB,[1/3,1/3]*Bn1); hold on;
    % xlabel('z');ylabel('r');
    % title('(d) n_s(r,z)');
    %         contour(zz,rr,squeeze(fns0(1,:,:)),100); hold on;
    plot(Cc1(1,2:Cc1(2,1)),Cc1(2,2:Cc1(2,1)),'--','linewidth',1); hold on;
%     plot(Cc2(1,2:Cc2(2,1)),Cc2(2,2:Cc2(2,1)),'--','linewidth',1); hold on;
    %         plot(Cc3(2,2:end),Cc3(1,2:end),'.','linewidth',1); hold on;
end
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print(gcf,'-dpng',[savepath,'boray_f=',num2str(f/1e9,3),...
    'GHz,nray=',num2str(nray),',dt=',num2str(dt,2),',nt=',num2str(nt),...
    ',S=',num2str(S),strtmp,'.png']);
end

%%
if(1==0)
figure;
plot(yyray(:,7),yyray(:,9),yyray(:,7),yyray(:,10),yyray(:,7),yyray(:,11),'linewidth',2);
xlabel('t');legend('v_{gr}','v_{g\phi}','v_{gz}');legend('boxoff');
end