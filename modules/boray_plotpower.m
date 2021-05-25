% 2021-05-05 15:26 Huasheng XIE, huashengxie@gmail.com, ENN
% boray_plotpower.m
% 22:18 to update

h=figure('unit','normalized','Position',[0.01 0.05 0.6 0.7],...
  'DefaultAxesFontSize',14);

subplot(231);
plot(tp,real(ww),tp,real(wwh(:,1)),'-','linewidth',2);
xlabel('t');ylabel('\omega_r');
ylim([0,1.5*max(real(ww))]);
legend('cold DR','hot DR','location','best');legend('boxoff');
for js=2:N_each
    hold on; plot(tp,real(wwh(:,js)),'--','linewidth',2);
end


subplot(232);
plot(tp,imag(wwh(:,1)),'-','linewidth',2);
if(N_each>=2)
for js=2:N_each
    hold on; plot(tp,imag(wwh(:,js)),'--','linewidth',2);
end
hold on; plot(tp,sum(imag(wwh(:,2:end)),2),':','linewidth',2);
legend('all','only s=1','only s=2','only s=3','sum','location','best');
legend('boxoff');
end

xlabel('t');ylabel('\omega_i');
title(['joutw=',num2str(joutw),', max(Te)=',num2str(max(Ts0(:,1)),3),'eV']);


subplot(233);
if(icase==1)
    plot(rp,1-Ptp(:,1),'linewidth',2); hold on; xlabel('R');
    
%     plot(yyid(1:80,1),yyid(1:80,end)/yyid(1,end),'--','linewidth',2);
    plot(yyid(:,1),yyid(:,end)/yyid(1,end),'--','linewidth',2);
    legend('boray','genray','loaction','best');legend('boxoff');
%     if(jeach==1)
%         damp2=2*cumsum(imag(squeeze(wwh(:,2)))*dt);
%         P2=exp(damp2);
%         hold on;
%         plot(rp,P2,':','linewidth',2); hold on;
%         legend('bo-ray','genray','bo-ray, only electron','location','best');legend('boxoff');
%     end
    
else
%     plot(zp,1-Ptp(:,1),'linewidth',2); hold on; xlabel('Z');
    plot(rp,1-Ptp(:,1),'linewidth',2); hold on; xlabel('R');
end
% ylim([0,1.1]);
ylabel('Power Absorb');


subplot(234);
psim=abs(min(min(fpsi)));
contour(rr,zz,fpsi,100); hold on;
% contour(zz,rr,fpsi,[0,0],'r--'); hold on;
plot(yy(1:end,1),yy(1:end,3),'-',yy(1,1),yy(1,3),'rx','linewidth',2); hold on;
xlabel('R');ylabel('Z');title('\psi and ray tracing');
axis tight; %colorbar;
if(icase==0)
patch([rp;NaN],[zp;zp(end)],[1-Ptp(:,1);NaN],...
    'EdgeColor','interp','MarkerFaceColor',...
    'flat','FaceColor','none','LineWidth',3,'LineStyle','-'); hold on;
end
% plot(yyid(1:100,1),yyid(1:100,3),'r--','linewidth',1); hold on;
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

subplot(235);
contour(rr,zz,fB,100);
xlabel('R');ylabel('Z');%title('B');
axis tight; colorbar;
title(['B, f=',num2str(f/1e6,4),'MHz, S=',num2str(S)]);

subplot(236);
contour(rr,zz,squeeze(fns0(1,:,:)),100); hold on;

plot(yyp(1:end,1),yyp(1:end,3),'-',yyp(1,1),yy(1,3),'rx',...
    'linewidth',2); hold on;
if(icase==1)
    plot(yyid(1:end,1),yyid(1:end,3),'r--','linewidth',2); hold on;
end
xlabel('R');ylabel('Z');title('n_{e}');
colorbar('eastoutside');

xlim([0,max(max(rr))]);ylim([zz(1,1),zz(end,end)]);
% axis equal;

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
  'PaperSize',[screenposition(3:4)]);

print(gcf,'-dpng',[savepath,'boray_power_f=',num2str(w(1,1)/(2*pi)/1e6,3),...
    'MHz,y0=',num2str(yy(1,1),3),',',num2str(yy(1,2),3),',',...
    num2str(yy(1,3),3),',',num2str(yy(1,4),3),',',num2str(yy(1,5),3),...
    ',',num2str(yy(1,6),3),',dt=',num2str(dt,2),',N=',num2str(N),',J=',num2str(J),'.png']);
%%
if(jray>nray) % 21-05-06 11:48 to plot all ray power absorb
h=figure('unit','normalized','Position',[0.01 0.05 0.6 0.6],...
    'DefaultAxesFontSize',14);

if(icase==1)
    contour(rr,zz,fpsi,100); hold on;
    
%     for j=1:nray
%         plot(yyray(:,1,j),yyray(:,3,j),...
%             yyray(1,1,j),yyray(1,3,j),'rx','linewidth',2); hold on;
%         
%         yyid=yyallray(ist:end,:,j); % 21-05-06 11:13 to update
%         plot(yyid(:,1),yyid(:,3),'r:','linewidth',1); hold on;
%     end
%     xlabel('R'); ylabel('Z');
    
elseif(icase==2)
    
%     tmp=[tp,rp,zp,wwh,Ptp];
%     yrayp(1:size(tmp,1),1:size(tmp,2),jray)=tmp;
    subplot(121);
    contour(rr,zz,fpsi,100); hold on;
    xlabel('R'); ylabel('Z');
    title(['f=',num2str(f/1e9,4),'GHz, N_\phi=',...
        num2str(yray0(1,5)./yray0(1,1)*c/w,3),...
        ', N_z=',num2str(yray0(1,6)*c/w,3)]);
    
    for j=1:nray
%         patch([squeeze(yrayp(:,2,j));NaN],[squeeze(yrayp(:,3,j));...
%             squeeze(yrayp(end,3,j))],...
%             [1-squeeze(yrayp(:,3+N_each+1,j));NaN],...
%             'EdgeColor','interp','MarkerFaceColor',...
%             'flat','FaceColor','none','LineWidth',3,'LineStyle','-'); hold on;

        plot(yyray(:,1,j),yyray(:,3,j),...
            yyray(1,1,j),yyray(1,3,j),'rx','linewidth',2); hold on;
    end
    
    subplot(122);
    for j=1:nray
%         plot(yrayp(:,1,j),1-yrayp(:,3+N_each+1,j),'linewidth',2); hold on;
        plot(yrayp(:,2,j),1-yrayp(:,3+N_each+1,j),'linewidth',2); hold on;
    end
    xlabel('R'); ylabel('Power');
    title(['N_{r1}=',num2str(yyray(1,4,1)*c/w,3),...
        ', N_{r2}=',num2str(yyray(1,4,2)*c/w,3)]);
end

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print(gcf,'-dpng',[savepath,'boray_f=',num2str(f/1e9,3),...
    'GHz,nray=',num2str(nray),',dt=',num2str(dt,2),',nt=',num2str(nt),...
    strtmp,'.png']);
end

%%
if(1==0)
figure;
subplot(421);
plot(yy(:,1),yy(:,18),yyid(:,1),yyid(:,13),'--','linewidth',2);
xlabel('R');ylabel('n_e');
legend('boray','genray');legend('boxoff');
subplot(422);
plot(yy(:,1),yy(:,19),yyid(:,1),yyid(:,14),'--','linewidth',2);
xlabel('R');ylabel('T_e');
subplot(423);
plot(yy(:,1),yy(:,2),yyid(:,1),yyid(:,2),'--','linewidth',2);
xlabel('R');ylabel('Z');
subplot(424);
plot(yy(:,1),yy(:,4),yyid(:,1),yyid(:,4),'--','linewidth',2);
xlabel('R');ylabel('k_r');
subplot(425);
plot(yy(:,1),yy(:,5),yyid(:,1),yyid(:,5),'--','linewidth',2);
xlabel('R');ylabel('n_\phi');
subplot(426);
plot(yy(:,1),yy(:,6),yyid(:,1),yyid(:,6),'--','linewidth',2);
xlabel('R');ylabel('k_z');
subplot(427);
plot(yy(:,1),yy(:,9),yyid(:,1),yyid(:,15),'--','linewidth',2);
xlabel('R');ylabel('v_{gr}');
subplot(428);
% plot(yy(:,1),yy(:,11),yyid(:,1),yyid(:,17),'--','linewidth',2);
% xlabel('R');ylabel('v_{gz}');
plot(yy(:,1),yy(:,17),yyid(:,1),yyid(:,11),'--','linewidth',2);
xlabel('R');ylabel('B');
end

if(1==1)
if(icase==1)
save([savepath,'boray_f=',num2str(f/1e9,3),...
    'GHz,y0=',num2str(yy(1,1),3),',',num2str(yy(1,2),3),',',...
    num2str(yy(1,3),3),',',num2str(yy(1,4),3),',',num2str(yy(1,5),3),...
    ',',num2str(yy(1,6),3),',dt=',num2str(dt,2),',nt=',num2str(nt),...
    strtmp,'.mat'],'yyallray','yyray','yrayp','Ptp','yyp','f','ww','wwh','N_each');
elseif(icase==2 || icase==3)
save([savepath,'boray_f=',num2str(f/1e9,3),...
    'GHz,y0=',num2str(yy(1,1),3),',',num2str(yy(1,2),3),',',...
    num2str(yy(1,3),3),',',num2str(yy(1,4),3),',',num2str(yy(1,5),3),...
    ',',num2str(yy(1,6),3),',dt=',num2str(dt,2),',nt=',num2str(nt),...
    strtmp,'.mat'],'yyray','yrayp','Ptp','yyp','f','ww','wwh','N_each');
end
end
