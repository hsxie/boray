% 2021-05-05 15:13 Huasheng XIE, huashengxie@gmail.com, ENN
% boray_calpower.m

ntp0=200; % only  calculate ntp points
djp=floor(it/ntp0);
yyp=yy(1:djp:it,:);
% yyp=yy(1:2:100,:);

rp=yyp(:,1);
zp=yyp(:,3);
ns0=zeros(size(yyp,1),S);
Ts0=0.*ns0;
for s=1:S   
    ns0(:,s)=yyp(:,16+2*s);
    Ts0(:,s)=yyp(:,17+2*s);
end

B0=yyp(:,17);
kz=abs(yyp(:,15));
kx=abs(yyp(:,16));
tp=yyp(:,7);
ntp=length(tp);
% qs=qs.'; ms=ms.';

wt=f*2*pi+0.*tp;
%%
[ww,detD]=colddr(qs.',ms.',ns0,B0,kx,kz,wt);

% jeach=0; 
% joutw=1;
% N=3; 
% J=8;
warning('off');

% 21-05-01 22:33 calculate damping rate from each species
[ww_each,detD_each,N_each]=hotdreach(qs.',ms.',ns0,Ts0,B0,kx,kz,ww,jeach,joutw,N,J);
wwh=ww_each;
if(joutw==0) %  calculate omega_i using another method
    [ww_each2,detD_each2,N_each2]=hotdreach(qs.',ms.',ns0,Ts0,B0,kx,kz,...
        ww*0.99999,jeach,joutw,N,J);
    wi_each2=0.*ww_each2;
    for jp=1:length(ww)
        wi_each2(jp,:)=imag(detD_each(jp,:))./(real(detD_each2(jp,:)-...
            detD_each(jp,:))./(0.00001*ww(jp)));
    end
    wwh=real(ww_each)+1i*wi_each2;
end

dtp=tp(2)-tp(1);

Ptp=zeros(ntp,N_each);
for js=1:N_each
    damp=2*cumsum(imag(squeeze(wwh(:,js)))*dtp);
    Ptp(:,js)=1-exp(damp);  % asborb power along the ray for each species
end


% [tp.',rp.',zp.',wwh,Ptp]