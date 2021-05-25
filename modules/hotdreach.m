% Hua-sheng XIE, huashengxie@gmail.com, ENN, 2021-04-20 10:12
% kinetic Plasma Dispersion Relastion solver
% Maxwellian equilibrium distribution
% J-pole approximation for Z(zeta)=sum(b_j/(zeta-c_j))
% Rewrite for ray tracing
% 21-05-01 22:20 Calculate damping for each species
function [ww_each,detD_each,N_each]=hotdreach(qs,ms,fns0,fTs0,fB0,...
    fkx,fkz,fw,jeach,joutw,N,J)
% size:
%   =S. qs, ms
%   =npoint*S. ns0, Ts0
%   =npoint. B0, kx, kz, w

c2=(2.9979E8)^2; epsilon0=8.8542E-12;
% mu0=1/(c2*epsilon0);
kB=1.3807e-23;

% enlarge N to make sure the results are convergent
% global N;
% N=4; % Number harmonics
% J=12; % J-pole, usually J=8 is sufficient
if(J==8) % J-pole, usually J=8 is sufficient; other choice: 4, 12
  opt=3;
%   opt=2;
  if(opt==1)
    % Ronnmark1982, 8-pole for Z function, and Z', J=8, I=10
    bzj(1)=-1.734012457471826E-2-4.630639291680322E-2i;
    bzj(2)=-7.399169923225014E-1+8.395179978099844E-1i;
    bzj(3)=5.840628642184073+9.536009057643667E-1i;
    bzj(4)=-5.583371525286853-1.120854319126599E1i;
    czj(1)=2.237687789201900-1.625940856173727i;
    czj(2)=1.465234126106004-1.789620129162444i;
    czj(3)=.8392539817232638-1.891995045765206i;
    czj(4)=.2739362226285564-1.941786875844713i;
  elseif(opt==2) % J=8, I=8
    bzj(1)=  2.262756990456044 + 0.574863654552974i;
    bzj(2)=  0.007755112891899 + 0.435320280255620i;
    bzj(3)=  -0.017631167787210 - 0.001036701965463i;
    bzj(4)=  -2.752880935560734 - 4.080667522419399i;
    czj(1)=  -0.845779591068908 - 1.626106289905113i;
    czj(2)=  1.485551642694134 - 1.512332910786984i;
    czj(3)=  2.286671582218345 - 1.340185584073022i;
    czj(4)=  0.275217303800922 - 1.682797352333194i;
  else % new calculation, J=8, I=10
%     bzj(1)=  -0.0173401116032742 - 0.0463064419344598i;
%     bzj(2)=  -0.739917851897683 + 0.839518298070637i;
%     bzj(3)=  5.84063227513760 + 0.953602843950785i;
%     bzj(4)=  -5.58337431170864 - 11.2085508179677i;
%     czj(1)=   2.23768772215616 - 1.62594103256666i;
%     czj(2)=   1.46523409042510 - 1.78962030806222i;
%     czj(3)=   0.839253965702731 - 1.89199521968963i;
%     czj(4)=   0.273936217871668 - 1.94178704551807i;

    % 18-12-28 11:11
    bzj(1)=   -0.017340112270401 - 0.046306439626294i;
    bzj(2)=   -0.739917811220052 + 0.839518284620274i;
    bzj(3)=   5.840632105105495 + 0.953602751322040i;
    bzj(4)=   -5.583374181615043 -11.208550459628098i;
    czj(1)=   2.237687725134293 - 1.625941024120362i;
    czj(2)=   1.465234091939142 - 1.789620299603315i;
    czj(3)=   0.839253966367922 - 1.891995211531426i;
    czj(4)=   0.273936218055381 - 1.941787037576095i;
    
  end

  bzj(5:8)=conj(bzj(1:4));
  czj(5:8)=-conj(czj(1:4));
elseif(J==12) % from Cal_J_pole_bjcj.m
  opt=0;
  if(opt==1) % J=12; I=16; 2014;
   bzj(1)=  -0.00454786121654587 - 0.000621096230229454i;
   bzj(2)=    0.215155729087593 + 0.201505401672306i;
   bzj(3)=    0.439545042119629 + 4.16108468348292i;
   bzj(4)=  -20.2169673323552 - 12.8855035482440i;
   bzj(5)=    67.0814882450356 + 20.8463458499504i;
   bzj(6)=  -48.0146738250076 + 107.275614092570i;
  
   czj(1)=  -2.97842916245164 - 2.04969666644050i;
   czj(2)=    2.25678378396682 - 2.20861841189542i;
   czj(3)=  -1.67379985617161 - 2.32408519416336i;
   czj(4)=  -1.15903203380422 - 2.40673940954718i;
   czj(5)=    0.682287636603418 - 2.46036501461004i;
   czj(6)=  -0.225365375071350 - 2.48677941704753i;
  elseif(opt==2) % J=12; I=16; 2018-12-28 22:31
   bzj(1)= - 47.913598578418315281 - 106.98699311451399461i;
   bzj(2)= - 20.148858425809293248 + 12.874749056250453631i;
   bzj(3)= - 0.0045311004339957471789 + 0.00063311756354943215316i;
   bzj(4)= 0.2150040123642351701 + 0.20042340981056393122i;
   bzj(5)= 0.43131038679231352184 - 4.1505366661190555077i;
   bzj(6)= 66.920673705505055584 + 20.747375125403268524i;
   czj(1)= 0.2253670862838072698686 - 2.48625584284603285647646i;
   czj(2)= 1.1590491549279069691485 - 2.406192125704074076408834i;
   czj(3)= 2.9785703941315209703937799933362 - 2.0490809954949754985i;
   czj(4)= 2.25685878923092272938 - 2.208022912648570057162637i;
   czj(5)= 1.6738373878120108270775693017312 - 2.3235155478934783777i;
   czj(6)= 0.6822944098171246799968 - 2.459833442261711494628i;
  else % J=12; I=12; 18-12-28 22:35
   bzj(1)=  - 10.020983259474214017 - 14.728932929429874883i;
   bzj(2)= - 0.58878169153449514493 + 0.19067303610080007359i;
   bzj(3)= - 0.27475707659732384029 + 3.617920717493884482i;
   bzj(4)= 0.00045713742777499515344 + 0.00027155393843737098852i;
   bzj(5)=  0.017940627032508378515 - 0.036436053276701248142i;
   bzj(6)=  10.366124263145749629 - 2.5069048649816145967i;
   czj(1)= 0.22660012611958088507627 - 2.0716877594897791206264i;
   czj(2)= - 1.70029215163003500750575 - 1.8822474221612724460388i;
   czj(3)= 1.17139325085601178534269 - 1.97725033192085410977458i;
   czj(4)= 3.0666201126826972102007 - 1.59002082593259971758095i;
   czj(5)= 2.307327490410578276422 - 1.7546732543728200653674i;
   czj(6)= 0.687200524906019065672977 - 2.040288525975844018682i;
  end

  bzj(7:12)=conj(bzj(1:6));
  czj(7:12)=-conj(czj(1:6));
elseif(J==16)
  opt=1;
  if(opt==1) % J=16; I=18;  
   bzj(1)= - 86.416592794839804566 - 147.57960545984972964i;
   bzj(2)= - 22.962540986214500398 + 46.211318219085729914i;
   bzj(3)= - 8.8757833558787660662 - 11.561957978688249474i;
   bzj(4)= - 0.025134802434111256483 + 0.19730442150379382482i;
   bzj(5)= - 0.0056462830661756538039 - 0.0027884991898011769583i;
   bzj(6)= 0.000028262945845046458372 + 0.000026335348714810255537i;
   bzj(7)= 2.3290098166119338312 - 0.57238325918028725167i;
   bzj(8)= 115.45666014287557906 - 2.8617578808752183449i; 
   czj(1)= 0.1966439744113664608458045976392 - 2.5854046363167904820930238267552i;
   czj(2)= 1.000427687089304511157923736374 - 2.5277610669350594581215470753126i;
   czj(3)= 1.4263380087098663428834704281261 - 2.4694803409658086505344546718783i;
   czj(4)= 2.382753075769737513956410751299 - 2.2903917960623787648467068236658i;
   czj(5)= 2.9566517643704010426779572638885 - 2.1658992556376956216621262314559i;
   czj(6)= - 3.6699741330155866185481740497527 - 2.008727613312046260114172119472i;
   czj(7)= 1.8818356204685089975461092960437 - 2.3907395820644127767937911780402i;
   czj(8)= 0.5933003629474285223202174828712 - 2.5662607006180515205167080595386i;
  end
  bzj(9:16)=conj(bzj(1:8));
  czj(9:16)=-conj(czj(1:8));
elseif(J==24) 
  opt=1;
  if(opt==1)% J=24,I=24
    bzj(1)= - 579.77656932346560644 - 844.01436313629880827i;
    bzj(2)= - 179.52530851977905732 - 86.660002027244731382i;
    bzj(3)= - 52.107235029274485215 + 453.3246806707749413i;
    bzj(4)= - 2.1607927691932962178 + 0.63681255371973499384i;
    bzj(5)= - 0.018283386874895507814 - 0.21941582055233427677i;
    bzj(6)= - 0.00006819511737162705016 + 0.00032026091897256872621i;
    bzj(7)= - 0.0000028986123310445793648 - 0.00000099510625011385493369i;
    bzj(8)= 0.0000000023382228949223867744 - 0.0000000040404517369565098657i;
    bzj(9)= 0.01221466589423530596 + 0.00097890737323377354166i;
    bzj(10)= 7.3718296773233126912 - 12.575687057120635407i;
    bzj(11)= 44.078424019374375065 - 46.322124026599601416i;
    bzj(12)= 761.62579175738689742 + 185.11797721443392707i;
    czj(1)= 0.16167711630587375808393823760988 - 2.9424665391729649010502939606152i;
    czj(2)= 1.1509135876493567244599398043479 - 2.8745542965490153159866506667543i;
    czj(3)= 0.81513635269214329286824152984179 - 2.9085569383176322446978082849749i;
    czj(4)= 2.2362950589041724110736073820844 - 2.7033607074680388479084431872604i;
    czj(5)= 2.6403561313404041541230494846625 - 2.6228400297078984516779261304916i;
    czj(6)= 3.5620497451197056657834990483967 - 2.4245607245823420555878190731282i;
    czj(7)= 4.116925125710675393072860873751 - 2.3036541720854573608940600179944i;
    czj(8)= 4.8034117493360317933109830717707 - 2.1592490859689535412501218722927i;
    czj(9)= 3.0778922349246567316482750461458 - 2.5301774598854448463007864644617i;
    czj(10)= - 1.8572088635240765003561090479193 - 2.7720571884094886583775397071469i;
    czj(11)= 1.4969881322466893380396663902149 - 2.8290855580900544693059801078858i;
    czj(12)= - 0.48636891219330428093331493852099 - 2.9311741817223824196339069754696i;
  end
  bzj(13:24)=conj(bzj(1:12));
  czj(13:24)=-conj(czj(1:12));
  
elseif(J==4)
  opt=2;
  if(opt==1)
  % Martin1980, 4-pole, J=4, I=5
    bzj(1)=0.5468-0.0372i;
    bzj(2)=-1.0468+2.1018i;
    czj(1)=-1.2359-1.2150i;
    czj(2)=-0.3786-1.3509i;
  else % new calculation, J=4, I=5
    bzj(1)=0.546796859834032 + 0.037196505239277i;
    bzj(2)=-1.046796859834027 + 2.101852568038518i;
    czj(1)=1.23588765343592 - 1.21498213255731i;
    czj(2)=-0.378611612386277 - 1.350943585432730i;
  end
  
  bzj(3:4)=conj(bzj(1:2));
  czj(3:4)=-conj(czj(1:2));
elseif(J==3)
  % Martin1980, 3-pole, J=3, I=3
  bzj(1)=0.1822+0.5756i;
  bzj(2)=-1.3643;
  czj(1)=-0.9217-0.9091i;
  czj(2)=-1.0204i;
  
  bzj(3)=conj(bzj(1));
  czj(3)=-conj(czj(1));
elseif(J==2)
  % Huba2009, 2-pole
  bzj(1)=-(0.5+0.81i);
  czj(1)=0.51-0.81i;
  bzj(2)=conj(bzj(1));
  czj(2)=-conj(czj(1));
else
  % Martin1979, 2-pole, J=2, I=3
  bzj(1)=-(0.5+1.2891i);
  czj(1)=0.5138-1.0324i;
  bzj(2)=conj(bzj(1));
  czj(2)=-conj(czj(1));
end

% [S,npoint]=size(fns0);
[npoint,S]=size(fns0);

if(jeach==0)
    ww_each=zeros(npoint,1);
    detD_each=zeros(npoint,1);
    N_each=1;
elseif(jeach==1)
    ww_each=zeros(npoint,1+S);
    detD_each=zeros(npoint,1+S);
    N_each=1+S;
end

for jp=1:npoint
    
    ns0=squeeze(fns0(jp,:));
    Ts0=squeeze(fTs0(jp,:));
    B0=fB0(jp);
    kx=fkx(jp);
    kz=fkz(jp);
    w=fw(jp);
for js=1:N_each % to calculate the omega_i from each species, with others be cold
    if(jp>1) % 21-04-27 11:02
%         w=ww(jp-1);
        w=ww_each(jp-1,js);
    end
%     Ts=Ts0*1.6022e-19/kB; % Ts, eV -> K
    if(js==1)
        Ts=Ts0*1.6022e-19/kB; % Ts, eV -> K
    else %if(js>1)
        Ts=0.*Ts0+1.0*1.6022e-19/kB; % default 1.0eV
        Ts(js-1)=Ts0(js-1)*1.6022e-19/kB; % Ts, eV -> K
    end
%     Ts
    
    vtzs=sqrt(2*kB*Ts./ms); % para thermal velocity, note the sqrt(2)
    wps=sqrt(ns0.*qs.^2./ms/epsilon0); % plasma frequency
    wcs=B0*qs./ms; % cyclotron frequency
    rhocs=sqrt(kB*Ts./ms)./wcs; % cyclotron radius, 2018-06-13 21:47
    
    wps2=wps.^2;
    
    SNJ=S*(2*N+1)*J;
%     SNJ1=SNJ+1;
    SNJ1=SNJ+0;
    SNJ3=3*SNJ1;
    NN=SNJ3+6;

    % % -- main program begin to set the matrix elements --
%     k=sqrt(kz^2+kx^2);
    bs=kx*rhocs;
    bs(abs(bs)<1e-50)=1e-50;  % to avoid singular when k_perp=0

    bs2=bs.^2;
    M=sparse(NN,NN);
    snj=0;

    % initialize
    csnj=zeros(1,3*SNJ);
    b11snj=csnj.*0; b12snj=csnj.*0; b13snj=csnj.*0;
    b21snj=csnj.*0; b22snj=csnj.*0; b23snj=csnj.*0;
    b31snj=csnj.*0; b32snj=csnj.*0; b33snj=csnj.*0;

    for s=1:S % species
      for n=-N:N % Bessel function
        Gamn=besseli(n,bs2(s),1); % 2014-10-13 12:51
        Gamnp=(besseli(n+1,bs2(s),1)+...
          besseli(n-1,bs2(s),1)-2*besseli(n,bs2(s),1))/2;
        for j=1:J % poles of Z(zeta)
          snj=snj+1;

          csnj(snj)=czj(j)*kz*vtzs(s)+n*wcs(s);

          tmp=wps2(s)*bzj(j);

          b11snj(snj)=tmp*n^2*Gamn/bs2(s);

          b12snj(snj)=tmp*1i*n*Gamnp;
          b21snj(snj)=-b12snj(snj);

          b22snj(snj)=tmp*(n^2*Gamn/bs2(s)-2*bs2(s)*Gamnp);

          %
          bnj2=czj(j);

          b13snj(snj)=tmp*bnj2*n*sqrt(2)*Gamn/bs(s);
          b31snj(snj)=b13snj(snj);

          b23snj(snj)=-1i*tmp*bnj2*sqrt(2)*Gamnp*bs(s);
          b32snj(snj)=-b23snj(snj);

          bnj2=czj(j)*czj(j);

          b33snj(snj)=tmp*bnj2*2*Gamn;
        end
      end
    end

    if(joutw==1)
    for snj=1:SNJ % set the eigen matrix
      jjx=snj+0*SNJ1;
      jjy=snj+1*SNJ1;
      jjz=snj+2*SNJ1;
      % v_snjx
      M=M+sparse(jjx,jjx,csnj(snj),NN,NN)+...
        sparse(jjx,SNJ3+1,b11snj(snj),NN,NN)+...
        sparse(jjx,SNJ3+2,b12snj(snj),NN,NN)+...
        sparse(jjx,SNJ3+3,b13snj(snj),NN,NN);

      % v_snjy
      M=M+sparse(jjy,jjy,csnj(snj),NN,NN)+...
        sparse(jjy,SNJ3+1,b21snj(snj),NN,NN)+...
        sparse(jjy,SNJ3+2,b22snj(snj),NN,NN)+...
        sparse(jjy,SNJ3+3,b23snj(snj),NN,NN);

      % v_snjz
      M=M+sparse(jjz,jjz,csnj(snj),NN,NN)+...
        sparse(jjz,SNJ3+1,b31snj(snj),NN,NN)+...
        sparse(jjz,SNJ3+2,b32snj(snj),NN,NN)+...
        sparse(jjz,SNJ3+3,b33snj(snj),NN,NN);

    end

    % E(J), J_{x,y,z}=sum(v_snj{x,y,z})
    tp=-1;
    jj=(0*SNJ1+1):(1*SNJ1); ii=0.*jj+SNJ3+1; M=M+sparse(ii,jj,tp,NN,NN);
    jj=(1*SNJ1+1):(2*SNJ1); ii=0.*jj+SNJ3+2; M=M+sparse(ii,jj,tp,NN,NN);
    jj=(2*SNJ1+1):(3*SNJ1); ii=0.*jj+SNJ3+3; M=M+sparse(ii,jj,tp,NN,NN);

    % E(B)
    M=M+sparse(SNJ3+1,SNJ3+5,c2*kz,NN,NN)+...
      sparse(SNJ3+2,SNJ3+4,-c2*kz,NN,NN)+...
      sparse(SNJ3+2,SNJ3+6,c2*kx,NN,NN)+...
      sparse(SNJ3+3,SNJ3+5,-c2*kx,NN,NN);

    % B(E)
    M=M+sparse(SNJ3+4,SNJ3+2,-kz,NN,NN)+...
      sparse(SNJ3+5,SNJ3+1,kz,NN,NN)+...
      sparse(SNJ3+5,SNJ3+3,-kx,NN,NN)+...
      sparse(SNJ3+6,SNJ3+2,kx,NN,NN);
  
%     d0=vpa(eigs(M,1,w),16);
%     ww(jp)=double(d0);
%     options=optimoptions('eigs','Display','off');
    wtmp=eigs(M,1,w);
    ww_each(jp,js)=wtmp(1);
    elseif(joutw==0) % do not calcualte omega from kinetic DR
        ww_each(jp,js)=w;
    end
    
    % % calculate tensor
    wd=w;
    sigmawk=zeros(3,3);
    sigmawk(1,1)=sum(b11snj./(wd-csnj));
    sigmawk(1,2)=sum(b12snj./(wd-csnj));
    sigmawk(1,3)=sum(b13snj./(wd-csnj));
    sigmawk(2,1)=sum(b21snj./(wd-csnj));
    sigmawk(2,2)=sum(b22snj./(wd-csnj));
    sigmawk(2,3)=sum(b23snj./(wd-csnj));
    sigmawk(3,1)=sum(b31snj./(wd-csnj));
    sigmawk(3,2)=sum(b32snj./(wd-csnj));
    sigmawk(3,3)=sum(b33snj./(wd-csnj));
    sigmawk=-1i*epsilon0*sigmawk; % conductivity tensor of each species s
    
    Qwk=-sigmawk/(1i*wd*epsilon0);
    Kwk=eye(3)+Qwk; % dielectric tensor
    Dwk=Kwk+(kron([kx,0,kz],[kx;0;kz])/(kx^2+kz^2)-eye(3))*((kx^2+...
        kz^2)*c2/wd^2); % dispersion tensor, D(w,k){\cdot}E=0
            
    detD_each(jp,js)=det(Dwk);
end
end
end
