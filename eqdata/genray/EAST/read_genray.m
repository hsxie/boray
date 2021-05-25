% 21-05-06 22:07 Huasheng Xie, huashengxie@gmail.com, ENN
% read genray netcdf file

ncfile='./EAST_compare_lowtemp/east_X/east_den5e19_temp500eV_phi20.nc';
% ncfile='./EAST_compare_lowtemp/east_O/east_den5e19_temp500eV_phi20.nc';
% ncfile='./lhw/east_lh_multiray.nc';

wr                = ncread(ncfile,'wr')/100;             % r (m)
wz                = ncread(ncfile,'wz')/100;             % z (m)
wphi              = ncread(ncfile,'wphi');               % toroidal angle phi (rad)
wn_r              = ncread(ncfile,'wn_r');               % N_r refractive index component
wn_z              = ncread(ncfile,'wn_z');               % N_z refractive index component
wn_phi            = ncread(ncfile,'wn_phi');             % N_phi refractive index component
w_theta_pol       = ncread(ncfile,'w_theta_pol')*pi/180; % poloidal angle theta (rad)
wnpar       = ncread(ncfile,'wnpar'); % parallel refractive index
wnper       = ncread(ncfile,'wnper'); % perpendicular refractive index
delpwr            = ncread(ncfile,'delpwr')*1e-7;        % power (W)
te                = ncread(ncfile,'ste')*1e3;                % temperature (eV)
ne                = ncread(ncfile,'sene')*1e6;           % plasmas density (m^-3)
br                = ncread(ncfile,'sb_r')*1e-4;          % B_r (T)
bz                = ncread(ncfile,'sb_z')*1e-4;          % B_z (T)
bphi              = ncread(ncfile,'sb_phi')*1e-4;        % B_phi (T)
btot              = ncread(ncfile,'sbtot')*1e-4;        % B_tot (T)
c=2.99792458e8; %  m/s
vgr_r             = ncread(ncfile,'vgr_r')*c;        % vgroup_r (m/s)
vgr_z             = ncread(ncfile,'vgr_z')*c;        % vgroup_z (m/s)
vgr_phi           = ncread(ncfile,'vgr_phi')*c;        % vgroup_phi (m/s)

flux_r            = ncread(ncfile,'flux_r')*c;     % flux_r (m/s)
flux_z            = ncread(ncfile,'flux_z')*c;     % flux_z (m/s)
flux_phi          = ncread(ncfile,'flux_phi')*c;     % flux_phi (m/s)

freqcy          = ncread(ncfile,'freqcy');  % Wave frequency, Hz

dmas          = ncread(ncfile,'dmas');  % plasma species mass: electrons, then ions. Normalized to electron mass
charge           = ncread(ncfile,'charge');  % plasma species charge: electrons, then ions. Normalized to electronic charge
eqdsk_r          = ncread(ncfile,'eqdsk_r');  % eqdsk r array, m
eqdsk_z          = ncread(ncfile,'eqdsk_z');  % eqdsk z array, m
eqdsk_psi          = ncread(ncfile,'eqdsk_psi');  % eqdsk psi array

indexrho          = ncread(ncfile,'indexrho');  % Radial coord type
psifactr          = ncread(ncfile,'psifactr');  % Reduces Psi-Value of LCFS
binvol          = ncread(ncfile,'binvol');  % Volumes of radial bins, cm^3
binarea          = ncread(ncfile,'binarea');  % Areas of radial bins, cm^2
rho_bin          = ncread(ncfile,'rho_bin');  % normalized small radius bin boundaries
rho_bin_center          = ncread(ncfile,'rho_bin_center');  % normalized small radius bin centers
densprof          = ncread(ncfile,'densprof')*1e6;  % plasma density at bin boundaries, e and ions. particles/m^3
temprof          = ncread(ncfile,'temprof')*1e3;  % plasma temperatures at bin boundaries, e and ions. eV
zefprof           = ncread(ncfile,'zefprof');  % plasma Zeff at bin boundaries, e and ions

w_r_densprof_nc          = ncread(ncfile,'w_r_densprof_nc');  % major radius r mesh for plasma profiles. m
w_dens_vs_r_nc          = ncread(ncfile,'w_dens_vs_r_nc');  % plasma density at bin bndries, e and ions vs major radius. particles/m^3
w_temp_vs_r_nc          = ncread(ncfile,'w_temp_vs_r_nc')*1e3;  % plasma temperatures at bin bndries, e and ions vs major radius. eV
w_zeff_vs_r_nc          = ncread(ncfile,'w_zeff_vs_r_nc');  % plasma Zeff at bin bndries vs major radius

iabsorp          = ncread(ncfile,'iabsorp');

f=freqcy;
w=2*pi*f;
% c=2.99792458e8;

wkr=wn_r/c*w;
wkz=wn_z/c*w;
wkphi=wn_phi/c*w;
