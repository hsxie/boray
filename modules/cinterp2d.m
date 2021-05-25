% Hua-sheng XIE, huashengxie@gmail.com, ENN, 2021-05-04 21:20 
% Bi-linear interp2d, (r,z) uniform grid
% Calculate the coefficients fc(1:4,:)
% jr=floor((ra-r(1))/dr)+1; jz=floor((za-z(1))/dz)+1;
% hr=ra-r(jr); hz=za-z(jz);
% ya=fc(1,jr,jz)+fc(2,jr,jz)*hr+fc(3,jr,jz)*hz+fc(4,jr,jz)*hz*hr
function fc=cinterp2d(rg,zg,frz)

[nr,nz]=size(frz);
dr=rg(2)-rg(1);
dz=zg(2)-zg(1);

fc=zeros(4,nr,nz);
fc(1,:,:)=frz;
fc(2,1:(nr-1),:)=(fc(1,2:nr,:)-fc(1,1:(nr-1),:))/dr;
fc(3,:,1:(nz-1))=(fc(1,:,2:nz)-fc(1,:,1:(nz-1)))/dz;

fc(4,1:(nr-1),1:(nz-1))=(fc(1,2:nr,2:nz)-fc(1,1:(nr-1),1:(nz-1))-...
    fc(2,1:(nr-1),1:(nz-1))*dr-fc(3,1:(nr-1),1:(nz-1))*dz)/(dr*dz);

end
