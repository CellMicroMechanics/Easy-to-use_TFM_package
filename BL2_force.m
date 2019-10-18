% Copyright (C) 2010 - 2019, Sabass Lab
%
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the Free
% Software Foundation, either version 3 of the License, or (at your option) 
% any later version. This program is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
% Public License for more details. You should have received a copy of the 
% GNU General Public License along with this program.
% If not, see <http://www.gnu.org/licenses/>.

function [ Ftfx, Ftfy]= BL2_force(Ftux,Ftuy,L,E,s,cluster_size,i_max, j_max) 

    V = 2*(1+s)/E;
    kx_vec = 2*pi/i_max/cluster_size.*[0:(i_max/2-1) (-i_max/2:-1)];
    ky_vec = 2*pi/j_max/cluster_size.*[0:(j_max/2-1) (-j_max/2:-1)];
    kx = repmat(kx_vec',1,j_max);
    ky = repmat(ky_vec,i_max,1);
    
    kx(1,1) = 1;
    ky(1,1) = 1;

    Ginv_xx = (kx.^2+ky.^2).^(-1/2).*V.*(kx.^2.*L+ky.^2.*L+V.^2).^(-1).*(kx.^2.* ...
               L+ky.^2.*L+((-1)+s).^2.*V.^2).^(-1).*(kx.^4.*(L+(-1).*L.*s)+ ...
               kx.^2.*((-1).*ky.^2.*L.*((-2)+s)+(-1).*((-1)+s).*V.^2)+ky.^2.*( ...
               ky.^2.*L+((-1)+s).^2.*V.^2));
    Ginv_yy =  (kx.^2+ky.^2).^(-1/2).*V.*(kx.^2.*L+ky.^2.*L+V.^2).^(-1).*(kx.^2.* ...
               L+ky.^2.*L+((-1)+s).^2.*V.^2).^(-1).*(kx.^4.*L+(-1).*ky.^2.*((-1)+ ...
               s).*(ky.^2.*L+V.^2)+kx.^2.*((-1).*ky.^2.*L.*((-2)+s)+((-1)+s).^2.* ...
               V.^2));
    Ginv_xy =  (-1).*kx.*ky.*(kx.^2+ky.^2).^(-1/2).*s.*V.*(kx.^2.*L+ky.^2.*L+ ...
               V.^2).^(-1).*(kx.^2.*L+ky.^2.*L+((-1)+s).*V.^2).*(kx.^2.*L+ky.^2.* ...
               L+((-1)+s).^2.*V.^2).^(-1);

    Ginv_xx(1,1) = 0;
    Ginv_yy(1,1) = 0;
    Ginv_xy(1,1) = 0;

    Ginv_xy(i_max/2+1,:) = 0;
    Ginv_xy(:,j_max/2+1) = 0;

    Ftfx = Ginv_xx.*Ftux + Ginv_xy.*Ftuy;
    Ftfy = Ginv_xy.*Ftux+ Ginv_yy.*Ftuy;
 
end