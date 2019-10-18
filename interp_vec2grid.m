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

function [grid_mat,u, i_max, j_max] = interp_vec2grid(pos, vec, cluster_size, grid_mat)
    if nargin == 3 || isempty(grid_mat)
        max_eck(1:2) = [max(pos(:,1)), max(pos(:,2))];
        min_eck(1:2) = [min(pos(:,1)), min(pos(:,2))];

        i_max = floor((max_eck(1)-min_eck(1))/cluster_size);
        j_max = floor((max_eck(2)-min_eck(2))/cluster_size);
        i_max = i_max - mod(i_max,2);
        j_max = j_max - mod(j_max,2);
        
        [X,Y] = meshgrid(min_eck(1)+(1/2:1:(i_max))*cluster_size, min_eck(2)+(1/2:1:(j_max))*cluster_size);
        grid_mat(:,:,1) = X';
        grid_mat(:,:,2) = Y';
        clear X Y;
      
    else
        i_max = size(grid_mat,1);
        j_max = size(grid_mat,2);
        cluster_size = grid_mat(1,1,1) - grid_mat(2,2,1);
    end
    
    if any(isnan(vec))
        disp('Warning: original data contains NAN values. Removing these values!');
        pos(isnan(vec(:,1)) | isnan(vec(:,2)),:) = [];
        vec(isnan(vec(:,1)) | isnan(vec(:,2)),:) = [];
    end
    
    u(1:i_max,1:j_max,1) = griddata(pos(:,1),pos(:,2),vec(:,1),grid_mat(:,:,1),grid_mat(:,:,2),'cubic');
    u(1:i_max,1:j_max,2) = griddata(pos(:,1),pos(:,2),vec(:,2),grid_mat(:,:,1),grid_mat(:,:,2),'cubic');
    u(isnan(u)) = 0;    
end