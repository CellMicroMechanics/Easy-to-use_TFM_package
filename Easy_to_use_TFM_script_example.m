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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%% An example showing how to use Easy-to-use_TFM from the command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%clear the variables that we will use
clear input_data E s pix_durch_mu meshsize frame input_pos input_vec ...
    input_noise pos traction traction_magnitude f_n_m fnorm grid_mat i_max ...
    j_max X fuu Ftux Ftuy u Ftu regparam L

%load the test data
load('test_data/input_data.mat');

%some variables
frame =1;               %frame number in movie data sequence
meshsize = 10;          %grid spacing in pix
pix_durch_mu = 0.1;     %size of one pixel in micrometers 

E = 10000;      %Young's modulus in Pa
s = 0.3;       %Poisson's ratio

%extract the input data we need
input_pos   = input_data.displacement(frame).pos;
input_vec   = input_data.displacement(frame).vec;
input_noise = input_data.noise(frame).vec;


%% Calculate traction with Bayesian TFM
%NOTE: Bayesian TFM is a method to reconstruct traction forces in a
%well-defined manner without requiring low-pass data filtering or
%regularization of the solution. For Bayesian TFM, one requires an estimate 
%of the variance of the measurement noise. The noise variance can be either estimated 
%with knowledge of the experimental conditions or one can provide a displacement sample that
%contains only noise. Such a sample could, e.g., be obtained by tracking
%stationary fiducial markers in a TFM setup.

%use the provided noise sample to calculate its variance [pix^2] 
varnoise = var(input_noise(:));
beta = 1/varnoise;

%calculate with E=1 and rescale Young's modulus at the end
%prepare Fourier-tranformed Green's function and displacement
[grid_mat,i_max,j_max, X,fuu,Ftux,Ftuy,u] = fourier_X_u(input_pos,input_vec, meshsize, 1, s,[]);
            
%Calculate the optimal regularization parameter
[L,~,~] = optimal_lambda(beta,fuu,Ftux,Ftuy,1,s,meshsize,i_max, j_max,X,1);

%calculate traction forces with optimal regularization parameter
[pos,traction,traction_magnitude,f_n_m,~,~] = reg_fourier_TFM(Ftux,Ftuy,L,1,s, meshsize,...
    i_max,j_max, grid_mat,pix_durch_mu, 0);

%rescale traction with proper Young's modulus
traction = E*traction;
traction_magnitude = E*traction_magnitude;
f_n_m = E*f_n_m;

%display a heatmap of the spatial traction distribution
fnorm = (f_n_m(:,:,2).^2 + f_n_m(:,:,1).^2).^0.5;
figure; colormap jet; hold on;
surf(grid_mat(:,:,1), grid_mat(:,:,2),fnorm),view(0,90),shading interp, axis equal;
set(gca, 'DataAspectRatio', [1,1,50],'YDir','reverse');
colorbar;
title('Traction magnitude calculated with Bayesian TFM');
hold off;

%% Calculate traction with Regularized TFM
%NOTE: Guessing the correct regularization parameter can be difficult.
%Nevertheless, there may be instances where one wants to use Regularized
%TFM with a given regularization parameter.

%provide a regularization parameter
regparam=140;

%interpolate data onto a regular rectangular grid
[grid_mat,u, i_max,j_max] = interp_vec2grid(input_pos,input_vec,meshsize);
%FFT of displacement field
Ftu(:,:,1) = fft2(u(:,:,1));
Ftu(:,:,2) = fft2(u(:,:,2));

%calculate traction forces
[pos,traction,traction_magnitude,f_n_m,~,~] = reg_fourier_TFM(Ftu(:,:,1), Ftu(:,:,2), regparam, 1, s,...
                  meshsize, i_max, j_max, grid_mat, pix_durch_mu, 0);                  

%rescale Young's modulus
traction = E*traction;
traction_magnitude = E*traction_magnitude;
f_n_m = E*f_n_m;

%display a heatmap of the spatial traction distribution
fnorm = (f_n_m(:,:,2).^2 + f_n_m(:,:,1).^2).^0.5;
figure; colormap jet; hold on;
surf(grid_mat(:,:,1), grid_mat(:,:,2),fnorm),view(0,90),shading interp, axis equal;
set(gca, 'DataAspectRatio', [1,1,50],'YDir','reverse');
colorbar;
title('Traction magnitude calculated with Regularized TFM');
hold off;
