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
%DESCRIPTION
%Function for calculating regularization parameter using Bayesian method

%------------------
%FUNCTION ARGUMENTS 
%beta: 1/variance of noise 
%fuu: displacement vector in Fourior space
%Ftux: x component of displacement matrix in Fourior space
%Ftuy: y component of displacement matrix in Fourior space
%E: Young's modulus
%s: Poisson's ratio
%cluster_size: grid spacing in pixels
%grid_mat: regular grid with size i_max*j_max 
%u: displacement vectors on grid
%i_max, j_max: sizes of grid
%X: matrix between displacement and force in Fourior space
%sequence: set to 1 if only maximum evidence parameter should be returned
%------------------

%------------------
%FUNCTION OUTPUTS
%lambda_2: optimal regularization parameter
%evidencep: matrix for regularization parameter and its value of
%          logevidence around the optimal regularization parameter 
%evidence_one: value of logevidence at the optimal regularization parameter
%------------------


function [lambda_2 evidencep evidence_one]  = optimal_lambda(beta,fuu,Ftux,Ftuy,E,s,cluster_size,i_max, j_max,X,sequence)

aa = size(X); 
c = ones(aa(2),1);
C = spdiags(c, 0:0,aa(2),aa(2));
XX = sparse(X)'*sparse(X);
BX_a = beta*sparse(XX)/aa(1)*2;
C_a = C/aa(2)*2;
constant = aa(1)*log(beta)-aa(1)*log(2*pi);


%%% Golden section search method to find alpha at minimum of -log(Evidence)
%%% 
%setting the range of parameter search. Change if maximum can not be found in your data
alpha1 =1e-8; 
alpha2 =1e8; 

%search optimal parameter
alpha_opt = fminbnd(@minus_logevidence,alpha1,alpha2);

    if nargin ==10 || ~sequence   %%%produce data for plotting the evidence function
           plot_alpha = (alpha_opt*0.2:alpha_opt*0.12:alpha_opt*2);
           a = size(plot_alpha);
           lambda_p = plot_alpha./beta;
           for i = 1:a(2)
             evidence(i) = -minus_logevidence(plot_alpha(i));
           end
           evidencep = [lambda_p;evidence];
    else
           evidencep = 0; %only calculate the maximum of the evidence curve
    end

  evidence_one = -minus_logevidence(alpha_opt);
  lambda_2 = alpha_opt/beta;
 
 
 
%%%Nested function for calculating -log(Evidence)  
      function evidence_value= minus_logevidence(alpha)
        LL = alpha/beta;
        [~,~,~,~,Ftfx, Ftfy] = reg_fourier_TFM(Ftux,Ftuy,LL,E,s,cluster_size,i_max, j_max);
        fxx = reshape(Ftfx,i_max*j_max,1);    
        fyy = reshape(Ftfy,i_max*j_max,1);
        f(1:2:size(fxx)*2,1) = fxx;
        f(2:2:size(fyy)*2,1) = fyy;

         A = alpha*sparse(C_a) + BX_a;
         L = chol(sparse(A));
         logdetA = 2*sum(log(diag(L)));

         Xf_u = X*f-fuu;
         Ftux1= Xf_u(1:2:end);
         Ftuy1= Xf_u(2:2:end);
         ff = sum(sum(Ftfx.*conj(Ftfx) + Ftfy.*conj(Ftfy)))/(0.5*aa(2));
         uu = sum(sum(Ftux1.*conj(Ftux1) + Ftuy1.*conj(Ftuy1)))/(0.5*aa(1));

         evidence_value = -0.5*(-alpha*ff-beta*uu ...
                          -logdetA +aa(2)*log(alpha)+constant);
      end
 
end


