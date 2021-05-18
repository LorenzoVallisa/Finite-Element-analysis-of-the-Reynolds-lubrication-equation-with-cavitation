function [K_loc]=C_lap_loc(Grad,w_2D,nln,BJ,pphys_2D,nu_fun)
%% [K_loc]=C_lap_loc(Grad,w_2D,nln,BJ)
%==========================================================================
% Build the local stiffness matrix for the term grad(u)grad(v)
%==========================================================================
%    called in C_matrix2D.m
%
%    INPUT:
%          Grad        : (array real) evaluation of the gradient on
%                        quadrature nodes
%          w_2D        : (array real) quadrature weights
%          nln         : (integer) number of local unknowns
%          BJ          : (array real) Jacobian of the map 
%
%    OUTPUT:
%          K_loc       :  (array real) Local stiffness matrix


K_loc=zeros(nln,nln);
% nu = @(x) (12e-6*(1+0.4*cos(x))/0.000017).^3;
% nu_vect = nu([0.5 0.5 0]);
x = pphys_2D(:,1);
y =pphys_2D(:,2);
% nu_vect = nu(x);
nu_vect = eval(nu_fun);

for i=1:nln
    for j=1:nln
        for k=1:length(w_2D)
            Binv = inv(BJ(:,:,k));   % inverse
            Jdet = det(BJ(:,:,k));   % determinant 
            K_loc(i,j) = K_loc(i,j) + (Jdet.*w_2D(k)) .* nu_vect(k).*( (Grad(k,:,i) * Binv) * (Grad(k,:,j) * Binv )');
        end
    end
end



                                              
                                              

