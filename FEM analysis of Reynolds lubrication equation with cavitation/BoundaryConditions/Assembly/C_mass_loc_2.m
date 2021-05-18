function [M_loc]=C_mass_loc_2(dphiq,w_2D,nln,BJ,pphys_2D,nu_fun,obstacle_fun)
%% [M_loc]=C_mass_loc(dphiq,w_2D,nln,BJ)
%==========================================================================
% Build the local mass matrix for the term (uv)
%==========================================================================
%    called in C_matrix2D.m
%
%    INPUT:
%          dphiq       : (array real) evaluation of the basis function on
%                        quadrature nodes
%          w_2D        : (array real) quadrature weights
%          nln         : (integer) number of local unknowns
%          BJ          : (array real) Jacobian of the map 
%
%    OUTPUT:
%          M_loc       :  (array real) Local mass matrix

    M_loc=zeros(nln,nln);
    K_loc=zeros(nln,nln);

    
    x = pphys_2D(:,1);
    y = pphys_2D(:,2);
    nu_vect = eval(nu_fun);
    cavitation_pressure = eval(obstacle_fun);

    
    for i=1:nln
        for j=1:nln
            for k=1:length(w_2D)
                Binv = inv(BJ(:,:,k));      % inverse
                Jdet = det(BJ(:,:,k));      % determinant
                M_loc(i,j) = M_loc(i,j) + (Jdet.*w_2D(k))  .* nu_vect(k) .* (dphiq(1,k,i)).* (dphiq(1,k,j)).*(cavitation_pressure(k,:));
            end
        end
    end
    
end



    


                                              
                                              

