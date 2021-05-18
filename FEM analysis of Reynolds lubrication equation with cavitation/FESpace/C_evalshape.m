function [dphiq,Grad] = C_evalshape(basis,node_2D)
%% [dphiq,Grad] = C_evalshape(basis,node_2D)
%==========================================================================
% Evaluation of basis functions and gradients on quadrature nodes 
%==========================================================================
%    called in C_matrix2D.m
%
%    INPUT:
%          basis    : (struct)  see C_shape_basis.m
%          node2D   : (array real) [nqn x 2] (x,y) coord of quadrature nodes   
%
%    OUTPUT:
%          dphi     : dphiq(:,:,i) vecotr containing the evaluation of the
%                     i-th basis function on the quadrature nodes 
%          Grad     : dphiq(:,1:2,i) vecotr containing the evaluation of the
%                     gradinet of dphi(:,:,i) on the quadrature nodes 


nln = length(basis);

%==========================================================================
% EVALUATION OF BASIS FUNCTIONS AND GRADIENTS
%==========================================================================

for s=1:nln
    % coordinates
    csi = node_2D(:,1);
    eta = node_2D(:,2);
    % evaluation of basis functions
    dphiq(1,:,s) = eval(basis(s).fbases);
    % evaluation of gradients
    Grad(:,1,s) = eval(basis(s).Gbases_1);
    Grad(:,2,s) = eval(basis(s).Gbases_2);
end
