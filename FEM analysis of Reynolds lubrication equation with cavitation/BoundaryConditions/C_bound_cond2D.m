%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Imposizione delle condizioni al contorno per il problema di
% Poisson non omogeneo omogeneo.
% INPUT: A matrice di stiffness
%        b vettore termine noto
%        femregion struttura contenente la matrice metrica
% OUTPUT: A,b modificate secondo le condizioni al bordo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,b,u_g,M]=C_bound_cond2D(A,b,M,femregion,Dati,bc,TestName)
%% [A,b,u_g,M]=C_bound_cond2D(A,b,M,femregion,Dati)
%==========================================================================
% Assign Dirchlet boundary conditions
%==========================================================================
%    called in C_main2D.m
%
%    INPUT:
%          A           : (sparse(ndof,ndof) real) stiffness matrix
%          b           : (sparse(ndof,1) real) rhs vector
%          femregion   : (struct)  see C_create_femregion.m
%          Dati        : (struct)  see C_dati.m

%
%    OUTPUT:
%          A           : (sparse(ndof,ndof) real) stiffness matrix
%          b           : (sparse(ndof,1) real) rhs vector
%          u_g         : (sparse(ndof,1) real) evaluation of Dirichlet conditions
%          M           : sparse(ndof,ndof) real mass matrix
%

fprintf('============================================================\n')
fprintf('Assign Dirichlet boundary conditions ... \n');
fprintf('============================================================\n')
if strcmpi('non_sym',bc)
    
    % Dirichlet Boundary conditions
    
    
    bp = femregion.boundary_points;
    
    
    
    
    for k = 1:length(bp)
        A(bp(k),:) = 0;
        A(bp(k),bp(k)) = 1;
        b(bp(k)) = 0;
    end
    
    
    
    if strcmpi('Test1',TestName)
        
        % Periodic Boundary conditions
        
        per_dx = femregion.per_dx;
        per_sx = femregion.per_sx;
        
        for jj = 1:length(per_dx)
            A(per_sx(jj),:) = A(per_sx(jj),:) + A(per_dx(jj),:);
            b(per_sx(jj)) =  b(per_sx(jj)) +  b(per_dx(jj));
        end
        
        
        
        for jj = 1:length(per_dx)
            A(per_dx(jj),:) = 0;
            A(per_dx(jj),per_dx(jj)) = 1;
            A(per_dx(jj),per_sx(jj)) = -1;
            b(per_dx(jj)) = 0;
        end
        
    end
    
    ndof = length(b);
    u_g = sparse(ndof,1);
    
    
else
    
    
    
    ndof = length(b);
    u_g = sparse(ndof,1);
    
    x = femregion.dof(boundary_points,1);
    y = femregion.dof(boundary_points,2);
    u_g(boundary_points) = eval(Dati.exact_sol);
    
    if max(abs(bc))~=0
        u_g(period_points) = bc;
    end
    x_g = sparse(ndof,1);
    A_0 = A;
    M_0 = M;
    
    b_0 = b-A*u_g;
    
    
    bp = boundary_points;
    
    
    for k = 1:length(bp)
        A_0(bp(k),:) = 0;
        A_0(:,bp(k)) = 0;
        A_0(bp(k),bp(k)) = 1;
        M_0(bp(k),:) = 0;
        M_0(:,bp(k)) = 0;
        M_0(bp(k),bp(k)) = 1;
        b_0(bp(k)) = 0;
    end
    
    b = b_0;
    A = A_0;
    M = M_0;
end

