function [A,f,M]=C_matrix2D(Dati,femregion,dummy,xx,u)
%% [A,f,M] = C_matrix2D(Dati,femregion)
%==========================================================================
% Assembly of the stiffness matrix A, rhs f and mass matrix M
%==========================================================================
%    called in C_main2D.m
%
%    INPUT:
%          Dati        : (struct)  see C_dati.m
%          femregion   : (struct)  see C_create_femregion.m
%
%    OUTPUT:
%          A           : (sparse(ndof,ndof) real) stiffnes matrix
%          f           : (sparse(ndof,1) real) rhs vector
%          M           : (sparse(ndof,ndof) real) mass matrix
%

warning off;
addpath FESpace
addpath Assembly
addpath BoundaryConditions

fprintf('============================================================\n')
fprintf('Assembling matrices and right hand side ... \n');
fprintf('============================================================\n')


% connectivity infos
ndof         = femregion.ndof;
nln          = femregion.nln;
ne           = femregion.ne;
connectivity = femregion.connectivity;

epsilon = femregion.epsilon; 
%(Check epsilon)

% shape functions
[basis] = C_shape_basis(Dati);

% quadrature nodes and weights for integrals
[nodes_2D, w_2D] = C_quadrature(Dati);

% evaluation of shape bases 
[dphiq,Grad] = C_evalshape(basis,nodes_2D);


% Assembly begin ...
A = sparse(ndof,ndof);
M = sparse(ndof,ndof);
ADV = sparse(ndof,ndof);

f = sparse(ndof,1);

cavitation_vector = sparse(ndof,1);
cav_fun = Dati.obstacle;
nu_fun = Dati.nu_fun;

for ie = 1 : ne
     
    iglo = connectivity(1:nln,ie);

    % Determinant of the Jacobian of the map and physical coordinates
    [BJ, pphys_2D] = C_get_Jacobian(femregion.coord(iglo,:), nodes_2D, Dati.MeshType);
    
    % Ogni tris di coordinate globali mi defnisce una matrice BJ diversa
   
    %=============================================================%
    % STIFFNESS AND MASS MATRIX
    %=============================================================%
    [A_loc] = C_lap_loc(Grad,w_2D,nln,BJ,pphys_2D,nu_fun);
    [M_loc_1] = C_mass_loc_1(dphiq,w_2D,nln,BJ,pphys_2D,nu_fun);
    [ADV_loc]=C_adv_loc(Grad,dphiq, Dati.beta,w_2D,nln,BJ);
            
    A(iglo,iglo) = A(iglo,iglo) + A_loc;
    M(iglo,iglo) = M(iglo,iglo) + M_loc_1;
    ADV(iglo,iglo) = ADV(iglo,iglo) + ADV_loc;
    
    %==============================================
    % FORCING TERM --RHS
    %==============================================
    
    
    [load] = C_loc_rhs2D(Dati.force,dphiq,BJ,w_2D,pphys_2D,nln); 
    f(iglo) = f(iglo) + load;
end

M(:,xx) = 0;

x = femregion.coord(:,1);
y = femregion.coord(:,2);
cavitation_vector = eval(cav_fun);


if nargin == 5
    if dummy~=0
        % Se è attivato obstacle
        
        A = A + ADV + M/epsilon;
        f = f + (M*cavitation_vector)/epsilon +u;
        
    else
        A = A;
        f = f;
        
    end
else
    if dummy~=0
        % Se è attivato obstacle
        
        A = A + ADV + M/epsilon;
        f = f + (M*cavitation_vector)/epsilon;
        
    else
        A = A;
        f = f;
        
    end
end







    
