function [solutions,femregion,Matrices,Dati,A,M_no_bc]=C_main2D(TestName,nRef,bc,dummy,xx,u)
%% [errors,solutions,femregion,Matrices,Dati]=C_main2D(TestName,nRef)
%==========================================================================
% Solution of the Poisson's problem with linear finite elements
% (non homogeneous Dirichlet boundary conditions)
%==========================================================================
%
%    INPUT:
%          Dati        : (struct)  see C_dati.m
%          nRef        : (int)     refinement level
%
%    OUTPUT:
%          errors      : (struct) contains the computed errors
%          solutions   : (sparse) nodal values of the computed and exact
%                        solution
%          femregion   : (struct) infos about finite elements
%                        discretization
%          Matrices    : (struct) fe stiffness and mass matrices
%          Dati        : (struc)  see C_dati.m
%          

warning off;
addpath MeshGeneration
addpath FESpace
addpath Assembly
addpath BoundaryConditions
addpath Postprocessing
addpath Errors


%==========================================================================
% LOAD DATA FOR TEST CASE
%==========================================================================

Dati = C_dati(TestName);
Dati.nRefinement = nRef;

%==========================================================================
% MESH GENERATION
%==========================================================================

[Region] = C_create_mesh(Dati);

%==========================================================================
% FINITE ELEMENT REGION
%==========================================================================

[femregion] = C_create_femregion(Dati,Region); 

%==========================================================================
% BUILD FINITE ELEMENT MATRICES and RIGHT-HAND SIDE
%==========================================================================
if nargin == 6
    [A_no_bc,b_no_bc,M_no_bc] = C_matrix2D(Dati,femregion,dummy,xx,u);
else
    [A_no_bc,b_no_bc,M_no_bc] = C_matrix2D(Dati,femregion,dummy,xx);
end
%(Check dummy)

    

%==========================================================================
% COMPUTE BOUNDARY CONDITIONS
%==========================================================================

[A,b,u_g,M] = C_bound_cond2D(A_no_bc,b_no_bc,M_no_bc,femregion,Dati,bc,TestName);
%(Check A)

%==========================================================================
% SAVE DATA in THE STRUCT Matrices
%==========================================================================

[Matrices] = C_build_out_matrices(A,A_no_bc,M,M_no_bc,b,b_no_bc,u_g);

%==========================================================================
% SOLVE THE LINEAR SYSTEM
%==========================================================================
    
    uh = A\b;

%==========================================================================
% ASSIGN DIRICHLET BOUNDARY CONDITIONS
%==========================================================================

uh = uh + u_g;

%==========================================================================
% POST-PROCESSING OF THE SOLUTION
%==========================================================================

[solutions] = C_postprocessing(Dati,femregion,uh);


end



