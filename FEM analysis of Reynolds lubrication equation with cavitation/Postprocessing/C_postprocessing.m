function [solutions]=C_postprocessing(Dati,femregion,uh)
%% [solutions]=C_postprocessing(Dati,femregion,uh)
%==========================================================================
% POST PROCESSING OF THE SOLUTION
%==========================================================================
%    called in C_main2D.m
%
%    INPUT:
%          Dati        : (struct)  see C_dati.m
%          femregion   : (struct)  see C_create_femregion.m
%          uh          : (sparse(ndof,1) real) solution vector
%
%    OUTPUT:
%          solutions   : (struct) containg solution vector uh and  
%                        analytical solution u_ex
%

fprintf('============================================================\n')
fprintf('Post-processing the solution ... \n');
fprintf('============================================================\n')


%==========================================================================
% EVALUATION OF THE EXACT SOLUTION
%==========================================================================

[u_ex] = C_eval_exact_sol(femregion,Dati.exact_sol);

%==========================================================================
% PLOT MESH
%==========================================================================

%C_plot_mesh(femregion);

%==========================================================================
% PLOT SOLUTION
%==========================================================================

C_pointwise_sol(femregion, uh, u_ex, Dati);

%==========================================================================
% SAVE SOLUTIONS
%==========================================================================

solutions = struct('u_ex',u_ex,'uh',uh);