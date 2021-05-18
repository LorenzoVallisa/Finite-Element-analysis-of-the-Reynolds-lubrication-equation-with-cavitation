function [errors]=C_compute_errors(Dati,femregion,solutions)
%% [errors]=C_compute_errors(Dati,femregion,solutions)
%==========================================================================
% Compute L2, semi-H1, H1 and L-inf errors
%==========================================================================
%    called in C_main2D.m
%
%    INPUT:
%          Dati        : (struct)  see C_dati.m
%          femregion   : (struct)  see C_create_femregion.m
%          solutions   : (struct)  see C_postprocessing.m
%
%    OUTPUT:
%          errors      : (struct)  


err = solutions.uh-solutions.u_ex;
[E_L2, E_SEMI_H1] = C_error_L2_H1(femregion, solutions.uh, Dati);

errors = struct('Error_L2',E_L2,...
              'Error_SEMI_H1',E_SEMI_H1,...
              'Error_H1', sqrt(E_L2.^2 + E_SEMI_H1.^2),...
              'Error_inf', norm (err,inf));