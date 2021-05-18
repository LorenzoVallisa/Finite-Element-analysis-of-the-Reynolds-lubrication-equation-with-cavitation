function C_write_parameters_test(Dati,fid)
%% C_write_parameters_test(Dati,fid)
%==========================================================================
% Write parameters test case 
%==========================================================================
%    called in C_mesh_tet.m
%
%    INPUT:
%          Dati         : (stuct)   see C_dati.m
%          fid          : (int)     file ID
%

if nargin == 1
    fid=1;
else
    cd Numerics
end
fprintf(fid,'%%=============================================================\n');
fprintf(fid,'%%                PARAMETERS               \n')
fprintf(fid,'%%=============================================================\n');
fprintf(fid,'%%CONFORMING FINITE ELEMENTS\n');    
fprintf(fid,'%%Triangular Grids: %s \n', Dati.MeshType);    
fprintf(fid,'%%Domain=[%2.2f,%2.2f]x[%2.2f,%2.2f]\n',Dati.domain(1,1),Dati.domain(1,2),Dati.domain(2,1),Dati.domain(2,2));
fprintf(fid,'%%Exact Solution: %s \n', Dati.exact_sol);
fprintf(fid,'%%Load: %s \n', Dati.force);
fprintf(fid,'%%FEM: %s\n',Dati.fem);
fprintf(fid,'%%nqn 1D: %2.0d\n',Dati.nqn_1D);
fprintf(fid,'%%nqn 2D: %2.0d\n',Dati.nqn_2D);
fprintf(fid,'%%Linear Solver: %s\n','Direct solver');
fprintf(fid,'%%=============================================================\n');


if nargin == 1
else
    cd ..
end
