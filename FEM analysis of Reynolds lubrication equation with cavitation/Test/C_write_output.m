function C_write_output(Dati,errors_table,rates,fid)
%% C_write_output(Dati,errors_table,rates,fid)
%==========================================================================
% Write outputs 
%==========================================================================
%    called in C_mesh_tet.m
%
%    INPUT:
%          Dati         : (struct)  see C_dati.m
%          error_table  : (struct)  containing computed errors
%          rates        : (struct)  containint computed rates of conv.
%          fid          : (int)     file ID
%
%
%


num_test=length(errors_table.ne);
if nargin ==3
    fid=1;
else
    cd Numerics
end

fprintf(fid,'%%=============================================================\n');
fprintf(fid,'%%                Errors Table                              \n');
fprintf(fid,'%%=============================================================\n');
fprintf(fid,'%%\t N \t \t L2 \t \t inf\n')
fprintf(fid,'%%-------------------------------------------------------------\n');
for i=1:num_test
    fprintf(fid,'%5.0i  \t %1.5e \t%1.5e \n',errors_table.ne(i),errors_table.Error_L2(i),errors_table.Error_inf(i));
end
fprintf(fid,'%%-------------------------------------------------------------\n');
fprintf(fid,'%%\t N \t \t H1 \t \t 1,h \n')
fprintf(fid,'%%-------------------------------------------------------------\n');
for i=1:num_test
    fprintf(fid,'%5.0i \t %1.5e \t%1.5e \n',errors_table.ne(i),errors_table.Error_H1(i),errors_table.Error_SEMI_H1(i));
end
fprintf(fid,'%%-------------------------------------------------------------\n');
fprintf(fid,'\n');
fprintf(fid,'%%=============================================================\n');
fprintf(fid,'%%                 rates                             \n');
fprintf(fid,'%%=============================================================\n');
fprintf(fid,'%%L2 \t \t inf\n')
fprintf(fid,'%%-------------------------------------------------------------\n');
for i=1:num_test-1
    fprintf(fid,'%1.5f \t%1.5f \n',rates.rate_L2(i),rates.rate_inf(i));
end
fprintf(fid,'%%-------------------------------------------------------------\n');
fprintf(fid,'%%H1 \t \t SEMI_H1 \n')
fprintf(fid,'%%-------------------------------------------------------------\n');
for i=1:num_test-1
    fprintf(fid,'%1.5f \t%1.5f \n',rates.rate_H1(i),rates.rate_SEMI_H1(i));
end
fprintf(fid,'%%-------------------------------------------------------------\n');

if nargin == 3
else
    cd ..
end
