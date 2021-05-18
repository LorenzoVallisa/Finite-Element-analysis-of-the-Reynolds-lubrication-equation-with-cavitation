function [errors_table,rates]=C_mesh_test(test_name,ts_num)
%% [errors_table,rates]=C_mesh_test(test_name,ts_num)
%==========================================================================
% Error analysis varying the mesh size h 
%==========================================================================
% Example of usage: [errors_table,rates] = C_mesh_test('Test1',3)
%
%    INPUT:
%          test_name    : (string)  test case name, see C_dati.m
%          ts_num       : (int)     number of refinements
%
%    OUTPUT:
%          errors_table : (struct) containing the computed errors
%          rates        : (struct) containing the computed rates


warning off;
addpath ../.
addpath ../MeshGeneration
addpath ../MeshGeneration/Tria
addpath ../FESpace
addpath ../Assembly
addpath ../BoundaryConditions
addpath ../Postprocessing
addpath ../Errors


Dati=C_dati(test_name);

refinement_vector=Dati.refinement_vector;
num_test=length(refinement_vector);

for k=1:num_test
    [errors,solutions,femregion,Matrices,Dati]=C_main2D(test_name,refinement_vector(k));
    Error_L2(k)=errors.Error_L2;
    Error_SEMI_H1(k)=errors.Error_SEMI_H1;
    Error_H1(k)=errors.Error_H1;
    Error_inf(k)=errors.Error_inf;
    ne(k)=femregion.ne;    
    h(k)=femregion.h;
    fprintf('==========================================\n');    
    fprintf('End test %i\n',k);
    fprintf('==========================================\n');
end
p=femregion.degree;

%ERROR TABLE
errors_table=struct('ne',ne,...
                   'h',h,...
                   'Error_L2', Error_L2,...
                   'Error_SEMI_H1', Error_SEMI_H1,...
                   'Error_H1', Error_H1,...
                   'Error_inf',Error_inf);


               
%TABLE
rate_L2=log10(Error_L2(2:num_test)./Error_L2(1:num_test-1))./log10(h(2:num_test)./h(1:num_test-1));
rate_SEMI_H1=log10( Error_SEMI_H1(2:num_test)./ Error_SEMI_H1(1:num_test-1))./log10(h(2:num_test)./h(1:num_test-1));
rate_H1=log10( Error_H1(2:num_test)./ Error_H1(1:num_test-1))./log10(h(2:num_test)./h(1:num_test-1));
rate_inf=log10( Error_inf(2:num_test)./ Error_inf(1:num_test-1))./log10(h(2:num_test)./h(1:num_test-1));

rates=struct('rate_L2',rate_L2,...
             'rate_SEMI_H1',rate_SEMI_H1,...
             'rate_H1',rate_H1,...
             'rate_inf',rate_inf);
         

% ERROR PLOTS
subplot(2,1,1)
loglog(1./h,h.^(p+1),'-+b');
hold on
loglog(1./h,Error_L2,'-or');
hold on
loglog(1./h,Error_inf,'-*k');
legend(sprintf('1/h^%i',p+1),'||.||_{L^2}','||.||_{inf}');
title ('Errors');
xlabel('1/h');

subplot(2,1,2)
loglog(1./h,h.^(p),'-+b');
hold on
loglog(1./h,Error_SEMI_H1,'-or');
hold on
loglog(1./h,Error_H1,'-*k');
legend(sprintf('1/h^%i',p),'|.|_{H^1}','||.||_{H^1}');
xlabel('1/h');
hold off
 
%-----------------------------------------
% WRITE OUTPUT on file
%----------------------------------------
if Dati.print_out=='Y'
    cd Numerics
    out_name=sprintf('CONF_ERR_%s_%s_Test_%i.m',Dati.MeshType,Dati.fem,ts_num);
    fid = fopen(out_name,'w');  
    cd ..
	C_write_parameters_test(Dati,fid);
    C_write_output(Dati,errors_table,rates,fid);
	fclose(fid);
end

%s-----------------------------------------
% WRITE OUTPUT on DISPLAY
%----------------------------------------
C_write_parameters_test(Dati)
C_write_output(Dati,errors_table,rates)
