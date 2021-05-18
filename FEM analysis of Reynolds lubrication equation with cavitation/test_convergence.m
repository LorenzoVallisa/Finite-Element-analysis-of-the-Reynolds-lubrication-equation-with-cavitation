function [table_error,rates]= test_convergence(sol_struct_vect,final_ref)

for i=3:final_ref
    u_ex = sol_struct_vect(i-2).u_ex;
    uh = sol_struct_vect(i-2).uh;
    solutions = struct('u_ex',u_ex,'uh',uh);
    [unused,femregion,Matrices,Dati,A,M]=C_main2D('Test2',i,'non_sym',0,[]);
    [errors]=C_compute_errors(Dati,femregion,solutions)
    E_L2(i)=errors.Error_L2;
    E_H1(i)=errors.Error_H1;
    E_inf(i)=errors.Error_inf;
%     E_DG(i) =errors.E_DG;
    ne(i)=femregion.ne;
    hh(i)=femregion.h;
    fprintf('End test %d\n',i);
end

% close all

table_error=[ne',hh',E_L2',E_H1', E_inf'];
rates=[log10(E_L2(2:end)./E_L2(1:end-1))', log10(E_H1(2:end)./E_H1(1:end-1))', log10(E_inf(2:end)./E_inf(1:end-1))']./log10(hh(2:end)/hh(1:end-1));

figure
grid on
loglog(1./hh,hh.^(femregion.degree+1),'-*k');
hold on
loglog(1./hh,hh.^femregion.degree,'-sr');
hold on
loglog(1./hh,E_L2,'-oc');
hold on
loglog(1./hh,E_H1,'-dm');
hold on
loglog(1./hh,E_inf,'-pg');
% hold on
% loglog(1./hh,E_DG,'-pb');

legend(sprintf('h^%d',femregion.degree+1),...
       sprintf('h^%d',femregion.degree),...
       '||u-u_h||_{0,\Omega}',...
       '||u-u_h||_{1,\Omega}',...
       '||u-u_h||_{inf}',...
       'Location','BestOutside')
end


