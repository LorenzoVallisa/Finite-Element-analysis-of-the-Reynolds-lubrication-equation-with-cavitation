function [hh,E_H1]= test_convergence_2(sol_struct_vect,final_ref)

for i=3:final_ref
    u_ex = sol_struct_vect(i-2).u_ex;
    uh = sol_struct_vect(i-2).uh;
    solutions = struct('u_ex',u_ex,'uh',uh);
    [unused,femregion,Matrices,Dati,A,M]=C_main2D('Test2',i,'non_sym',0,[]);
    [errors]=C_compute_errors(Dati,femregion,solutions)
%     E_L2(i)=errors.Error_L2;
    E_H1(i)=errors.Error_H1;
%     E_inf(i)=errors.Error_inf;
%     E_DG(i) =errors.E_DG;
    ne(i)=femregion.ne;
    hh(i)=femregion.h;
    fprintf('End test %d\n',i);
end




end