clear all
close all 
clc
%in order to determine the numerical convergence of the method, to be faster, values of pressure
%are saved in files px_x.mat, same for femregions.
p3_3 = load ('p3_3.mat');
p3_4 = load ('p3_4.mat');
p3_5 = load ('p3_5.mat');
p3_6 = load ('p3_6.mat');
p3_8 = load ('p3_8.mat');

fem3 = load ('dati3.mat');
fem4 = load ('dati4.mat');
fem5 = load ('dati5.mat');
fem6 = load ('dati6.mat');
fem8 = load ('dati8.mat');

coord_3 = fem3.femregion.coord;
coord_4 = fem4.femregion.coord;
coord_5 = fem5.femregion.coord;
coord_6 = fem6.femregion.coord;
coord_8 = fem8.femregion.coord;

% linear extensions from xcoordinate to 8coordinate are made through
% griddata function!
ext_3 = griddata(coord_3(:,1),coord_3(:,2),full(p3_3.p3),coord_8(:,1),coord_8(:,2));
ext_4 = griddata(coord_4(:,1),coord_4(:,2),full(p3_4.p3),coord_8(:,1),coord_8(:,2));
ext_5 = griddata(coord_5(:,1),coord_5(:,2),full(p3_5.p3),coord_8(:,1),coord_8(:,2));
ext_6 = griddata(coord_6(:,1),coord_6(:,2),full(p3_6.p3),coord_8(:,1),coord_8(:,2));

err_3 = full(p3_8.p3) - ext_3;
err_4 = full(p3_8.p3) - ext_4;
err_5 = full(p3_8.p3) - ext_5;
err_6 = full(p3_8.p3) - ext_6;

%errors are determined through A and M matrixes, in which (in particular for
%M) the obstacle problem is not acting!

A = load ('A_pure.mat');
A = A.A_pure;
M = load ('M_pure.mat');
M = M.M_pure;

K =  A + M;

err(1) = err_3' * K * err_3;
err(2) = err_4' * K * err_4;
err(3) = err_5' * K * err_5;
err(4) = err_6' * K * err_6;

h = 2.^(-[3:6]);

loglog(h,sqrt(err),'o-',h,h,'o-')
legend('||p-p_{nRef=8}||_{H^1}','h^1','Location','NorthWest')
grid on
xlabel('h')
ylabel('Error')