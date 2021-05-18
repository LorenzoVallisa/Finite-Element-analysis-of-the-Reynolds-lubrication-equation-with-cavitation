function [f]=Cav_fun(force,pphys_2D,nln)



f = zeros(nln,1);
% evaluation of the right hand side on the physical nodes
x = pphys_2D(:,1);
y = pphys_2D(:,2);

f = eval(force);

end






                                              
                                              


