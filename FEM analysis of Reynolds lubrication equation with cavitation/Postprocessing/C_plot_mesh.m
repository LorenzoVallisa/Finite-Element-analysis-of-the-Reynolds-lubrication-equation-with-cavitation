function C_plot_mesh(femregion)
%% C_plot_mesh(femregion)
%==========================================================================
% PLOT of THE MESH
%==========================================================================
%    called in C_postprocessing.m
%
%    INPUT:
%          femregion   : (struct)  see C_create_femregion.m

ne = femregion.ne;
connectivity = femregion.connectivity;
coord = femregion.coord;

for i=1:ne
  if length(connectivity(:,i)) == 4 
    fill([coord(connectivity(1,i),1) coord(connectivity(2,i),1) ...
        coord(connectivity(3,i),1) coord(connectivity(1,i),1)], ...
       [coord(connectivity(1,i),2) coord(connectivity(2,i),2) ...
        coord(connectivity(3,i),2) coord(connectivity(1,i),2)],'w')

  else

    fill([coord(connectivity(1,i),1) coord(connectivity(2,i),1) ...
        coord(connectivity(3,i),1) coord(connectivity(4,i),1) ...
        coord(connectivity(1,i),1)], ...
       [coord(connectivity(1,i),2) coord(connectivity(2,i),2) ...
        coord(connectivity(3,i),2) coord(connectivity(4,i),2) ...
        coord(connectivity(1,i),2)],'w')
  end
  hold on
end