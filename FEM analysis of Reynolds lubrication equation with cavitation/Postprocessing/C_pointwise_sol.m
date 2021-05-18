function C_pointwise_sol(femregion, uh, u_ex, Dati)
%% C_pointwise_sol(femregion, uh, u_ex,Dati)
%==========================================================================
% PLOT THE EXACT SOLUTION ON THE DOFS
%==========================================================================
%    called in C_postprocessing.m
%
%    INPUT:
%          femregion   : (struct)  see C_create_femregion.m
%          uh          : (sparse(ndof,1) real) solution vector
%          u_ex        : (sparse(ndof,1) real) exact solution vector
%          Dati        : (struct) see C_Dati.m
%



dof=femregion.dof;

x1=femregion.domain(1,1);
x2=femregion.domain(1,2);
y1=femregion.domain(2,1);
y2=femregion.domain(2,2);

M=max(u_ex);
m=min(u_ex);
if (abs(m-M) < 0.1)
    M=m+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% VISUALIZZAZIONE DELLE SOLUZIONI MESHATE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% trisurf(femregion.connectivity(1:3,:)',femregion.coord(:,1),femregion.coord(:,2),full(uh));
% title('u_h(x,y)'); xlabel('x-axis'); ylabel('y-axis');

%axis([x1,x2,y1,y2,m,M]); 
% colorbar;



% figure;
% stem3(dof(:,1),dof(:,2),u_ex,'-r*');
% hold on
% stem3(dof(:,1),dof(:,2),uh,'-ko');
% title(Dati.name);
% %axis([x1,x2,y1,y2,m,M]);
% xlabel('x-axis'); ylabel('y-axis');
% legend('exact', 'computed')

