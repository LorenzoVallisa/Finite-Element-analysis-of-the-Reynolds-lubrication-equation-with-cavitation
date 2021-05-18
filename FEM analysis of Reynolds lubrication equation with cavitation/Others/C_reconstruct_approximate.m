function [uh_rec] = C_reconstruct_approximate(uh,femregion,Dati)
%% [uh_rec] = C_reconstruct_approximate(uh,femregion,Dati)
%==========================================================================
% RECONSTRUCT THE SOLUTION ON THE GRID NODES
%==========================================================================
%    called in C_postprocessing.m
%
%    INPUT:
%          uh          : (spare(ndof,1) real) solution vector
%          femregion   : (struct) see C_create_femregion.m
%          Dati        : (struct) see C_Dati.m
%
%    OUTPUT:
%          uh_rec      : (sparse(ndof,1) real) solution vector


%fem=Dati.fem;
%shape_base=Dati.shape_base;
%[base]=T_HP_bases2D_HP(fem,shape_base);


[base] = C_shape_basis(Dati);

dof_Q1 = femregion.dof_Q1;
dof = femregion.dof;

nln = femregion.nln;
ne = femregion.ne(1)*femregion.ne(2);

uh_rec=[];

for ie=1:ne
    index=(ie-1)*nln*ones(nln,1) + [1:nln]';
    i_Q1=[4*(ie-1)+1,4*(ie-1)+2,4*(ie-1)+3,4*(ie-1)+4];
    
    % determinante dello Jacobiano della trasformazione
    a=dof_Q1(i_Q1(1),1);
    b=dof_Q1(i_Q1(2),1);
    c=dof_Q1(i_Q1(1),2);
    d=dof_Q1(i_Q1(3),2);
    
    len_x=(b-a)/2;
    len_y=(d-c)/2;
    Jdet=len_x*len_y;
    
    uh_loc=uh(index);
    dof_loc=dof(index,:);
    
    % this translate the phisical nodes to the reference ones
    x_hat=(2.*dof_loc(:,1) -a - b)./(b-a);
    y_hat=(2.*dof_loc(:,2) -d - c)./(d-c);
    
    for k=1:nln
        temp=0;
        csi= x_hat(k);
        eta= y_hat(k);
        for ell=1:nln
            temp=temp + uh_loc(ell).*eval(base(ell).fbases);
        end
        uh_rec = [uh_rec; temp];
    end
end

fprintf('\n========================================================================\n');
fprintf(' RECONSTRUCTER OUTPUT:\n');
fprintf('I have reconstructed the approximate solution in the nodal vertex \n');
fprintf('========================================================================\n');
