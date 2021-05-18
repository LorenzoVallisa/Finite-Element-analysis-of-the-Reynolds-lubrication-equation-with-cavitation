% quadrature nodes and whieghs for faces integrals
function  [node_2D,w_2D] = C_quadrature(Dati)
%% function [node_2D,w_2D] = C_quadrature(Dati)
%==========================================================================
% Compute quadrature nodes and weights triangular and quadrilateral elements
%==========================================================================
%    called in C_matrix2D.m
%
%    INPUT:
%          Dati        : (struct)  see C_dati.m
%
%    OUTPUT:
%          node2D      : (struct) .num  (int) number of basis functions
%                                 .nedge (int) number of edges
%                                 .fbases (string) basis functions
%                                 .Gbases_1 (string) df/dx
%                                 .Gbases_2 (string) df/dy
%


switch  Dati.MeshType
    
    case{'TS', 'TU'}
        
        [node_2D,w_2D] = C_Tria_int_2D(Dati.nqn_2D);
        
        
    case{'QS', 'QU'}
        
        
        [csi_x,w_x] = C_gauleg(-1,1,Dati.nqn_1D);
        [csi_y,w_y] = C_gauleg(-1,1,Dati.nqn_1D);
        
        k=1;
        for j=1:Dati.nqn_1D
            for i=1:Dati.nqn_1D
                w_2D(k)=w_x(i).*w_y(j);
                node_2D(k,1)=csi_x(i);
                node_2D(k,2)=csi_y(j);
                k=k+1;
            end
        end
        
        
end
