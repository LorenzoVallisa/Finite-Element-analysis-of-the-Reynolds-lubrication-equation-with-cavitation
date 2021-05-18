function [Region] = C_create_mesh(Dati)
%% [Region] = C_create_mesh(Dati)
%==========================================================================
% Creates triangular or quadrilateral mesh
%==========================================================================
%    called in C_main2D.m
%
%    INPUT:
%          Dati        : (struct)  see C_dati.m
%
%    OUTPUT:
%          Region      : (struct) having fields: dimension
%                                                mesh type 
%                                                domain 
%                                                mesh size
%                                                number of vertices
%                                                number of elements
%                                                coordinates
%                                                boundary edges 
%                                                connectivity


addpath MeshGeneration/Quad
addpath MeshGeneration/Tria



switch Dati.MeshType
    % g = describes the geometry of the PDE problem (dummy)
    % p = points of the grid
    % be = boundary egdes
    % t = triangulation
    
    case{'TS', 'TU'}
        
         [g, p, be, t] = C_create_mesh_tria(Dati);
    case{'QS', 'QU'}
        
         [g, p, be, t] = C_create_mesh_quad(Dati);
         
    otherwise

        error('The mesh could only be of triangular elements or quadrilateral elements');
end


%================================================
% GEOMETRICAL INFO
 nEl = size(t,2);
 nVert =size(p,2);
 MeshSize = 1/sqrt(nEl);
%================================================

% struttura dati della mesh
Region = struct('dim',2,...
               'MeshType',Dati.MeshType,...
               'domain',Dati.domain,...
               'h',MeshSize,...
               'nvert',nVert,...
               'ne',nEl,...
               'coord_x',p(1,:)',...
               'coord_y',p(2,:)',...
               'coord',p',...
               'boundary_edges',be,...
               'connectivity',t);
           
           
