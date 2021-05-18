function [g, p, be, t] = C_create_mesh_quad(Dati)
%% [g, p, be, t] = C_create_mesh_quad(Dati)
%==========================================================================
% Creates triangular mesh
%==========================================================================
%    called in C_create_mesh.m
%
%    INPUT:
%          Dati        : (struct)  see C_dati.m
%
%    OUTPUT:
%          g           : (array real) geometry of the grid
%          p           : (array real  [2 x ...] list of the grid points
%          be          : (array int)  [2 x ...] the first and second rows
%                                     contain indices of the starting and
%                                     ending point of the boundary edges
%          t           : (array int)  [5 x Nel] t(1:4,:) indices to
%                                     the corner points, given in counter
%                                     clockwise order, t(5,:)contains the
%                                     subdomain number.

be =[]; % not implemented
warning('boundary elements structure for quadrilateral meshes non implemented');
p=[];
t=[];

%============================================
% Domain
%============================================
x0=Dati.domain(1,1);
x1=Dati.domain(1,2);
y0=Dati.domain(2,1);
y1=Dati.domain(2,2);

g=[
    2     2     2     2
    x0    x1    x1    x0
    x1    x1    x0    x0
    y1    y1    y0    y0
    y1    y0    y0    y1
    0     0     0     0
    1     1     1     1
    ];


if Dati.nRefinement == 0
    
    %============================================
    % points
    %============================================
    p = [
        x1     x1     x0    x0
        y0     y1     y0    y1
        ];
    %============================================
    % boundary
    %============================================
    be = [  ];  
    
    
    
    %============================================
    % SQUARE
    %============================================
   
    t =[
    1
    2
    4
    3
    1
    ];

else
    
    nx = 2^Dati.nRefinement +1;
    ny = 2^Dati.nRefinement +1;
    ne = 4^Dati.nRefinement  ;
    
    % coordinates
    x = x0 : 1/(nx-1): x1;
    y = y0 : 1/(ny-1): y1;
    
    for i = 1 : nx
        for j=1: ny
            p = [p; x(end+1-i), y(j)];
        end
    end
    p = transpose(p);
    
    % connectivity
    k=1;
    counter=0;
    for i=1:ny-1
        for j=1:nx-1
            counter=counter+1;
            t(:,counter)=[j+k-1,j+k,k+ny-1+j+1,k+ny-1+j,1];
        end
        k=k+nx;
    end
    
    
end



switch Dati.MeshType
    case{'QU'}
        p=C_distortmesh_quad(Dati.domain, p, t);
end


