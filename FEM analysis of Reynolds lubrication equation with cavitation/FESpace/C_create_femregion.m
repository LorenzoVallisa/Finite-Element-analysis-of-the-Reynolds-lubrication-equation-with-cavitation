function [femregion] = C_create_femregion(Dati,Region) 
%% [femregion] = C_create_femregion(Dati,Region)
%==========================================================================
% Creates conforming finite element space
%==========================================================================
%    called in C_main2D.m
%
%    INPUT:
%          Dati        : (struct)  see C_dati.m
%          Region      : (struct)  see C_create_mesh.m
%
%    OUTPUT:
%          femregion    : (struct) 

fprintf('============================================================\n')
fprintf('Creating finite element space ... \n');
fprintf('============================================================\n')


switch Dati.MeshType
    
    case{'TS', 'TU'}
         %nln = number of local nodes
         degree=1;
         nln=(degree+1).*(degree+2)./2;  
        
        
    case{'QS', 'QU'}
         %nln = number of local nodes
         degree=1;
         nln=(degree+1).^2;
         
    otherwise

    error('The mesh could only be of triangular elements or quadrilateral elements');
end

[bound_pts,per_sx,per_dx] = C_create_bound_pts(Dati,Region);


%==========================================================================
% COSTRUZIONE STRUTTURA FEMREGION
%==========================================================================
femregion=struct('fem',Dati.fem,...
                'domain',Region.domain,...
                'type_mesh',Dati.MeshType,...
                'h',Region.h,...
                'nln',nln,...
                'ndof',length(Region.coord),...
                'ne',Region.ne,...
                'dof',Region.coord,...
                'nqn_1D',Dati.nqn_1D,...
                'nqn_2D',Dati.nqn_2D,...
                'degree',degree,...
                'coord',Region.coord,...
                'connectivity',Region.connectivity,...
                'boundary_points',bound_pts,...
                'per_sx',per_sx,...
                'per_dx',per_dx,....
                'epsilon',0.05*(Region.h)^(2));
            
            
            
function [bound_pts,per_sx,per_dx] = C_create_bound_pts(Dati, Region)
%% [bounds_pts] = C_create_bound_pts(domain, coord)
%==========================================================================
% Creates boundary point list
%==========================================================================
%    called in C_create_femregion.m
%
%    INPUT:
%          domain     : [2 x 2] (array real)  extrema of the rectangle
%          coord      : [Np x 2] (array real) (x,y) coordinates of the grid
%                                             points
%
%    OUTPUT:
%          bound_pts  : [Nbp x 1] (array int) indices of the boundary
%                                             points
coord = Region.coord;
domain = Dati.domain;
x0 = domain(1,1);
x1 = domain(1,2);
y0 = domain(2,1);
y1 = domain(2,2);
bound_pts = ones(length(coord),1);


% Scelgo qui dove metter Dirichlet

 bound_pts =find(coord(:,2)==y0 | coord(:,2)==y1);

 
%Scelgo qui invece i lati sui quali impostare le periodic boundary
% ATTENZIONE: NON basta metter y0 ed y1 per cambiar lato ma vanno invertiti
% tutti gli indici (in alternativa lascia queste definizioni ed inverti il
% problema come abbiamo fatto noi (cambio variabile x-y nelle funzioni in
% dati))

index_x0 = find(coord(:,1)==x0);
index_x1 = find(coord(:,1)==x1);


sorted_x0_var = sort(coord((index_x0),2));
sorted_x1_var = sort(coord((index_x1),2));
new_indexes_x0 =[];
new_indexes_x1= [];
for i = 1:length(sorted_x0_var)
new_indexes_x0 = [new_indexes_x0 find(coord(:,2)==sorted_x0_var(i) & coord(:,1)==x0)];
new_indexes_x1 = [new_indexes_x1 find(coord(:,2)==sorted_x1_var(i) & coord(:,1)==x1)];
end
per_sx = new_indexes_x0';
per_dx = new_indexes_x1';






      
            
            
