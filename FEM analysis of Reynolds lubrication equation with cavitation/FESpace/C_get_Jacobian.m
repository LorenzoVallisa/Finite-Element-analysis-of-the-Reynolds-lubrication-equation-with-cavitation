function [BJ, pphys_2D] = C_get_Jacobian(loc_coord, nodes_2D, type_mesh)
%% [BJ, pphys_2D] = C_get_Jacobian(loc_coord, nodes_2D, type_mesh)
%==========================================================================
% Compute Jacobian of the map transformation and physical nodes
%==========================================================================
%    called in C_matrix2D.m
%
%    INPUT:
%          loc_coord   : (array real) coordinates of grid nodes
%          nodes_2D    : (struct) see C_quadrature.m
%          type_mesh   : (string) mesh type 
%
%    OUTPUT:
%          BJ          :  (array real) [2 x 2 x nNode] Jacobian
%          pphys_2D    :  (array real) map of the reference element nodes     



for k=1:length(nodes_2D)


    switch  type_mesh

        case{'TS', 'TU'}

            x0=loc_coord(1,1);   % x-coordinates of vertices
            x1=loc_coord(2,1);
            x2=loc_coord(3,1);

            y0=loc_coord(1,2);   % y-coordinates of vertices
            y1=loc_coord(2,2);
            y2=loc_coord(3,2);

            BJ(:,:,k) = [x1-x0,x2-x0;y1-y0,y2-y0];   % Jacobian of elemental map

            trans=[x0;y0];                      % translation vector

            pphys_2D(k,:) = transpose((BJ(:,:,k)*transpose(nodes_2D(k,:))+trans));

        case{'QS', 'QU'}
            
            x0=loc_coord(1,1);   % x-coordinates of vertices
            x1=loc_coord(2,1);
            x2=loc_coord(3,1);
			x3=loc_coord(4,1);
	
            y0=loc_coord(1,2);   % y-coordinates of vertices
            y1=loc_coord(2,2);
            y2=loc_coord(3,2);
			y3=loc_coord(4,2);
	
	
			S = (0.25) .* [-x0  - x1 + x2  + x3,  x0 - x1 - x2  + x3; -y0 - y1 + y2 + y3  , y0 - y1  - y2 + y3];
			
			r     = (0.25) .* [-x0 + x1 - x2  + x3 ; -y0 + y1  - y2 + y3];                       % bilinear vector
			trans = (0.25) .* [ x0 + x1 + x2  + x3 ;  y0 + y1  + y2 + y3];                       % translation vector
			 
            BJ(:,:,k) = S + r*[nodes_2D(k,2) , nodes_2D(k,1)] ;       % Jacobian of elemental map

			pphys_2D(k,:) = transpose( S * transpose(nodes_2D(k,:))  + r * nodes_2D(k,1)*nodes_2D(k,2) + trans);
    end

    
end


 