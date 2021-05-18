function [phys_points]=C_Map_phys_points_edge(edge, nqn_1D, node_1D, n_edge, x0, x1,x2,y0,y1,y2)


%================================================================
%  map from the reference quadrure ponit to the physical ones
%================================================================

switch edge
case{1} 
    %================================
    % EDGE 1
    %================================
    if (y1==y0)
        x_e1=((x1-x0).*node_1D + x0)';
        y_e1=ones(nqn_1D,1).*y0;
    elseif (x1==x0)
        x_e1=ones(nqn_1D,1).*x0;
        y_e1=((y1-y0).*node_1D + y0)';
    else ;
        x_e1=zeros(nqn_1D,1);
        y_e1=zeros(nqn_1D,1);
    end
    phys_points(:,1)=x_e1;
    phys_points(:,2)=y_e1;
  
case{2} 
    
    %================================
    % EDGE 2
    %================================
    if (y2==y1)
        x_e2=((x2-x1).*node_1D + x1)';
        y_e2=ones(nqn_1D,1).*y1;
    elseif (x2==x1)
        x_e2=ones(nqn_1D,1).*x1;
        y_e2=((y2-y1).*node_1D + y1)';
    else ;
        x_e2=zeros(nqn_1D,1);
        y_e2=zeros(nqn_1D,1);
    end
    phys_points(:,1)=x_e2;
    phys_points(:,2)=y_e2;
    
case{3} 
    %================================
    % EDGE 3
    %================================
    if (y2==y0)
        x_e3=((x0-x2).*node_1D + x2)';
        y_e3=ones(nqn_1D,1).*y2;
    elseif (x2==x0)
        x_e3=ones(nqn_1D,1).*x2;
        y_e3=((y0-y2).*node_1D + y2)';
    else ;
        x_e3=zeros(nqn_1D,1);
        y_e3=zeros(nqn_1D,1);
    end
    phys_points(:,1)=x_e3;
    phys_points(:,2)=y_e3;
end