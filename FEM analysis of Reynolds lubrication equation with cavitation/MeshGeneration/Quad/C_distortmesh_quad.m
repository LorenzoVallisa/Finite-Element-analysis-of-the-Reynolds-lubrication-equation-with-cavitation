% generate a distorted quadrilateral mesh
% by giggle the internal nodes


function p=C_distortmesh_quad(domain, p, t)


%============================================
% Domain
%============================================
x0=domain(1,1);
x1=domain(1,2);
y0=domain(2,1);
y1=domain(2,2);

h = 1 ./ sqrt(length(t(1,:)));

factor = 0.5;

for ip = 1 : length(p)
    
    x = p(1,ip);
    y = p(2,ip);
    
    
    dx = (-1).^ip.*(h./5).*rand(1);
    dy = (-1).^ip.*(h./5).*rand(1);
    %dx = (-1).^ip.*(h./5).*factor;
    %dy = (-1).^ip.*(h./5).*factor;
    %dx=0;
    %dy=0;   
    
    if (x  ~= x0 & x ~= x1 &  y  ~= y0 & y ~= y1 )
        
        
        x = x + dx;
        y = y + dy;
    end
    
    % boundary point
    if ((y  == y0 | y  == y1) & (x  ~= x0 & x ~= x1))
        
        x = x + dx;
    end
    
    % boundary point
    if ((x  == x0 | x  == x1) & (y  ~= y0 & y ~= y1))
        
        y = y + dy; 
    end
    
   p(1,ip) = x;
   p(2,ip) = y;
end

    
    
        
        
   
        