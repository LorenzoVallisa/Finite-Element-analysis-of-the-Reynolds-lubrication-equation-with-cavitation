function [node_2D,w_2D,nqn_2D]=C_Tria_int_2D(nqn_2D)
%% [node_2D,w_2D,nqn_2D]=C_Tria_int_2D(nqn_2D)
%==========================================================================
% Compute Gauss-Legendre nodes and weights on the reference triangle
%==========================================================================
%    called in C_quadrature.m
%
%    INPUT:
%          nqn_2D      : (integer) number of quadrature nodes
%
%    OUTPUT:
%          node_2D     : (array) [nqn_2D x 2] list of quadrature nodes
%          w_2D        : (array) [nqn_2D x 1] list of quadrature weights        
%          nqn_2D      : (integer) number of quadrature nodes
       

mis_T=0.5;
third=0.333333333333333333333333333333333333;  
switch nqn_2D 
    
case{1}
    %===================================================
    % 1 node quadrature formula on the reference triangle
    % degree of precision: 1
    %===================================================
    
    xnod=third;
    ynod=third;
    w=1;
    
    w=mis_T.*w;
    
case{3}
    
    %===================================================
    % 3 node quadrature formula on the reference triangle
    % degree of precision: 2
    %===================================================
    xnod=zeros(1,3);
    ynod=zeros(1,3);
    w=zeros(1,3);
    
    w(1)=third;
    w(2)=third;
    w(3)=third;
    
    xnod(1)=0.5;
    xnod(2)=0.5;
    xnod(3)=0;
    
    ynod(1)=0;
    ynod(2)=0.5;
    ynod(3)=0.5;
    
    w=mis_T.*w;
    
case{4}
    %===================================================
    % 4 node quadrature formula on the reference triangle
    % degree of precision: 3
    %===================================================
    xnod=zeros(1,4);
    ynod=zeros(1,4);
    w=zeros(1,4);
        
    a=25/48;
    b=9/16;
    
    
    w(1)=a;
    w(2)=a;
    w(3)=a;
    w(4)=-b;
    
    xnod(1)=0.2;
    xnod(2)=0.6;
    xnod(3)=0.2;
    xnod(4)=third;
    
    ynod(1)=0.2;
    ynod(2)=0.2;
    ynod(3)=0.6;
    ynod(4)=third;
    
    w=mis_T.*w;
    
case{7}
    %===================================================
    %7 node quadrature formula on the reference triangle
    % degree of precision: 4
    %===================================================
    % 
    
    xnod=zeros(1,7);
    ynod=zeros(1,7);
    w=zeros(1,7);
    
    
    a=0.0597158717;
    b=0.4701420641;
    c=0.7974269853;
    d=0.1012865073d0;
    
    w(1)=0.225;
    w(2)=0.1323941527;
    w(3)=0.1323941527;
    w(4)=0.1323941527;
    w(5)=0.1259391805;
    w(6)=0.1259391805;
    w(7)=0.1259391805;
    
    lam=zeros(3,7);
    
    lam(1,1)=third;
    lam(2,1)=third;
    lam(3,1)=third;
    
    lam(1,2)=a;
    lam(2,2)=b;
    lam(3,2)=b;  
    
    lam(1,3)=b;
    lam(2,3)=b;
    lam(3,3)=a;  
    
    lam(1,4)=b;
    lam(2,4)=a;
    lam(3,4)=b;
    
    lam(1,5)=c;
    lam(2,5)=d;
    lam(3,5)=d;
    
    lam(1,6)=d;
    lam(2,6)=d;
    lam(3,6)=c;
    
    lam(1,6)=d;
    lam(2,6)=d;
    lam(3,6)=c;
    
    lam(1,7)=d;
    lam(2,7)=c;
    lam(3,7)=d;
    
    xnod=lam(1,:);
    ynod=lam(2,:); 
    
    w=mis_T.*w;

   
    case{10}
    %===================================================
    %10 node quadrature formula on the reference triangle
    % degree of precision: 5
    %===================================================
    % 
    
    xnod=zeros(1,10);
    ynod=zeros(1,10);
    w=zeros(1,10);
    
    
    a=0;
    b=0.5;
    c=(1/3);
    d=(1/7);
    e=(5/7);
    f=1;
    
   w=[1/90, 1/90, 1/90, 16/225 , 16/225, 16/225, (49/120).^2 , (49/120).^2 , (49/120).^2 , 81/320];
   
   xnod=[a,f,a,b,b,a,d,e,d,c];
   ynod=[a,a,f,a,b,b,d,d,e,c];
    
   w=mis_T.*w;
   
   
    case{16}
    %===================================================
    % 16 node quadrature formula on the reference triangle
    % degree of precision: 7
    %===================================================
    % 
    
    xnod=[];
    ynod=[];
    w=[];
    
    % points are (s(j),r(i)(1-s(j))) with i,j=1,..,4
    % respective weights are A(i)B(j) with i,j=1,..,4   
   
    [r,A]=C_gauleg(0,1,4);
   
    s=[0.0571041961,...
       0.2768430136,...
       0.5835904324,...
       0.8602401357];

    
    B=[0.1355069134,...
       0.2034645680,...
       0.1298475476,...
       0.0311809709];
        
    
    for i=1:4
        for j=1:4
            xnod=[xnod,s(j)];
            ynod=[ynod,r(i).*(1-s(j))];
            w=[w, A(i)*B(j)];
        end
    end
    
    
    case{64}
    %===================================================
    % 64 node quadrature formula on the reference triangle
    % degree of precision: 15
    %===================================================
      
    xnod=[];
    ynod=[];
    w=[];
    
    % points are (s(j),r(i)(1-s(j))) with i,j=1,..,8
    % respective weights are A(i)B(j) with i,j=1,..,8   
   [r,A]=C_gauleg(0,1,8);

   s=[     0.9644401697,...
           0.8173527842,...
           0.5713830412,...
           0.2561356708,...
           -0.09037336960,...
           -0.4263504857,...
           -0.7112674859,...
           -0.9107320894];
   
   
   B=[     0.1782032174,...
           0.3644760945,...
           0.4500231978,...
           0.4241894377,...
           0.3167983979,....
           0.1817572780,...
           0.07137161062...
           0.01318076576];

   
s=0.5.*(ones(1,8)-s);
B=B./4;

%    s=[0.04463395529,...
%       0.14436625704,...
%       0.28682475714,...
%       0.45481331520,...
%       0.62806783542,...
%       0.78569152060,...
%       0.90867639210,...
%       0.98222008485];
%    
%     B=[0.05061426814,...
%        0.11119051723,...
%        0.15685332294,...
%        0.18134189169,...
%        0.18134189169,...
%        0.15685332294,...
%        0.11119051723,... 
%        0.05061426814];

     for i=1:8
        for j=1:8
            xnod=[xnod,s(j)];
            ynod=[ynod,r(i).*(1-s(j))];
            w=[w, A(i)*B(j)];
        end
    end
end



node_2D=[xnod',ynod'];
w_2D=w;    
