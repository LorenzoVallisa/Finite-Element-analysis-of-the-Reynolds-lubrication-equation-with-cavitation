function [E_L2, E_SEMI_H1] = C_error_L2_H1(femregion,uh,Dati)

% Prendo soluzione con numero di refinement molto alto, la traslo su mesh
% con nRefinement più piccolo di modo da poter comparare le due soluzioni


nln=femregion.nln;
ne=femregion.ne;

%type mesh
type_mesh = Dati.MeshType;



% shape functions
[basis] = C_shape_basis(Dati);

% quadrature nodes and whieghs for  integrals
[nodes_2D,w_2D] = C_quadrature(Dati);

% evaluation of shape bases 
[dphiq,Grad] = C_evalshape(basis, nodes_2D);

E_SEMI_H1_LOC = zeros(ne,1);
E_L2_LOC = zeros(ne,1);
%%%%%%%%%%%%%%%%% CICLO SUGLI ELEMENTI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for ie=1:ne
         
    i=femregion.connectivity(1:nln,ie);   
    
    
    % determinante dello Jacobiano della trasformazione
    % and the physical coordinates 
    [BJ, pphys_2D] = C_get_Jacobian(femregion.coord(i,:), nodes_2D, type_mesh);
        
    local_uh=uh(i);
    
    % stima della soluzione esatta e del gradiente della soluzione esatta nei nodi fisici
	x=pphys_2D(:,1);
    y=pphys_2D(:,2);
%     epsilon = Dati.epsilon;
    local_exact=eval(Dati.exact_sol)';
    local_grad_exact=[eval(Dati.grad_exact_1),eval(Dati.grad_exact_2)];
       
    % ricostruzione soluzione approssimata e gradiente della soluzione approssimata nei nodi di quadratura
    local_grad_aprox=zeros(length(w_2D),2);
    local_aprox=zeros(1,length(w_2D));
    
    for k=1:length(w_2D)
        for s=1:nln
             local_aprox(k)=local_aprox(k) + dphiq(1,k,s).*local_uh(s);
             local_grad_aprox(k,:)=local_grad_aprox(k,:)+ Grad(k,:,s).*local_uh(s);
        end
    end
    
    % Stima integrali
    for k=1:length(w_2D)
        Jdet=det(BJ(:,:,k));    % determinant
        dx = abs(Jdet).*w_2D(k); % wheight
        Binv=inv(BJ(:,:,k));    % inverse
 
        pointwise_diff(k,:)=(local_grad_exact(k,:))-(local_grad_aprox(k,:)*Binv);
          
            
        E_SEMI_H1_LOC(ie)= E_SEMI_H1_LOC(ie) + (pointwise_diff(k,:) * transpose(pointwise_diff(k,:))).*dx;
        E_L2_LOC(ie)=E_L2_LOC(ie) + ((local_aprox(k)-local_exact(k)).^2).*dx;
    end
end

% assegnazione variabili OUTPUT
E_SEMI_H1=sqrt(sum(E_SEMI_H1_LOC));
E_L2=sqrt(sum(E_L2_LOC));
