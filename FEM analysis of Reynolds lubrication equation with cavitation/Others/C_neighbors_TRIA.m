%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Crea la matrice topologica dei vicini: specifica per ogni 
% elemento i suoi vicini attraverso ciascuno dei suoi lati
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Determines neighboring elements 
% and number of neighboring edges of 
% a list of elements given
% by connectivity list t
% 
% neigh(i,j) = jth neighbor of element i
% neighedges(i,j) = number of edge of jth neighbor of element i
%  
% 1 \leq i \leq number of elements
% 1 \leq j \leq 3
%
% value: -1 if boundary edge...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [neighbor]=T_HP_neighbors_TRIA(femregion)

fprintf('============================================\n')
fprintf('>>>>>>>>> determine neighbors\n')
fprintf('============================================\n')
%

NEL=femregion.ne;
tria_matrix=femregion.tria_matrix;

neigh=-ones(NEL,3);
neighedges=-ones(NEL,3);

for i=1:(NEL-1)
   
  v1=tria_matrix(1,i);
  v2=tria_matrix(2,i);
  v3=tria_matrix(3,i);

  for j=(i+1):NEL

      vn1=tria_matrix(1,j);
      vn2=tria_matrix(2,j);
      vn3=tria_matrix(3,j);
      
      flag=0;
      if (v1==vn2 & v2==vn1)
        neigh(i,1)=j;
        neigh(j,1)=i;
        neighedges(i,1)=1;
        neighedges(j,1)=1;
      elseif (v1==vn1 & v2==vn3) 
        neigh(i,1)=j;
        neigh(j,3)=i;
        neighedges(i,1)=3;
        neighedges(j,3)=1;
      elseif (v1==vn3 & v2==vn2)
        neigh(i,1)=j;
        neigh(j,2)=i;
        neighedges(i,1)=2;
        neighedges(j,2)=1;
      elseif (v2==vn2 & v3==vn1)
        neigh(i,2)=j;
        neigh(j,1)=i;
        neighedges(i,2)=1;
        neighedges(j,1)=2;
      elseif (v2==vn1 & v3==vn3) 
        neigh(i,2)=j; 
        neigh(j,3)=i; 
        neighedges(i,2)=3;
        neighedges(j,3)=2;
      elseif (v2==vn3 & v3==vn2)
        neigh(i,2)=j;
        neigh(j,2)=i;
        neighedges(i,2)=2;
        neighedges(j,2)=2;
      elseif (v3==vn2 & v1==vn1)
        neigh(i,3)=j;
        neigh(j,1)=i;
        neighedges(i,3)=1;
        neighedges(j,1)=3;
      elseif (v3==vn1 & v1==vn3) 
        neigh(i,3)=j;
        neigh(j,3)=i;
        neighedges(i,3)=3;
        neighedges(j,3)=3;
      elseif (v3==vn3 & v1==vn2)
        neigh(i,3)=j;
        neigh(j,2)=i;
        neighedges(i,3)=2;
        neighedges(j,2)=3;
      end

   end  
end


%===================================================================================
% COSTRUZIONE STRUTTURA NEIGHBOUR
%===================================================================================
neighbor=struct('neigh',neigh,...
                 'neighedges',neighedges,...
                 'nedges',3);
                
%===================================================================================
% END COSTRUZIONE STRUTTURA NEIGHBOUR
%===================================================================================