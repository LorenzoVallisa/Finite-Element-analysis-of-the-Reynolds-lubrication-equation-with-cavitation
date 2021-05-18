
function [ADV_loc]=C_adv_loc(Grad,dphiq,beta,w_2D,nln,BJ)

ADV_loc=sparse(nln,nln);

for i=1:nln
    for j=1:nln
        for k=1:length(w_2D)
            Binv=inv(BJ(:,:,k));                       % inverse
            Jdet=det(BJ(:,:,k));                       % determinant 
            ADV_loc(i,j)=ADV_loc(i,j) + (Jdet.*w_2D(k)) .*  dphiq(1,k,i) *( (beta)*(Grad(k,:,j) * Binv )');
        end
    end
end



                                              
                                              

