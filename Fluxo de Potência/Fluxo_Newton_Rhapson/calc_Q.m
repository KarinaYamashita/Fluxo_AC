%---------------------------------------------------------------
% Calcula_Potencia_reativa
%---------------------------------------------------------------
 function [Q]=calc_Q(Num_Barras,G,B,V,theta)
Q=zeros(Num_Barras,1);
    for k = 1:Num_Barras
        for m = find(B(k,:)~=0)
           Q(k) = Q(k) + V(m) * ( G(k,m)*sin(theta(k)-theta(m)) - B(k,m)*cos(theta(k)-theta(m)) );  
        end
        Q(k) = V(k) * Q(k);
    end
 end