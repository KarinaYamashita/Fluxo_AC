%---------------------------------------------------------------
 %Calcula_Potencia_ativa
%---------------------------------------------------------------
function [P]= calc_P(Num_Barras,G,B,V,theta)
P=zeros(Num_Barras,1);
la=zeros(Num_Barras,Num_Barras);

    for k =1:Num_Barras
     for m = find(B(k,:)~=0)
            la(k,m)=G(k,m);
           P(k) = P(k) + V(m) * ( G(k,m)*cos(theta(k)-theta(m)) + B(k,m)*sin(theta(k)-theta(m)) )  ;
     end
        P(k) = V(k) * P(k);
    end
end