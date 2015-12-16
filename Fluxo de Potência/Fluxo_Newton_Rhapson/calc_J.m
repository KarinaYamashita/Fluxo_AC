function [J]=calc_J(Num_Barras,G,B,V,theta,Slack,PV,Q,P)
%---------------------------------------------------------------
% Calcula H
%---------------------------------------------------------------
H=zeros(Num_Barras,Num_Barras);  
for k = 1:Num_Barras
        for m = 1:Num_Barras%Tind
            if k==m
                H(k,k) = -B(k,k)*V(k)^2 - Q(k);
            elseif k<m
                H(k,m) = +V(k)*V(m) * ( G(k,m)*sin(theta(k)-theta(m)) - B(k,m)*cos(theta(k)-theta(m)) );
                H(m,k) = -V(k)*V(m) * ( G(k,m)*sin(theta(k)-theta(m)) + B(k,m)*cos(theta(k)-theta(m)) );
            end
        end
end
      H(Slack,Slack)=1e10;
%---------------------------------------------------------------
% Calcula N
%---------------------------------------------------------------
N=zeros(Num_Barras,Num_Barras);
for k = 1:Num_Barras%Pind
        for m = 1:Num_Barras%Vind
            if k==m 
                N(k,k) = G(k,k)*V(k) + P(k)/V(k)  ;
            elseif k<m
                N(k,m) = V(k) * ( G(k,m)*cos(theta(k)-theta(m)) + B(k,m)*sin(theta(k)-theta(m)) );
                N(m,k) = V(m) * ( G(k,m)*cos(theta(k)-theta(m)) - B(k,m)*sin(theta(k)-theta(m)) );
            end 
        end
end
%---------------------------------------------------------------
% Calcula M
%---------------------------------------------------------------
 M=zeros(Num_Barras,Num_Barras);   
for k = 1:Num_Barras%Qind
        for m = 1:Num_Barras%Tind
            if k==m 
                M(k,k) = -G(k,k)*V(k)^2 + P(k);
            elseif k<m
                M(k,m) = -V(k)*V(m) * ( G(k,m)*cos(theta(k)-theta(m)) + B(k,m)*sin(theta(k)-theta(m)) );
                M(m,k) = -V(m)*V(k) * ( G(k,m)*cos(theta(k)-theta(m)) - B(k,m)*sin(theta(k)-theta(m)) );
            end
        end
 end
%---------------------------------------------------------------
% Calcula L
%---------------------------------------------------------------
L=zeros(Num_Barras,Num_Barras);   
for k = 1:Num_Barras%Qind
        for m = 1:Num_Barras%Vind
            if k==m 
                L(k,k) = -B(k,k)*V(k) + Q(k)/V(k)  ;
            elseif k<m
                L(k,m) = +V(k) * ( G(k,m)*sin(theta(k)-theta(m)) - B(k,m)*cos(theta(k)-theta(m)) );
                L(m,k) = -V(m) * ( G(k,m)*sin(theta(k)-theta(m)) + B(k,m)*cos(theta(k)-theta(m)) );    
            end
        end
end
	 L(Slack,Slack)=1e10;
     for i=1:length(PV)
         L(PV(i),PV(i))=1e10;
     end
 %---------------------------------------------------------------
% Calcula monta a matriz jacobiana
%---------------------------------------------------------------   
 J = [ H N;
       M L ];
    J=inv(J);
end
