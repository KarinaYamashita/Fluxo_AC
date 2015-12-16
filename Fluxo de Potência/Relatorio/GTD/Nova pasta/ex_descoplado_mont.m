function ex_descoplado_mont
clc
clear all

% Número de barras
    NB = 2;
    
% Ligações e impedância entre linhas
    z = [1 2 0.2 1 0.02];
    
% Classificacão de barras por tipo
    PV = [];
    PQ = [2];
    VT = [1];
    
% Tensões nas barras
    V  = [1 1].';
    T  = [0 0].';
    
% Potências nas barras
    Pe = [1 -.3].';
    Qe = [1 0.07].';
    
% Erro
    erro = 0.003;


G = zeros(NB);
B = zeros(NB);
P = zeros(length(Pe),1);
Q = zeros(length(Pe),1);
it   = 0;
Pind = sort([PQ PV]);
Qind = sort(PQ);
Vind = sort(PQ);
Tind = sort([PQ PV]);


% Subsistema 1
    calc_Y;
    calc_P;
    calc_Q;
    while max(abs([Pe(Pind)-P(Pind);Qe(Pind)-Q(Pind)]))>erro

        fprintf('Iteração %d ---> Erro máximo = %f\n',it,max(abs([Pe(Pind)-P(Pind);Qe(Pind)-Q(Pind)])) )
        disp('    V         T')
        disp([V  T])
        it = it +1;
        calc_H;
        T(Tind) = T(Tind) + H \ (Pe(Pind)-P(Pind));
        calc_Q;
        calc_L;
        V(Vind) = V(Vind) + L \ (Qe(Qind)-Q(Qind));
        calc_P;
        calc_Q;
    end
        fprintf('Iteração %d ---> Erro máximo = %f\n',it,max(abs([Pe(Pind)-P(Pind);Qe(Pind)-Q(Pind)])) )
        disp('    V         T')
        disp([V  T])

% Subsistema 2
    Pind = 1:NB;
    Qind = 1:NB;
    Vind = 1:NB;
    Tind = 1:NB;
    calc_P;
    calc_Q;
    disp('    P        Q         V         T')
    disp([P  Q  V  T])


function calc_H
    for k = Pind
        for m = Tind
            if k==m
                H(k,k) = -B(k,k)*V(k)^2 - Q(k);
            elseif k>m
                H(k,m) = +V(k)*V(m) * ( G(k,m)*sin(T(k)-T(m)) - B(k,m)*cos(T(k)-T(m)) );
            elseif k<m
                H(k,m) = -V(k)*V(m) * ( G(k,m)*sin(T(m)-T(k)) + B(k,m)*cos(T(m)-T(k)) );
            end
        end
    end
	H = H(Pind,Tind);
end

function calc_L
    for k = Qind
        for m = Vind
            if k==m 
                L(k,k) = -B(k,k)*V(k) + Q(k)/V(k)  ;
            elseif k>m
                L(k,m) = +V(k) * ( G(k,m)*sin(T(k)-T(m)) - B(k,m)*cos(T(k)-T(m)) );
            elseif k<m
                L(k,m) = -V(k) * ( G(k,m)*sin(T(m)-T(k)) + B(k,m)*cos(T(m)-T(k)) );        
            end
        end
    end
	L = L(Qind,Vind);
end

function calc_P
    for k = Pind
        P(k) = 0;
        for m = find(G(k,:)~=0)
           P(k) = P(k) + V(m) * ( G(k,m)*cos(T(k)-T(m)) + B(k,m)*sin(T(k)-T(m)) );  
        end
        P(k) = V(k) * P(k);
    end
end

function calc_Q
    for k = Qind
        Q(k) = 0;
        for m = find(G(k,:)~=0)
           Q(k) = Q(k) + V(m) * ( G(k,m)*sin(T(k)-T(m)) - B(k,m)*cos(T(k)-T(m)) );  
        end
        Q(k) = V(k) * Q(k);
    end
end

function calc_Y
     YB = zeros(NB,NB);
     for i = 1:size(z,1)
         YB(z(i,1),z(i,2)) = (z(i,3) + 1i*z(i,4))^-1;
         YB(z(i,2),z(i,1)) = (z(i,3) + 1i*z(i,4))^-1;
     end
     for i = 1:NB
         YB(i,i) = -sum(YB(:,i)) - 1i* (sum(z(z(:,1)==i,5)) + sum(z(z(:,2)==i,5)));
     end
     YB=-YB;
     G = real(YB);
     B = imag(YB);
end





end