function ex_newton2_mont
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
        calc_J;
        J
        D  = J \ [Pe(Pind)-P(Pind);Qe(Qind)-Q(Qind)];
        T(Tind) = T(Tind) + D(1:length(Tind));
        V(Vind) = V(Vind) + D(length(Tind)+1:end);
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


function calc_J
    J = [ calc_H calc_N;
          calc_M calc_L ];
end

function H = calc_H
    for k = 1:length(Pind)
        for m = 1:length(Tind)
            if k==m
                H(k,k) = -B(Pind(k),Pind(k))*V(Pind(k))^2 - Q(Pind(k));
            elseif k>m
                H(k,m) = +V(Pind(k))*V(Tind(m)) * ( G(Pind(k),Tind(m))*sin(T(Pind(k))-T(Tind(m))) - B(Pind(k),Tind(m))*cos(T(Pind(k))-T(Tind(m))) );
            elseif k<m
                H(k,m) = -V(Pind(k))*V(Tind(m)) * ( G(Pind(k),Tind(m))*sin(T(Tind(m))-T(Pind(k))) + B(Pind(k),Tind(m))*cos(T(Tind(m))-T(Pind(k))) );
            end
        end
    end
end

function N = calc_N
    for k = 1:length(Pind)
        for m = 1:length(Vind)
            if k==m 
                N(k,k) = G(Pind(k),Pind(k))*V(Pind(k)) + P(Pind(k))/V(Pind(k))  ;
            elseif k>m
                N(k,m) = V(Pind(k)) * ( G(Pind(k),Vind(m))*cos(T(Pind(k))-T(Vind(m))) + B(Pind(k),Vind(m))*sin(T(Pind(k))-T(Vind(m))) );
            elseif k<m
                N(k,m) = V(Pind(k)) * ( G(Pind(k),Vind(m))*cos(T(Vind(m))-T(Pind(k))) - B(Pind(k),Vind(m))*sin(T(Vind(m))-T(Pind(k))) );
            end 
        end
    end
end

function M = calc_M
    for k = 1:length(Qind)
        for m = 1:length(Tind)
            if k==m 
                M(k,k) = -G(Qind(k),Qind(k))*V(Qind(k))^2 + P(Qind(k));
            elseif k>m
                M(k,m) = -V(Qind(k))*V(Tind(m)) * ( G(Qind(k),Tind(m))*cos(T(Qind(k))-T(Tind(m))) + B(Qind(k),Tind(m))*sin(T(Qind(k))-T(Tind(m))) );
            elseif k<m
                M(k,m) = -V(Qind(k))*V(Tind(m)) * ( G(Qind(k),Tind(m))*cos(T(Tind(m))-T(Qind(k))) - B(Qind(k),Tind(m))*sin(T(Tind(m))-T(Qind(k))) );
            end
        end
    end
end

function L = calc_L
    for k = 1:length(Qind)
        for m = 1:length(Vind)
            if k==m 
                L(k,k) = -B(Qind(k),Qind(k))*V(Qind(k)) + Q(Qind(k))/V(Qind(k))  ;
            elseif k>m
                L(k,m) = +V(Qind(k)) * ( G(Qind(k),Vind(m))*sin(T(Qind(k))-T(Vind(m))) - B(Qind(k),Vind(m))*cos(T(Qind(k))-T(Vind(m))) );
            elseif k<m
                L(k,m) = -V(Qind(k)) * ( G(Qind(k),Vind(m))*sin(T(Vind(m))-T(Qind(k))) + B(Qind(k),Vind(m))*cos(T(Vind(m))-T(Qind(k))) );        
            end
        end
    end
% 	L = L(Qind,Vind);
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