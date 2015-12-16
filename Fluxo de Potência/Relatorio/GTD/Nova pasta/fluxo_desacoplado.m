function fluxo_desacoplado
clc
clear all

% Número de barras
    NB = 4;

% Ligações e impedância entre linhas
    z = [1 2 5.33232/529 26.6616/529 9.6883e-5;
         1 3 3.93576/529 19.6788/529 7.3253e-5;
         2 4 3.93576/529 19.6788/529 7.3253e-5;
         3 4 6.72888/529 33.6444/529 0.00012051];
% Classificacão de barras por tipo
    PV = [4];
    PQ = [2 3];
    VT = [1];
    
% Tensões nas barras
    V  = [1 1 1 1.02].';
    T  = [0 0 0 0].';
    
% Potências nas barras
    Pe = [-50 -170 -200 -238].'/100;
    Qe = Pe * tan(acos(0.85));
    
% Erro
    erro = 0.01;


G = zeros(NB);
B = zeros(NB);
P  = Pe;
Q  = Qe;
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
        calc_L;
        T(Tind) = T(Tind) + H \ (Pe(Pind)-P(Pind));
        V(Vind) = V(Vind) + L \ (Qe(Qind)-Q(Qind));
        calc_P;
        calc_Q;
    end


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
                H(k,m) = -V(k)*V(m) * ( G(k,m)*sin(T(k)-T(m)) + B(k,m)*cos(T(k)-T(m)) );
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
                L(k,m) = -V(m) * ( G(k,m)*sin(T(k)-T(m)) + B(k,m)*cos(T(k)-T(m)) );        
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
         YB(z(i,1),z(i,2)) = ((z(i,3) + 1i*z(i,4)))^-1;
         YB(z(i,2),z(i,1)) = ((z(i,3) + 1i*z(i,4)))^-1;
     end
     for i = 1:NB
         YB(i,i) = -sum(YB(:,i)) - 1i* (sum(z(z(:,1)==i,5)) + sum(z(z(:,2)==i,5)));
     end
     YB=-YB;
     G = real(YB);
     B = imag(YB);
end




% Pind = [2 3 4];
% Qind = [2 3];
% 
% Vind = [2 3];
% Tind = [2 3 4];

% con = [ 1 1 1 0 ;
%         1 1 0 1 ;
%         1 0 1 1 ;
%         0 1 1 1 ];
% 
% G = [ +8.951900 -3.815630 -5.169561 0         ;
%       -3.815630 +8.985130 0         -5.169561 ;
%       -5.169561 0         +8.163267 -3.023705 ;
%       0         -5.169561 -3.023705 +8.193264 ];
%    
% B = [ -44.83595 +19.07814 +25.84780 0         ;
%       +19.07814 +44.83395 0         +25.84809 ;
%       +25.84780 0         -40.86384 +15.11852 ;
%       0         +25.84809 +15.11853 -40.86384 ];



end