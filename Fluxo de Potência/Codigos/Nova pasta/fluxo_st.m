function fluxo_st
clc
clear all

% Número de barras
    NB = 4;

% Ligações e impedância entre linhas
    z = [1 2 5.33232/529 26.6616/529 9.6883e-5;
         1 3 3.93576/529 19.6788/529 7.3253e-5;
         2 4 3.93576/529 19.6788/529 7.3253e-5;
         3 4 6.72888/529 33.6444/529 0.00012051];
     
    z = [1 2 0.01008 0.05040 0.05125;
         1 3 0.00744 0.03720 0.03875;
         2 4 0.00744 0.03720 0.03875;
         3 4 0.01272 0.06360 0.06375 ];

     
     
% Classificacão de barras por tipo
    PV = [4];
    PQ = [2 3];
    VT = [1];
    
% Tensões nas barras
    V  = [1 1 1 1.02].';
    T  = [0 0 0 0].';
    
% Potências nas barras
    Pe = [-50 -170 -200 +238].'/100;
%     Qe = Pe * tan(acos(0.85));
    Qe = [-30.99 -105.35 -123.94 -49.58]'/100;
    
% Erro
    erro = 0.001;

G = zeros(NB);
B = zeros(NB);
g = zeros(NB);
b = zeros(NB);
BSH = zeros(NB);
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
%         disp('    V         T')
%         disp([V  T])
        it = it +1;
        calc_J;
        D  = J \ [Pe(Pind)-P(Pind);Qe(Qind)-Q(Qind)];
        T(Tind) = T(Tind) + D(1:length(Tind));
        V(Vind) = V(Vind) + D(length(Tind)+1:end);
        calc_P;
        calc_Q;
    end


% Subsistema 2
    Pind = sort(VT)';
    Qind = sort([VT;PV])';
    Pind = 1:NB;
    Qind = 1:NB;
    calc_P;
    calc_Q;
    it
    disp([ [P Q]*100 V  T*180/pi])

    
    calc_fluxo



function calc_J
    J = [ calc_H calc_N;
          calc_M calc_L ];
end

function H = calc_H
    for k = 1:NB%Pind
        for m = 1:NB%Tind
            if k==m
                H(k,k) = -B(k,k)*V(k)^2 - Q(k);
            elseif k<m
                H(k,m) = +V(k)*V(m) * ( G(k,m)*sin(T(k)-T(m)) - B(k,m)*cos(T(k)-T(m)) );
                H(m,k) = -V(k)*V(m) * ( G(k,m)*sin(T(k)-T(m)) + B(k,m)*cos(T(k)-T(m)) );
            end
        end
    end
	H = H(Pind,Tind);
end

function N = calc_N
    for k = 1:NB%Pind
        for m = 1:NB%Vind
            if k==m 
                N(k,k) = G(k,k)*V(k) + P(k)/V(k)  ;
            elseif k<m
                N(k,m) = V(k) * ( G(k,m)*cos(T(k)-T(m)) + B(k,m)*sin(T(k)-T(m)) );
                N(m,k) = V(m) * ( G(k,m)*cos(T(k)-T(m)) - B(k,m)*sin(T(k)-T(m)) );
            end 
        end
    end
    N = N(Pind,Vind);
end

function M = calc_M
    for k = 1:NB%Qind
        for m = 1:NB%Tind
            if k==m 
                M(k,k) = -G(k,k)*V(k)^2 + P(k);
            elseif k<m
                M(k,m) = -V(k)*V(m) * ( G(k,m)*cos(T(k)-T(m)) + B(k,m)*sin(T(k)-T(m)) );
                M(m,k) = -V(m)*V(k) * ( G(k,m)*cos(T(k)-T(m)) - B(k,m)*sin(T(k)-T(m)) );
            end
        end
    end
	M = M(Qind,Tind);
end

function L = calc_L
    for k = 1:NB%Qind
        for m = 1:NB%Vind
            if k==m 
                L(k,k) = -B(k,k)*V(k) + Q(k)/V(k)  ;
            elseif k<m
                L(k,m) = +V(k) * ( G(k,m)*sin(T(k)-T(m)) - B(k,m)*cos(T(k)-T(m)) );
                L(m,k) = -V(m) * ( G(k,m)*sin(T(k)-T(m)) + B(k,m)*cos(T(k)-T(m)) );    
            end
        end
    end
	L = L(Qind,Vind);
end

function calc_P
    for k = Pind
        P(k) = 0;
        for m = find(B(k,:)~=0)
           P(k) = P(k) + V(m) * ( G(k,m)*cos(T(k)-T(m)) + B(k,m)*sin(T(k)-T(m)) );  
        end
        P(k) = V(k) * P(k);
    end
end

function calc_Q
    for k = Qind
        Q(k) = 0;
        for m = find(B(k,:)~=0)
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
     g = real(YB);
     b = imag(YB);
     for i = 1:NB
         YB(i,i) = -sum(YB(:,i)) - 1i* (sum(z(z(:,1)==i,5)) + sum(z(z(:,2)==i,5)));
     end
     YB=-YB;
     G = real(YB);
     B = imag(YB);
     for i = 1:length(z(:,1))
        BSH(z(i,1),z(i,2)) = z(i,end);
        BSH(z(i,2),z(i,1)) = z(i,end);
     end
end

function calc_fluxo
    hl = '\\\hline';
    fprintf('De  Para     Pkm       Pmk       Qkm       Qmk           Perdas\n')
    for i = 1:length(z(:,1))
        k = z(i,1);
        m = z(i,2);
        Pkm = +V(k)^2 * g(k,m) - V(k)*V(m) * ( g(k,m)*cos(T(k)-T(m)) + b(k,m)*sin(T(k)-T(m)) );
        Qkm = -V(k)^2 * (  b(k,m) + BSH(k,m) ) - V(k)*V(m) * ( g(k,m)*sin(T(k)-T(m)) - b(k,m)*cos(T(k)-T(m)) );

        Pmk = +V(m)^2 * g(k,m) - V(k)*V(m) * ( g(k,m)*cos(T(k)-T(m)) - b(k,m)*sin(T(k)-T(m)) );
        Qmk = -V(m)^2 * (  b(k,m) + BSH(k,m) ) + V(k)*V(m) * ( g(k,m)*sin(T(k)-T(m)) + b(k,m)*cos(T(k)-T(m)) );
        %Pt = Pt + Pkm + Pmk;
        %Qt = Qt + Qkm + Qmk;
        %fprintf('%2d   %2d    %+4.3f    %+4.3f    %+4.3f    %+4.3f     %+4.3f    %+4.3f \n',a(i),b(i),Pkm,Qkm,Pmk,Qmk,Pkm+Pmk,Qkm+Qmk)
        fprintf('%2d & %2d & %+4.3f & %+4.3f & %+4.3f & %+4.3f & %+4.3f & %+4.3f %s \n',z(i,1),z(i,2),Pkm*100,Pmk*100,Qkm*100,Qmk*100,Pkm+Pmk,Qkm+Qmk,hl)
    end
end

end

% Pind = [2  4];
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

