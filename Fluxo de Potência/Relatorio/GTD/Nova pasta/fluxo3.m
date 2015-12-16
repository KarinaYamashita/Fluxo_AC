function fluxo3
clc
clear all

con = [ 1 1
        1 1 ];

G = [ +0.1923 -0.1923 ;
      -0.1923 +0.1923 ];
   
B = [ -0.9415 +0.9615 ;
      +0.9615 -0.9415 ];

Pind = [2];
Qind = [2];

Vind = [2];
Tind = [2];

V  = [1 1].';
T  = [0 0].';

Pe = [1 -.3].';
Qe = [1 .07].';
P  = Pe;
Q  = Qe;

erro = 0.003;
it   = 1;

calc_P;
calc_Q;
calc_J;
while max(abs([Pe-P;Qe-Q]))>erro
    it = it +1;
    calc_J;
    D  = J \ [Pe(Pind)-P(Pind);Qe(Qind)-Q(Qind)];
    T(Tind) = T(Tind) + D(1:length(Tind));
    V(Vind) = V(Vind) + D(length(Tind)+1:end);
    calc_P;
    calc_Q;
    fprintf('Iteração %d ---> Erro máximo = %f\n',it,max(abs([P-Pe;Q-Qe])) )
%     waitforbuttonpress
end
disp('    Pe        P         Q         Qe         V         T')
disp([Pe  P  Qe  Q  V  T])




function H = calc_H
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

function N = calc_N
    for k = Pind
        for m = Vind
            if k==m 
                N(k,k) = G(k,k)*V(k) + P(k)/V(k)  ;
            elseif k>m
                N(k,m) = V(k) * ( G(k,m)*cos(T(k)-T(m)) + B(k,m)*sin(T(k)-T(m)) );
            elseif k<m
                N(k,m) = V(k) * ( G(k,m)*cos(T(k)-T(m)) - B(k,m)*sin(T(k)-T(m)) );
            end 
        end
    end
    N = N(Pind,Vind);
end

function M = calc_M
    for k = Qind
        for m = Tind
            if k==m 
                M(k,k) = -G(k,k)*V(k)^2 + P(k);
            elseif k>m
                M(k,m) = -V(k)*V(m) * ( G(k,m)*cos(T(k)-T(m)) + B(k,m)*sin(T(k)-T(m)) );
            elseif k<m
                M(k,m) = -V(k)*V(m) * ( G(k,m)*cos(T(k)-T(m)) - B(k,m)*sin(T(k)-T(m)) );
            end
        end
    end
	M = M(Qind,Tind);
end

function L = calc_L
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
        for m = find(con(k,:)~=0)
           P(k) = P(k) + V(m) * ( G(k,m)*cos(T(k)-T(m)) + B(k,m)*sin(T(k)-T(m)) );  
        end
        P(k) = V(k) * P(k);
    end
end

function calc_Q
    for k = Qind
        Q(k) = 0;
        for m = find(con(k,:)~=0)
           Q(k) = Q(k) + V(m) * ( G(k,m)*sin(T(k)-T(m)) - B(k,m)*cos(T(k)-T(m)) );  
        end
        Q(k) = V(k) * Q(k);
    end
end

function calc_J
    J = [ calc_H calc_N;
          calc_M calc_L ];
end

end