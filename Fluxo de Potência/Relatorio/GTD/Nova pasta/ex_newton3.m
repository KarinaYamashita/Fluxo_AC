function ex_newton3
clc
clear all

file = 'IEEE_014.txt';
% Número de barras
    NB = textread(file,'%d',1,'whitespace','numero_de_barras :    ','headerlines',0);

% Ligações e impedância entre linhas
    [a,b,c,d,e] = textread(file','%n%n%n%n%n','headerlines',NB+7);
    z=[a,b,c,d,e];


    
% Classificacão de barras por tipo
    [a,Tipo,V,Pd,Qd,Pg,Qg,h] = textread(file,'%n%n%n%n%n%n%n%n',NB,'headerlines',6);
    PQ = find(Tipo==0);
    PV = find(Tipo==2);
    VT = find(Tipo==3);
    
% Tensões nas barras
    V(PQ') = 1;
    T  = zeros(NB,1);
    
% Potências nas barras
    Spu = 100;
    Pe = -(Pd-Pg)/Spu;
    Qe = -(Qd-Qg)/Spu;
    
    
% Erro
    tol = 0.003;


G = zeros(NB);
B = zeros(NB);
P = Pe;
Q = Qe;
it   = 0;
erro = 100;
Pind = sort([PQ;PV])';
Qind = sort(PQ)';
Vind = sort(PQ)';
Tind = sort([PQ;PV])';
Pit = [];
Qit = [];
Vit = [];
Tit = [];

% Subsistema 1
    calc_Y; 
    calc_P;
    calc_Q;
    erro = max([abs(P(Pind)-Pe(Pind));abs(Q(Qind)-Qe(Qind))]);
    while erro>tol
%         Pit = [Pit P];
%         Qit = [Qit Q];
%         Vit = [Vit V];
%         Tit = [Tit T];
        erro = max([abs(P(Pind)-Pe(Pind));abs(Q(Qind)-Qe(Qind))]);
        calc_J;
%         waitforbuttonpress()
        D  = pinv(J) * [Pe(Pind)-P(Pind);Qe(Qind)-Q(Qind)];
        T(Tind) = T(Tind) + D(1:length(Tind));
        V(Vind) = V(Vind) + D(length(Tind)+1:end);
        calc_P;
        calc_Q;
        it = it +1;
    end



% Subsistema 2
    Pind = sort(VT)';
    Qind = sort(PV)';
    calc_P;
    calc_Q;
    disp('    P        Q         V         T')
    disp([P*Spu  Q*Spu  V  T])
% 
%     Pit
%     Qit
%     Vit
%     Tit
% clf
% plot(0:it-1,Vit)
% legend('V1')

function calc_J
    J = [ calc_H calc_N;
          calc_M calc_L ];
end

function H = calc_H
    for k = Pind
        for m = Tind
            if k==m
                H(k,k) = -B(k,k)*V(k)^2 - Q(k);
            elseif k>m
                H(k,m) = +V(k)*V(m) * ( G(k,m)*sin(T(k)-T(m)) - B(k,m)*cos(T(k)-T(m)) );
            elseif k<m
                H(k,m) = -V(k)*V(m) * ( G(m,k)*sin(T(m)-T(k)) + B(m,k)*cos(T(m)-T(k)) );
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
                N(k,m) = V(k) * ( G(m,k)*cos(T(m)-T(k)) - B(m,k)*sin(T(m)-T(k)) );
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
                M(k,m) = -V(k)*V(m) * ( G(m,k)*cos(T(m)-T(k)) - B(m,k)*sin(T(m)-T(k)) );
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
                L(k,m) = -V(k) * ( G(m,k)*sin(T(m)-T(k)) + B(m,k)*cos(T(m)-T(k)) );        
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
           %Q(k) = Q(k) + V(m) * ( G(k,m)*sin(T(k)-T(m)) - B(k,m)*cos(T(k)-T(m)) );  
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


% 
% function H = calc_H
%     for k = 1:length(Pind)
%         for m = 1:length(Tind)
%             if k==m
%                 H(k,k) = -B(Pind(k),Pind(k))*V(Pind(k))^2 - Q(Pind(k));
%             elseif k>m
%                 H(k,m) = +V(Pind(k))*V(Tind(m)) * ( G(Pind(k),Tind(m))*sin(T(Pind(k))-T(Tind(m))) - B(Pind(k),Tind(m))*cos(T(Pind(k))-T(Tind(m))) );
%             elseif k<m
%                H(k,m) = -V(Pind(k))*V(Tind(m)) * ( G(Pind(k),Tind(m))*sin(T(Tind(m))-T(Pind(k))) + B(Pind(k),Tind(m))*cos(T(Tind(m))-T(Pind(k))) );
%               % H(k,m) = -V(Pind(k))*V(Tind(m)) * ( G(Pind(k),Tind(m))*sin(T(Pind(k))-T(Tind(m))) + B(Pind(k),Tind(m))*cos(T(Pind(k))-T(Tind(m))) );
%             end
%         end
%     end
% end
% 
% function N = calc_N
%     for k = 1:length(Pind)
%         for m = 1:length(Vind)
%             if k==m 
%                 N(k,k) = G(Pind(k),Pind(k))*V(Pind(k)) + P(Pind(k))/V(Pind(k))  ;
%             elseif k>m
%                 N(k,m) = V(Pind(k)) * ( G(Pind(k),Vind(m))*cos(T(Pind(k))-T(Vind(m))) + B(Pind(k),Vind(m))*sin(T(Pind(k))-T(Vind(m))) );
%             elseif k<m
%                 N(k,m) = V(Pind(k)) * ( G(Pind(k),Vind(m))*cos(T(Vind(m))-T(Pind(k))) - B(Pind(k),Vind(m))*sin(T(Vind(m))-T(Pind(k))) );
%               % N(k,m) = V(Pind(k)) * ( G(Pind(k),Vind(m))*cos(T(Pind(k))-T(Vind(m))) - B(Pind(k),Vind(m))*sin(T(Pind(k))-T(Vind(m))) );
%             end 
%         end
%     end
% end
% 
% function M = calc_M
%     for k = 1:length(Qind)
%         for m = 1:length(Tind)
%             if k==m 
%                 M(k,k) = -G(Qind(k),Qind(k))*V(Qind(k))^2 + P(Qind(k));
%             elseif k>m
%                 M(k,m) = -V(Qind(k))*V(Tind(m)) * ( G(Qind(k),Tind(m))*cos(T(Qind(k))-T(Tind(m))) + B(Qind(k),Tind(m))*sin(T(Qind(k))-T(Tind(m))) );
%             elseif k<m
%                M(k,m) = -V(Qind(k))*V(Tind(m)) * ( G(Qind(k),Tind(m))*cos(T(Tind(m))-T(Qind(k))) - B(Qind(k),Tind(m))*sin(T(Tind(m))-T(Qind(k))) );
%               % M(k,m) = -V(Qind(k))*V(Tind(m)) * ( G(Qind(k),Tind(m))*cos(T(Qind(k))-T(Tind(m))) - B(Qind(k),Tind(m))*sin(T(Qind(k))-T(Tind(m))) );
%             end
%         end
%     end
% end
% 
% function L = calc_L
%     for k = 1:length(Qind)
%         for m = 1:length(Vind)
%             if k==m 
%                 L(k,k) = -B(Qind(k),Qind(k))*V(Qind(k)) + Q(Qind(k))/V(Qind(k))  ;
%             elseif k>m
%                 L(k,m) = +V(Qind(k)) * ( G(Qind(k),Vind(m))*sin(T(Qind(k))-T(Vind(m))) - B(Qind(k),Vind(m))*cos(T(Qind(k))-T(Vind(m))) );
%             elseif k<m
%                 L(k,m) = -V(Qind(k)) * ( G(Qind(k),Vind(m))*sin(T(Vind(m))-T(Qind(k))) + B(Qind(k),Vind(m))*cos(T(Vind(m))-T(Qind(k))) );   
%               % L(k,m) = -V(Qind(k)) * ( G(Qind(k),Vind(m))*sin(T(Qind(k))-T(Vind(m))) + B(Qind(k),Vind(m))*cos(T(Qind(k))-T(Vind(m))) );
%             end
%         end
%     end
% % 	L = L(Qind,Vind);
% end