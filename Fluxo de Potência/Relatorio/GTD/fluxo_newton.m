function fluxo_newton
clear all
file = 'IEEE_014.txt';

% Erro
    tol = 0.001;
    
% Número de barras
    NB = textread(file,'%d',1,'whitespace','numero_de_barras :    ','headerlines',0);

% Ligações e impedância entre linhas
    [a,b,c,d,e] = textread(file','%n%n%n%n%n','headerlines',NB+7);
    z=[a,b,c,d,e];

% Classificacão de barras por tipo
    [~,Tipo,V,Pg,Qg,Pd,Qd,~] = textread(file,'%n%n%n%n%n%n%n%n',NB,'headerlines',6);
    PQ = find(Tipo==0);
    PV = find(Tipo==2);
    VT = find(Tipo==3);
    
% Tensões nas barras
    V(PQ') = 1;
    T  = zeros(NB,1);
    
% Potências nas barras
    Sb = 100;
    Pe = (Pg-Pd)/Sb;
    Qe = (Qg-Qd)/Sb;
   
G = zeros(NB);
B = zeros(NB);
BSH = zeros(NB);
P = Pe;
Q = Qe;
Pkm = 0;
Qkm = 0;
Pmk = 0;
Qmk = 0;
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
DPit = [];
DQit = [];

% Subsistema 1
    calc_Y; 
    calc_P;
    calc_Q;
    erro = max([abs(P(Pind)-Pe(Pind));abs(Q(Qind)-Qe(Qind))]);
    tic
    while erro>tol
        calc_J;
        D  = J \ [Pe(Pind)-P(Pind);Qe(Qind)-Q(Qind)];
        T(Tind) = T(Tind) + D(1:length(Tind));
        V(Vind) = V(Vind) + D(length(Tind)+1:end);
        calc_P;
        calc_Q;
        it = it +1;
        erro = max([abs(P(Pind)-Pe(Pind));abs(Q(Qind)-Qe(Qind))]);
        DPit = [ DPit  Pe-P ];
        DQit = [ DQit  Qe-Q ];
    end
    toc
    
% Subsistema 2
    Pind = VT';
    Qind = [PV;VT]';
    calc_P;
    calc_Q;
    
    disp('     it')
    disp(it)
    disp('      V         T         P         Q')
    disp([V  T  P  Q])
    arr2tab([[1:NB]' V  T  P  Q])
%     [~,ind] = max([abs(DQit(:,end));abs(DPit(:,end))])
%     plot(DPit(ind,:),DQit(ind,:),'-v')
%     DQit
%     DPit
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
     for i = 1:NB
         YB(i,i) = -sum(YB(:,i)) - 1i* (sum(z(z(:,1)==i,5)) + sum(z(z(:,2)==i,5)));
     end
     YB=-YB;
     G = real(YB);
     B = imag(YB);
     for i = 1:length(a)
        BSH(a(i),b(i)) = e(i);
        BSH(b(i),a(i)) = e(i);
     end
end

function calc_fluxo
    hl = '\\\hline';
    fprintf('De  Para     Pkm       Qkm       Pmk       Qmk           Perdas\n')
    for i = 1:length(a)
        k = a(i);
        m = b(i);
        Pkm = +V(k)^2 * G(k,m) - V(k)*V(m)*G(k,m)*cos(T(k)-T(m)) - V(k)*V(m)*B(k,m)*sin(T(k)-T(m));
        Pmk = +V(m)^2 * G(k,m) - V(k)*V(m)*G(k,m)*cos(T(k)-T(m)) + V(k)*V(m)*B(k,m)*sin(T(k)-T(m));
        Qkm = -V(k)^2 * B(k,m) - V(k)^2 * BSH(k,m) + V(k)*V(m)*B(k,m)*cos(T(k)-T(m)) - V(k)*V(m)*G(k,m)*sin(T(k)-T(m));
        Qmk = -V(m)^2 * B(k,m) - V(k)^2 * BSH(k,m)+ V(k)*V(m)*B(k,m)*cos(T(k)-T(m)) + V(k)*V(m)*G(k,m)*sin(T(k)-T(m));
        %Pt = Pt + Pkm + Pmk;
        %Qt = Qt + Qkm + Qmk;
        %fprintf('%2d   %2d    %+4.3f    %+4.3f    %+4.3f    %+4.3f     %+4.3f    %+4.3f \n',a(i),b(i),Pkm,Qkm,Pmk,Qmk,Pkm+Pmk,Qkm+Qmk)
        fprintf('%2d & %2d & %+4.3f & %+4.3f & %+4.3f & %+4.3f & %+4.3f & %+4.3f %s \n',a(i),b(i),Pkm,Pmk,Qkm,Qmk,Pkm+Pmk,Qkm+Qmk,hl)
    end
end

end

% function H = calc_H
%     for k = Pind
%         for m = Tind
%             if k==m
%                 H(k,k) = -B(k,k)*V(k)^2 - Q(k);
%             elseif k>m
%                 H(k,m) = +V(k)*V(m) * ( G(k,m)*sin(T(k)-T(m)) - B(k,m)*cos(T(k)-T(m)) );
%             elseif k<m
%                 H(k,m) = -V(k)*V(m) * ( G(m,k)*sin(T(m)-T(k)) + B(m,k)*cos(T(m)-T(k)) );
%             end
%         end
%     end
% 	H = H(Pind,Tind);
% end
% 
% function N = calc_N
%     for k = Pind
%         for m = Vind
%             if k==m 
%                 N(k,k) = G(k,k)*V(k) + P(k)/V(k)  ;
%             elseif k>m
%                 N(k,m) = V(k) * ( G(k,m)*cos(T(k)-T(m)) + B(k,m)*sin(T(k)-T(m)) );
%             elseif k<m
%                 N(k,m) = V(k) * ( G(m,k)*cos(T(m)-T(k)) - B(m,k)*sin(T(m)-T(k)) );
%             end 
%         end
%     end
%     N = N(Pind,Vind);
% 
% end
% 
% function M = calc_M
%     for k = Qind
%         for m = Tind
%             if k==m 
%                 M(k,k) = -G(k,k)*V(k)^2 + P(k);
%             elseif k>m
%                 M(k,m) = -V(k)*V(m) * ( G(k,m)*cos(T(k)-T(m)) + B(k,m)*sin(T(k)-T(m)) );
%             elseif k<m
%                 M(k,m) = -V(k)*V(m) * ( G(m,k)*cos(T(m)-T(k)) - B(m,k)*sin(T(m)-T(k)) );
%             end
%         end
%     end
% 	M = M(Qind,Tind);
% end
% 
% function L = calc_L
%     for k = Qind
%         for m = Vind
%             if k==m 
%                 L(k,k) = -B(k,k)*V(k) + Q(k)/V(k)  ;
%             elseif k>m
%                 L(k,m) = +V(k) * ( G(k,m)*sin(T(k)-T(m)) - B(k,m)*cos(T(k)-T(m)) );
%             elseif k<m
%                 L(k,m) = -V(k) * ( G(m,k)*sin(T(m)-T(k)) + B(m,k)*cos(T(m)-T(k)) );        
%             end
%         end
%     end
% 	L = L(Qind,Vind);
% end


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