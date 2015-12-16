function fluxo_desacoplado_rapido
clear all
file = 'IEEE_014.txt';

% Erro
    tol = 0.003;
    
% Número de barras
    NB = textread(file,'%d',1,'whitespace','numero_de_barras :    ','headerlines',0);

% Ligações e impedância entre linhas
    [a,b,c,d,e] = textread(file','%n%n%n%n%n','headerlines',NB+7);
    z=[a,b,c,d,e];

% Classificacão de barras por tipo
    [~,Tipo,V,Pg,Qg,Pd,Qd,h] = textread(file,'%n%n%n%n%n%n%n%n',NB,'headerlines',6);
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
g = zeros(NB);
b = zeros(NB);
BSH = zeros(NB);
x = zeros(NB);
P = Pe;
Q = Qe;
Pkm = 0;
Qkm = 0;
Pmk = 0;
Qmk = 0;
p = 0;
q = 0;
KP = 1;
KQ = 1;
Pind = [PQ;PV]';
Qind = PQ';
Vind = PQ';
Tind = [PQ;PV]';
Pit = [];
Qit = [];
Vit = [];
Tit = [];
DPit = [];
DQit = [];

% Subsistema 1
    calc_Y; 
    calc_H;
    calc_L;
    tic
    while 1
        [max(abs(P(Pind)-Pe(Pind))) max(abs(Q(Qind)-Qe(Qind))./V(Qind))];
        calc_P;
        if max(abs(P(Pind)-Pe(Pind)))>tol
            T(Tind) = T(Tind) +  inv(H) * (Pe(Pind)-P(Pind))  ./V(Pind);
            p = p + 1;
            KQ = 1;
        else
            KP = 0;
            if KQ == 0
               break 
            end
        end
        calc_Q;
        if max(abs(Q(Qind)-Qe(Qind)))>tol
            V(Vind) = V(Vind) +  inv(L) * (Qe(Qind)-Q(Qind))  ./V(Qind);
            q = q + 1;
            KP = 1;
        else
            KQ = 0;
            if KP == 0
                break
            end
        end
    end
    toc
    erro = max([abs(Q(Qind)-Qe(Qind)) abs(Q(Qind)-Qe(Qind))]);
    
% Subsistema 2
    Pind = VT';
    Qind = [PV;VT]';
    calc_P;
    calc_Q;

    disp('    p     q')
    disp([p q])
    disp('                V         T         Pg        Qg        Pd        Qd')
    disp([[1:NB]' V  T*180/pi  Pd+P*Sb Qd+Q*Sb Pd  Qd])
%     arr2tab([[1:NB]' V  T*180/pi  Pd+P*Sb Qd+Q*Sb Pd  Qd])
%     [~,ind] = max([abs(DQit(:,end));abs(DPit(:,end))])
%     plot(DPit(ind,:),DQit(ind,:),'-v')
%     DQit
%     DPit
    calc_fluxo;

function calc_H
    for k = 1:NB%Pind
        for m = 1:NB%Tind
            if k==m
                H(k,k) = sum(1./x(k,:));
                %H(k,k) = - B(k,k);
            elseif k<m
                H(k,m) = -x(k,m)^-1;
                %H(k,m) = - B(k,m);
            elseif k>m
                H(k,m) = -x(k,m)^-1;
                %H(k,m) = - B(k,m);
            end
        end
    end
	H = H(Pind,Tind);
end

function calc_L
    for k = 1:NB%Qind
        for m = 1:NB%Vind
            if k==m 
                L(k,k) = -B(k,k);
            elseif k<m
                L(k,m) = -B(k,m);
            elseif k>m
                L(k,m) = -B(k,m);
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
         x(z(i,1),z(i,2))  = z(i,4);
         x(z(i,2),z(i,1))  = z(i,4);
         BSH(z(i,1),z(i,2)) = z(i,5);
         BSH(z(i,2),z(i,1)) = z(i,5);
     end
     g = real(YB);
     b = imag(YB);
     for i = 1:NB
         YB(i,i) = -sum(YB(:,i)) - 1i* (sum(z(z(:,1)==i,5)) + sum(z(z(:,2)==i,5)));
     end
     YB=-YB;
     G = real(YB);
     B = imag(YB);
     x(x==0) = 1e100;
end

function calc_fluxo
    hl = '\\\hline';
    fprintf('De  Para     Pkm       Pmk       Qkm       Qmk           Perdas\n')
    for i = 1:length(z(:,1))
        k = z(i,1);
        m = z(i,2);
        Pkm = +V(k)^2 * g(k,m) - V(k)* V(m) * ( g(k,m)*cos(T(k)-T(m)) + b(k,m)*sin(T(k)-T(m)) );
        Pmk = +V(m)^2 * g(k,m) - V(k)* V(m) * ( g(k,m)*cos(T(k)-T(m)) - b(k,m)*sin(T(k)-T(m)) );
        Qkm = -V(k)^2 * (  b(k,m) + BSH(k,m) ) - V(k)*V(m) * ( g(k,m)*sin(T(k)-T(m)) - b(k,m)*cos(T(k)-T(m)) );
        Qmk = -V(m)^2 * (  b(k,m) + BSH(k,m) ) + V(k)*V(m) * ( g(k,m)*sin(T(k)-T(m)) + b(k,m)*cos(T(k)-T(m)) );
        %fprintf('%2d   %2d    %+4.3f    %+4.3f    %+4.3f    %+4.3f     %+4.3f    %+4.3f \n',a(i),b(i),Pkm,Qkm,Pmk,Qmk,Pkm+Pmk,Qkm+Qmk)
        fprintf('%2d & %2d & %+4.3f & %+4.3f & %+4.3f & %+4.3f & %+4.3f & %+4.3f %s \n',z(i,1),z(i,2),Pkm*Sb,Pmk*Sb,Qkm*Sb,Qmk*Sb,(Pkm+Pmk)*Sb,(Qkm+Qmk)*Sb,hl)
    end
end

end
