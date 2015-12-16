function fluxo_DC
clear all
file = 'IEEE_014.txt';

% Número de barras
    NB = textread(file,'%d',1,'whitespace','numero_de_barras :    ','headerlines',0);

% Ligações e impedância entre linhas
    [a,b,c,d,e] = textread(file','%n%n%n%n%n','headerlines',NB+7);
    %c = zeros(size(b));
    z=[a,b,c,d,e];
    
% Classificacão de barras por tipo
    [~,Tipo,V,Pg,~,Pd,~,h] = textread(file,'%n%n%n%n%n%n%n%n',NB,'headerlines',6);
    PQ = find(Tipo==0);
    PV = find(Tipo==2);
    VT = find(Tipo==3);
    
% Tensões nas barras
	%V = ones(NB,1);
    V(V==0) = 1;
    T  = zeros(NB,1);
    
% Potências nas barras
    Sb = 100;
    Pe = (Pg-Pd)/Sb;

G = zeros(NB);
B = zeros(NB);
g = zeros(NB);
b = zeros(NB);
BSH = zeros(NB);
x = zeros(NB);
P = Pe;
Pind = sort([PQ;PV])';
Qind = sort(PQ)';
Vind = sort(PQ)';
Tind = sort([PQ;PV])';

% Subsistema 1
    calc_Y; 
    calc_H;
    tic
	calc_P;
	T(Tind) = T(Tind) + H \ (Pe(Pind)-P(Pind));

% Subsistema 2
    Pind = 1:NB;
    calc_P;

    disp('                V         T         Pg        Pd        ')
    disp([[1:NB]' V  T*180/pi  Pd+P*Sb Pd ])
    calc_fluxo;
    %arr2tab([[1:NB]' V  T*180/pi  Pd+P*Sb  Pd ])

function calc_H
    for k = 1:NB%Pind
        for m = 1:NB%Tind
            if k==m
                H(k,k) = sum(1./x(k,:));
            elseif k<m
                H(k,m) = -x(k,m)^-1;
            elseif k>m
                H(k,m) = -x(k,m)^-1;
            end
        end
    end
	H = H(Pind,Tind);
end

function calc_P
    for k = Pind
        P(k) = 0;
        for m = find(x(k,:)~=0)
           P(k) = P(k) +1/x(k,m)*(T(k)-T(m));
        end
    end
end

function calc_Y
     for i = 1:size(z,1)
         x(z(i,1),z(i,2))  = z(i,4);
         x(z(i,2),z(i,1))  = z(i,4);
     end
     x(x==0) = 1e100;
end

function calc_fluxo
    hl = '\\\hline';
    fprintf('De  Para    Pkm       Pmk     Perdas\n')
    for i = 1:length(a)
        k = z(i,1);
        m = z(i,2);
        Pkm =  + (T(k)-T(m))/x(k,m);
        Pmk =  - (T(k)-T(m))/x(k,m);
        
        %Pkm = +V(k)^2 * g(k,m) - V(k)* V(m) * ( g(k,m)*cos(T(k)-T(m)) + b(k,m)*sin(T(k)-T(m)) );
        %Pmk = +V(m)^2 * g(k,m) - V(k)* V(m) * ( g(k,m)*cos(T(k)-T(m)) - b(k,m)*sin(T(k)-T(m)) );
        %fprintf('%2d   %2d    %+4.3f    %+4.3f   %+4.3f \n',z(i,1),z(i,2),Pkm*Sb,Pmk*Sb,Pkm+Pmk)
        %fprintf('%2d & %2d & %+4.3f & %+4.3f & %+4.3f & %+4.3f & %+4.3f & %+4.3f %s \n',z(i,1),z(i,2),Pkm*Sb,Pmk*Sb,(Pkm+Pmk)*Sb,hl)
        fprintf('%2d & %2d & %+4.3f & %+4.3f & %+4.3f %s \n',k,m,Pkm*Sb,Pmk*Sb,(Pkm+Pmk)*Sb,hl)
    end
end

end