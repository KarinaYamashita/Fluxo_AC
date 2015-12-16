function fluxo_DC
clear all
file = 'IEEE_014.txt';

% Erro
    tol = 0.0001;
    
% Número de barras
    NB = textread(file,'%d',1,'whitespace','numero_de_barras :    ','headerlines',0);

% Ligações e impedância entre linhas
    [a,b,~,d,e] = textread(file','%n%n%n%n%n','headerlines',NB+7);c = zeros(size(b));
    z=[a,b,c,d,e];
    
% Classificacão de barras por tipo
    [~,Tipo,V,Pg,Qg,Pd,Qd,h] = textread(file,'%n%n%n%n%n%n%n%n',NB,'headerlines',6);
    PQ = find(Tipo==0);
    PV = find(Tipo==2);
    VT = find(Tipo==3);
    
% Tensões nas barras
	V = ones(NB,1);
    T  = zeros(NB,1);
    
% Potências nas barras
    Sb = 100;
    Pe = (Pg-Pd)/Sb;
    Qe = (Qg-Qd)/Sb;
    
G = zeros(NB);
B = zeros(NB);
P = Pe;
it   = 0;
erro = 100;
Pind = sort([PQ;PV])';
Qind = sort(PQ)';
Vind = sort(PQ)';
Tind = sort([PQ;PV])';
Pit = [];
Vit = [];
Tit = [];
DPit = [];
DQit = [];

% Subsistema 1
    calc_Y; 
    calc_H;
    tic
    while erro>tol
        calc_P;
        DPit = [ DPit  Pe-P ];
        T(Tind) = T(Tind) + H \ (Pe(Pind)-P(Pind));
        erro = max(abs(P(Pind)-Pe(Pind)));
        it = it+1;
    end
    toc

% Subsistema 2
    Pind = VT';
    calc_P;

    disp('     it')
    disp(it)
%     disp('      V         T         P         Q')
%     disp([V  T  P  Q])
%     [~,ind] = max([abs(DQit(:,end));abs(DPit(:,end))])
%     plot(DPit(ind,:),DQit(ind,:),'-v')
%     DQit
%     DPit
    calc_fluxo

function calc_H
    for k = Pind
        for m = Tind
            if k==m
                H(k,k) = - B(k,k);
            elseif k<m
                H(k,m) = - B(k,m);
                H(m,k) = - B(k,m);
            end
        end
    end
	H = H(Pind,Tind);
end

function calc_P
    for k = Pind
        P(k) = 0;
        for m = find(B(k,:)~=0)
           P(k) = P(k) + B(k,m)*sin(T(k)-T(m));  
        end
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

function calc_fluxo
    hl = '\\\hline';
    fprintf('De  Para    Pkm       Pmk     Perdas\n')
    for i = 1:length(a)
        k = a(i);
        m = b(i);
        Pkm =  - B(k,m)*sin(T(k)-T(m));
        Pmk =  + B(k,m)*sin(T(k)-T(m));
        %Pt = Pt + Pkm + Pmk;
        %Qt = Qt + Qkm + Qmk;
        fprintf('%2d   %2d    %+4.3f    %+4.3f   %+4.3f \n',a(i),b(i),Pkm,Pmk,Pkm+Pmk)
        %fprintf('%2d & %2d & %+4.3f & %+4.3f & %+4.3f & %+4.3f & %+4.3f & %+4.3f %s \n',a(i),b(i),Pkm,Pmk,Qkm,Qmk,Pkm+Pmk,Qkm+Qmk,hl)
    end
end

end

