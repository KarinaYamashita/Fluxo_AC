
clc
%-------------------------------------------------------------------    
%Determina o sistema que ser� utilizado para se realizar os c�lculos
%-------------------------------------------------------------------   
    file = 'IEEE_057.txt';
%     file = '02_PQ.txt';

%------------------------------
% Toler�ncia para o erro
%-------------------------------    
    tol = 0.003;

%-------------------------------    
% N�mero de barras
%-------------------------------
    Num_Barras = textread(file,'%d',1,'whitespace','numero_de_barras :    ','headerlines',0);
    Num_linhas = textread(file,'%d',1,'whitespace','numero_de_lineas :    ','headerlines',1);
%------------------------------------
% Liga��es e imped�ncia entre linhas
%------------------------------------
      [N_i,N_f,R,X,Y_shunt] = textread(file','%n%n%n%n%n','headerlines',Num_Barras+7);
      z=[N_i,N_f,R,X,Y_shunt];
      
%----------------------------------
% Classificac�o de barras por tipo
%----------------------------------
     [~,Tipo,V,Pg,Qg,Pd,Qd,~] = textread(file,'%n%n%n%n%n%n%n%n',Num_Barras,'headerlines',6);
    PQ = find(Tipo==0);
    PV = find(Tipo==2);
    Slack = find(Tipo==3);

%------------------------------
% N�mero de cada tipo de barra
%------------------------------
NPQ= length(PQ);
NPV= length(PV);
%---------------------------------------------------------------
% Inicializa��o do valor de m�dulo e fase das tens�es nas barras
%---------------------------------------------------------------
    V(PQ') = 1; %as barras PQ s�o as �nicas barras onde ir� variar o valor de tens�o
    theta = zeros(Num_Barras,1);
    
 %------------------------------    
% Potências nas barras
%------------------------------
     Sb = 100;%para sistemas que não estão em pu
%   Sb = 1;%para sistemas que estão em pu
Pe=zeros(Num_Barras,1);
Qe=zeros(Num_Barras,1);
for i=1:Num_Barras
    Pe(i) = (Pg(i)-Pd(i))/Sb;
    Qe(i) = (Qg(i)-Qd(i))/Sb;
end
%------------------------------    
% Inicializa��o de variaveis
%------------------------------
it   = 0;
Pind = sort([PQ;PV])';
Qind = sort(PQ)';
%---------------------------------------------------------------
% Cálculo do subsistema 1(Calcula V e theta)
%---------------------------------------------------------------
%---------------------------------------------------------------
% Calcula a admit�ncia do sistema
%---------------------------------------------------------------
        Y_Bus= zeros(Num_Barras,Num_Barras);
        Y_Sh= zeros(Num_Barras,1);
        for i = 1:Num_linhas
            Y_Bus(N_i(i),N_f(i)) = 1/(R(i)+ 1j*X(i));
            Y_Bus(N_f(i),N_i(i)) = 1/(R(i)+ 1j*X(i));
            Y_Sh(N_i(i),N_f(i))= Y_shunt(i);
            Y_Sh(N_f(i),N_i(i))= Y_shunt(i);
        end
       
        for i=1:Num_Barras   
            %Y_Bus(i,i)=sum(Y_Bus(i,:)+Y_shunt
            Y_Bus(i,i) = -sum(Y_Bus(i,:)) - 1i*(sum(Y_shunt(N_i(:)==i))+sum(Y_shunt(N_f(:)==i)));
        end
%  G = zeros(Num_Barras);
%  B = zeros(Num_Barras);
Y_Bus=-Y_Bus;
G = real(Y_Bus);
B = imag(Y_Bus);
% P =zeros(Num_Barras,1);
% Q = zeros(Num_Barras,1);
P=calc_P(Num_Barras,G,B,V,theta);
Q=calc_Q(Num_Barras,G,B,V,theta);
 J=calc_J(Num_Barras,G,B,V,theta,Slack,PV,Q,P);

if NPQ~=0
        erro_m=max(abs([Pe-P;Qe-Q]));
    else
        erro_m = max(abs(Pe-P));
end
tic
%   while it<5
while erro_m>0.003
    J=calc_J(Num_Barras,G,B,V,theta,Slack,PV,Q,P);
    Delta=[Pe-P;Qe-Q];
    b=size(Delta);
     D  = J*Delta
    theta = theta + D(1:length(theta));
    V = V+ D(length(theta)+1:end);
    P=calc_P(Num_Barras,G,B,V,theta);
    Q=calc_Q(Num_Barras,G,B,V,theta);
    it = it +1;
    if NPQ~=0
        erro_m=max(abs([Pe(Pind)-P(Pind);Qe(Qind)-Q(Qind)]));
    else
        erro_m = max(abs(Pe(Pind)-P(Pind)));
    end
 end
    toc
    erro_m
    V
    theta
%     J
%  %---------------------------------------------------------------
% % Cálculo do subsistema 2(Calcula V e theta)
% %---------------------------------------------------------------
%     Pind = Slack';
%     Qind = [PV;Slack]';
%     calc_P(Num_Barras,G,B,V,theta);
%     calc_Q(Num_Barras,G,B,V,theta);
 Q=calc_Q(Num_Barras,G,B,V,theta)
     disp('    it')
     disp((it))
     disp('Num_Barras|  V     |    T   |     Pg        Qg        Pd        Qd')
   for i=1:Num_Barras
     fprintf('   %2d     | %+4.3f | %+4.3f | %+4.3f | %+4.3f | %+4.3f | %+4.3f  \n',[i V(i)  theta(i)*180/pi  Pg(i) Qg(i)-Q(i)*Sb Pd(i)  Qd(i)])
   end
% %---------------------------------------------------------------
% % Cálculo do subsistema 3)
% %---------------------------------------------------------------
% %            calc_fluxo(G,B,V,theta,Y_Sh,N_i,N_f,Sb);
    
