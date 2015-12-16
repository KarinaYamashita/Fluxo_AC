% Cálculo do fluxo de potència pelo método de Newton-Raphson.
close all
clc

% Lista dos arquivos com os dados dos sistemas. Comentar a linha para remover o
% sistema do cálculo.
files = [];
files = [files; 'dados/IEEE_014.txt'];
files = [files; 'dados/IEEE_030.txt'];
files = [files; 'dados/IEEE_057.txt'];
files = [files; 'dados/IEEE_118.txt'];
% files = [files; 'dados/IEEE_300.txt'];

% Tolerância para o erro.
tol = 0.003;

% Mostrar resultados na tela?
mostrar_resultados = false;

% Guarda as informações da simulação na seguinte ordem:
% Num_Barras, Num_linhas, Iterações, Tempo;
informacoes = zeros(length(files(:,1)), 4);

% Calcula o fluxo para os vários sistemas.
for sistema=1:length(files(:,1))
  % Número de barras e linhas do sistema.
  Num_Barras = textread(files(sistema,:),'%d',1,'whitespace', ...
    'numero_de_barras :    ','headerlines',0);
  Num_linhas = textread(files(sistema,:),'%d',1,'whitespace', ...
    'numero_de_lineas :    ','headerlines',1);
  
  % Ligações e impedência entre linhas.
  [N_i,N_f,R,X,Y_shunt] = textread(files(sistema,:)', ...
    '%n%n%n%n%n','headerlines',Num_Barras+7);
  z=[N_i,N_f,R,X,Y_shunt];
  
  % Classificacão de barras por tipo.
  [~,Tipo,V,Pg,Qg,Pd,Qd,~] = textread(files(sistema,:),'%n%n%n%n%n%n%n%n', ...
    Num_Barras,'headerlines',6);
  PQ = find(Tipo==0);
  PV = find(Tipo==2);
  Slack = find(Tipo==3);
  
  % Número de cada tipo de barra.
  NPQ= length(PQ);
  NPV= length(PV);
  
  % Inicialização do valor de módulo e fase das tensões nas barras.
  % As barras PQ são as únicas barras onde irá variar o valor de tensão.
  V(PQ') = 1;
  theta = zeros(Num_Barras,1);
  
  % Potências nas barras
  Sb = 100;%para sistemas que não estão em pu
  %   Sb = 1;%para sistemas que estão em pu
  Pe=zeros(Num_Barras,1);
  Qe=zeros(Num_Barras,1);
  for i=1:Num_Barras
    Pe(i) = (Pg(i)-Pd(i))/Sb;
    Qe(i) = (Qg(i)-Qd(i))/Sb;
  end
  
  % Inicialização de variaveis
  it   = 0;
  Pind = sort([PQ;PV])';
  Qind = sort(PQ)';
  
  % Cálculo do subsistema 1(Calcula V e theta)
  % Calcula a admitância do sistema
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
  
  Y_Bus=-Y_Bus;
  G = real(Y_Bus);
  B = imag(Y_Bus);
  P = calc_P(Num_Barras,G,B,V,theta);
  Q = calc_Q(Num_Barras,G,B,V,theta);
  J = calc_J(Num_Barras,G,B,V,theta,Slack,PV,Q,P);
  
  if NPQ~=0
    erro_m=max(abs([Pe-P;Qe-Q]));
  else
    erro_m = max(abs(Pe-P));
  end
  
  tic
  while erro_m>tol
    J=calc_J(Num_Barras,G,B,V,theta,Slack,PV,Q,P);
    Delta=[Pe-P;Qe-Q];
    b=size(Delta);
    D  = J*Delta;
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
  
  % Salva as informações da simulação.
  % toc vem primeiro para não adicionar tempo extra ao relógio.
  informacoes(sistema, 4) = toc;
  informacoes(sistema, 1) = Num_Barras;
  informacoes(sistema, 2) = Num_linhas;
  informacoes(sistema, 3) = it;
  
  % Gera resultados do subsistema 1.
  if mostrar_resultados
    disp('Num_Barras|  V     |    T   |     Pg        Qg        Pd        Qd')
    for i=1:Num_Barras
      fprintf('   %2d     | %+4.3f | %+4.3f | %+4.3f | %+4.3f | %+4.3f | %+4.3f  \n',[i V(i)  theta(i)*180/pi  Pg(i) Qg(i)-Q(i)*Sb Pd(i)  Qd(i)])
    end
  end
  
  barFile = strcat('result/', num2str(Num_Barras), 'barra.tex');
  barLabel = {'Barra', 'V [V]', 'Fase [graus]', 'Pg [W]', 'Qg [VAR]', 'Pd [W]', 'Qd [VAR]'};
  barras = 1:1:Num_Barras;
  barResult = [barras' V  theta*180/pi  Pg Qg-Q*Sb Pd  Qd]
  matrix2latex(barResult, char(barFile), 'columnLabels', barLabel);
  
  % % Cálculo do subsistema 2(Calcula V e theta)
  %     Pind = Slack';
  %     Qind = [PV;Slack]';
  
  
  
  % Cálculo do subsistema 3)
  calc_fluxo(G,B,V,theta,Y_Sh,N_i,N_f,Sb);
end

infoLabel = {'Num. Barras', 'Num. Linhas', 'Iterações', 'Tempo [s]'};
matrix2latex(informacoes, 'result/info.tex', 'columnLabels', infoLabel);


