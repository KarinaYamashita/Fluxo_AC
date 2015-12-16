file = 'IEEE_300.txt';
% NB = textread(file,'%d',1,'whitespace','numero_de_barras :    ','headerlines',0);
% [a,Tipo,V,Pd,Qd,Pg,Qg,h] = textread(file,'%n%n%n%n%n%n%n%n',NB,'headerlines',6);
% [de,para,c,d,e] = textread(file','%n%n%n%n%n','headerlines',NB+7);
% z=[de,para,c,d,e];
% size(de)
% size([1:length(de)]')
% [de [1:length(de)]'] 


clc

a = [1 2 3 4 5 16 17 18 19]

b = [5 4 2 3 1 19 16 18 17]



[~,ind] = sort(b)
[~,ind] = sort(ind)
c = 1:length(b);









%     net = zeros(NB,NB);
%     for i = 1:length(a)
%     	net(a(i),b(i)) = 1;
%     end
%     net = net + net' + eye(NB,NB);
% %     G(abs(G)>0) = 1;
% %     B(abs(B)>0) = 1;
% %     sum(sum(net~=G))
% %     sum(sum(net~=B))