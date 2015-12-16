function calc_fluxo(G,B,V,theta,Y_Sh,N_i,N_f,Sb)
    fprintf('De  Para     Pkm       Pmk       Qkm       Qmk           Perdas\n')
    for i = 1:length(N_i)
        k = N_i(i);
        m = N_f(i);
        Pkm(k) = V(k)*V(k) * G(k,m) - V(k)*V(m)*  ( G(k,m)*cos(theta(k)-theta(m)) + B(k,m)*sin(theta(k)-theta(m)) );
        Qkm(k) = -V(k)^2 * (  B(k,m) + Y_Sh(k,m) ) - V(k)*V(m) * ( G(k,m)*sin(theta(k)-theta(m)) - B(k,m)*cos(theta(k)-theta(m)) );

        Pmk(k) = +V(m)^2 * G(k,m) - V(k)*V(m) * ( G(k,m)*cos(theta(k)-theta(m)) - B(k,m)*sin(theta(k)-theta(m)) );
        Qmk(k) = -V(m)^2 * (  B(k,m) + Y_Sh(k,m) ) + V(k)*V(m) * ( G(k,m)*sin(theta(k)-theta(m)) + B(k,m)*cos(theta(k)-theta(m)) );
        %Pt = Pt + Pkm + Pmk;
        %Qt = Qt + Qkm + Qmk;
        %fprintf('%2d   %2d    %+4.3f    %+4.3f    %+4.3f    %+4.3f     %+4.3f    %+4.3f \n',a(i),B(i),Pkm,Qkm,Pmk,Qmk,Pkm+Pmk,Qkm+Qmk)
        fprintf('%2d & %2d & %+4.3f & %+4.3f & %+4.3f & %+4.3f & %+4.3f & %+4.3f  \n',N_i(i),N_f(i),Pkm(k)*Sb,Pmk(k)*Sb,Qkm(k)*Sb,Qmk(k)*Sb,(Pkm(k)+Pmk(k))*Sb,(Qkm(k)+Qmk(k))*Sb)
    end
end