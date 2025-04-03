function [ris] = calcola_m_su_restrizioni(M,f,V)
    lavoro=[V(:,1:2);V(1,1:2)];
    valori=[];
    punti=[];
    for i=1:size(lavoro,1)-1
        dk=(lavoro(i+1,:)-lavoro(i,:))';
        %x1= r1+r2t
        r1= lavoro(i,1);
        r2= dk(1);
        %x2= r3+r4t
        r3= lavoro(i,2);
        r4= dk(2);

        phi_t2 = f(1)*r2*r2 + f(2)*r2*r4 + f(3)*r4*r4;
        phi_t = 2*f(1)*r1*r2 + f(2)*r1*r4 + f(2)*r2*r3 + 2*f(3)*r3*r4 + f(4)*r2 + f(5)*r4;
        phi_1= f(1)*r1*r1 + f(2)*r1*r3 + f(3)*r3*r3 + f(4)*r1 + f(5)*r3 + f(6);
        
        fprintf('\tφ(t)=f(%s%s%s%st, %s%s%s%st)=\t',c_segno(sign(sym(r1))),abs(sym(r1)),c_segno(sign(sym(r2))),abs(sym(r2)),c_segno(sign(sym(r3))),abs(sym(r3)),c_segno(sign(sym(r4))),abs(sym(r4)));
        fprintf('%s t^2 %s%s t %s%s\t\t',sym(phi_t2), c_segno(sign(sym(phi_t))), abs(sym(phi_t)), c_segno(sign(sym(phi_1))), abs(sym(phi_1)) );
        
        vert=-phi_t/(2*phi_t2);
        vect_f=[0,phi_1; 1, phi_t2+phi_t+phi_1];
        if vert>0 && vert<1
            vect_f =[vect_f; [vert, (phi_t2*vert^2+phi_t*vert+phi_1)] ];
        end
        f_m=NaN;
        if M==1
            [val, indice_massimo] = max(vect_f(:, 2));
            tk=vect_f(indice_massimo,1);
            f_m=val;
        else
            [val, indice_min] = min(vect_f(:, 2));
            tk=vect_f(indice_min,1);
            f_m=val;
        end
        valori=[valori,f_m];
        Xm=lavoro(i,:)'+tk*dk;
        punti=[punti;Xm'];
        fprintf("-> X=(%s,%s), f(x)= %s \n",sym(Xm(1)),sym(Xm(2)),sym(f_m));
    end

    if M==1
        [val,ind]=max(valori);
        X_max=punti(ind,:);
        ris=X_max;
    else
        [val,ind]=min(valori);
        X_min=punti(ind,:);
        ris=X_min;
    end

end
