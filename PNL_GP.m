function  PNL_GP(f,V,Xk,M,A_b)
%PNL Summary of this function goes here
%   PARAMETRI:
% f -> coeff di x1^2, X1X2, X2^2, X1, X2, 1
% V -> vettore num_vert*2 DEVE ESSERE ORDINATO (in senso antiorario) SENNO'
% NON FUNZIONA (vicino a ogni verice lo zero si riferisce al vincolo che
% parte da quel vertice (senzo antiorario) 0->lo zero è contenuto nel
% semipiano
% Xk ->punto di partenza
% M _> 1-MAX 0-min

    A=[];
    b=[];
    Max=M;
    if ~isempty(V) && isempty(A_b)

        fprintf("\nV={ ")
        for i=1:size(V,1)
            fprintf("(%d,%d) ",V(i,1),V(i,2));
        end
        fprintf("}\n")

        num_vert=size(V,1);
        fprintf("\n VINCOLI:\n")
        for i=1:num_vert
            [A(i,1),A(i,2),b(i)]= diseq([V(i,1:2),V(mod(i,num_vert)+1, 1:2), V(i,3)]);
        end
        b=b';
        A_b=[A,b];
        
    elseif isempty(V) && ~isempty(A_b)
        vertici=[];
        A=A_b(:,1:2);
        b=A_b(:,3);

        for i=1:size(A_b,1)
            for j=i+1:size(A_b,1)

                a1=A_b(i,1);
                a2=A_b(i,2);
                b1=A_b(i,3);

                a3=A_b(j,1);
                a4=A_b(j,2);
                b2=A_b(j,3);
                
                x1=[];
                x2=[];
                if (a1*a4 - a2*a3) == 0  %non è una base
                    continue
                else 
                    x=[a1 a2;a3 a4]\[b1;b2];
                    x1=x(1);
                    x2=x(2);

                %qualcosa non funziona in questa parte commentata
                %elseif a1 %posso dividere per a1
                %    x2=(a1*b2 - a3*b1)/(a1*a4 - a2*a3);
                %    x1=(b1-a2*x2)/a1;
                %elseif a2 %posso dividere per a2
                %    x1=(a4*b1 - a2*b2)/(a1*a4 - a2*a3);
                %    x2=(b1-a2*x1)/a2;
                end

                if all(A*[x1;x2] <= b) %se è ammissibile è un vertice
                    vertici=[vertici; x1 x2 NaN];
                end
            end
        end
        V=ordina_antiorario(vertici, mean(vertici)); %per ordinarli devo farlo rispetto a un punto interno uso il baricentro così e interno sicuro
        fprintf("\nV={ ")
        for i=1:size(vertici,1)
            fprintf("(%d,%d) ",V(i,1),V(i,2));
        end
        fprintf("}\n")
    else
        fprintf("A_b o V errato (uno deve essere vuoto e l'altro pieno)\n");
        return;
    end
    num_vert=size(V,1);
    display(sym(Xk),'Xk');
    Grad_f=[2*f(1),f(2),f(4);f(2),2*f(3),f(5)];
    grad_f_xk = Grad_f(:,1:2)*Xk + Grad_f(:,3);
    fprintf("\n%sf(x)=\t[%+d X1  %+d X2  +%d ]\n \t\t[%+d X1  %+d X2  %+d ]\n",char(8711),Grad_f');
    fprintf("\n%sf(xk)=\t[ %s ]\n \t\t[ %s  ]\n",char(8711),sym(grad_f_xk'));
    
    M=[];
    vert1=[]; %il primo (senso antiorario) dei due vertici fra cui sta xk
    for i=1:size(A,1)
        if A(i,:)*Xk==b(i)
            M=[M;A(i,:)];
            %trovo primo vertice del segmento dove sta Xk
            for vi=[flip(1:num_vert)] %prima assegna il secondo poi viene sovrascritto dal primo (l'ultimo serve a gestire il caso in cui il segmento sta fra 1 e num_vert (num_vert viene prima di 1 in senso antiorario)
                if M*(V(vi,1:2)')==b(i)
                    vert1=vi;
                end
            end
            if vert1==4 && M*(V(1,1:2)')==b(1)
                vert1=4;
            end
        end
    end
    display(M);
    %display(vert1);
    if ~isempty(M)
        H= eye(2)-(M*M')^(-1)* M'*M ;
    else
        H=eye(2);
    end

    display(sym(H),"H");
    
    if Max==1
        dk=H*grad_f_xk;
    else
        dk=H*(-grad_f_xk);
    end
    display(sym(dk),'dk');
    retta=[Xk,dk];

    %contollo contro quale dei due vertici è quello che limita

    %controllo se uno dei due indici è nullo perchè in caso non posso
    %dividere
    ind=1;
    if retta(1,2)==0
        ind=2;
    end

    tk1 =  (V(vert1,ind) - retta(ind,1))/retta(ind,2) ;
    tk2 =  (V( mod(vert1,size(V,1))+1 ,ind) - retta(ind,1))/retta(ind,2) ;
    
    tk_max=tk1;
    v=V(vert1,1:2);
    if  tk1<0
        tk_max=tk2;
        v=V(mod(vert1,size(V,1))+1,1:2);
    end

    display(sym(v'),'V');
    display(sym(tk_max), "tk_max");
    
    r1=retta(1,1);
    r2=retta(1,2);
    r3=retta(2,1);
    r4=retta(2,2);

    phi_t2 = f(1)*r2*r2 + f(2)*r2*r4 + f(3)*r4*r4;
    phi_t = 2*f(1)*r1*r2 + f(2)*r1*r4 + f(2)*r2*r3 + 2*f(3)*r3*r4 + f(4)*r2 + f(5)*r4;
    phi_1= f(1)*r1*r1 + f(2)*r1*r3 + f(3)*r3*r3 + f(4)*r1 + f(5)*r3 + f(6);
    
    
    fprintf('\n\tφ(t)=f(%s%s%s%st, %s%s%s%st)=\n',c_segno(sign(sym(r1))),abs(sym(r1)),c_segno(sign(sym(r2))),abs(sym(r2)),c_segno(sign(sym(r3))),abs(sym(r3)),c_segno(sign(sym(r4))),abs(sym(r4)));
    fprintf('\t\t%s t^2 %s%s t %s%s\n\n\n',sym(phi_t2), c_segno(sign(sym(phi_t))), abs(sym(phi_t)), c_segno(sign(sym(phi_1))), abs(sym(phi_1)) );
    
    fprintf("forma Parabola:\n\n")
    if phi_t2==0
        if phi_t>0
            fprintf('   /   \n  /    ');
        elseif phi_t<0
            fprintf('   \\   \n    \\   ');
        else
           fprintf('   ____   ');
        end
    end

    if phi_t2 > 0
        fprintf('\\       /\n \\     / \n   ***   ');
    elseif phi_t2 < 0
        fprintf('   ___   \n /     \\ \n/       \\');
    end
    fprintf('\n\n');

    vert=-phi_t/(2*phi_t2);    
    vect_f = [0, phi_1; tk_max, phi_t2*tk_max*tk_max + phi_t*tk_max + phi_1];

    if vert>0 && vert<tk_max
        vect_f =[vect_f; vert, phi_t2*vert^2 + phi_t*vert + phi_1];
    end
    
    if Max==1
        [val, indice_massimo] = max(vect_f(:, 2));
        tk=vect_f(indice_massimo,1);
    else
        [val, indice_min] = min(vect_f(:, 2));
        tk=vect_f(indice_min,1);
    end
    
    display(sym(vert),'-b/2a');
    display(sym(tk),'tk');
    Xnew=Xk+tk*dk;
    display(sym(Xnew),'X_k+1');
    display(sym(val),'f(X_k+1)');

end