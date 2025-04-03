function PNL_studio(f,V,A_b,P)
%
%se P=[] trova min e max assoluto su dominio
%se P contiene dei punti (x,y) n*2 studia quei punti controllandone
% moltiplicatori e cosa sono
    
%richiede funzione sistema_lkkt.m

    conv=false;
    conc=false;

    A=[];
    b=[];
    vertici=V;

    if ~isempty(V) && isempty(A_b)

        %fprintf("\nV={ ")
        %for i=1:size(vertici,1)
        %    fprintf("(%d,%d) ",V(i,1),V(i,2));
        %end
        %fprintf("}\n")

        num_vert=size(V,1);
        fprintf("\n VINCOLI:\n")
        for i=1:num_vert
            [A(i,1),A(i,2),b(i)]= diseq([V(i,1:2),V(mod(i,num_vert)+1, 1:2), V(i,3)]);
        end
        
        b=b';
        A_b=[A b];


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
        %fprintf("\nV={ ")
        %for i=1:size(vertici,1)
        %    fprintf("(%d,%d) ",V(i,1),V(i,2));
        %end
        %fprintf("}\n")
    else
        fprintf("A_b o V errato (uno deve essere vuoto e l'altro pieno)\n");
        return;
    end

    num_vert=size(V,1);
    num_moltiplicatori=size(A,1);

    grad_f=[2*f(1),f(2),f(4);f(2),2*f(3),f(5)];
    Hf=2*[f(1) f(2)/2; f(2)/2 f(3)];
    display(sym(Hf),"Hf")
    auto_val=eig(Hf);
    [num,den]=rat(auto_val);
    fprintf("\nautovalori:\n");
    for i=1:size(num,1)
        fprintf("%+d/%+d\n",num(i),den(i));
    end 
    %display(sym(auto_val),"autovalori");

    if all(auto_val >= 0)
        conv=true;
        fprintf("FUNZIONE CONVESSA -> MAX su un VERTICE\n");
    elseif all(auto_val <= 0)
        conc=true;
        fprintf("FUNZIONE CONCAVA -> MIN su un VERTICE\n");
    else
        fprintf("FUNZIONE NON CONCAVA NON CONVESSA\n");
    end

    fprintf("-------------------------------------------------------------------\n\n");
    %matrice_moltiplicatori=[]; % ogni riga contiene i moltiplicatori di uno dei punti studiati
    for i=1:size(P,1)
        fprintf("PUNTO P=[%s,%s]:\t",sym(P(i,1)),sym(P(i,2)) )
        p=(P(i,:))';
        interno = true;
        esterno = false;
        %verifico se il punto è interno
        for j=1:size(A,1)
            if A(j,:)*p > b(j)
                interno = false;
                esterno = true;
            elseif A(j,:)*p == b(j)
                interno = false;
            end
        end
        
        grad_f_p = grad_f(:,1:2)*p+grad_f(:,3);
        if interno
            fprintf("PUNTO INTERNO\n")
            if ~all(grad_f_p==0)
                display(grad_f_p, 'grad_f(p)');
                fprintf("INTERNO + NON ANNULLA IL GRADIENTE ->NON STAZIONARIO\n");
                continue;
            end
        elseif esterno
            fprintf("PUNTO ESTERNO\n");
            continue;
        elseif ismember(p',V(:,1:2),'rows')
            fprintf("VERTICE\n")
        else
            fprintf("PUNTO DI FRONTIERA, NON vertice\n")
        end
        moltiplicatori=sistema_lkkt(grad_f,A_b,p');
        
        fprintf("Classificazione: ")
        if all(moltiplicatori==0)
            fprintf("SELLA, MIN o MAX\n")
        elseif all(moltiplicatori>=0)
            fprintf("MIN (l o g) o SELLA\n");
        elseif all(moltiplicatori<=0)
            fprintf("MAX (l o g) o SELLA\n");
        else
            fprintf("SELLA\n");
        end
    

    end
    
    %studio max e min globali sul poliedro
    fprintf("-------------------------------------------------------------------\n\n");
    arr_val_vert=[];
    if conc
        fprintf("Funzione CONCAVA -> uno dei vertici è MIN GLOBALE\n")
        for j = 1: size(V,1)
            val= f*[(V(j,1))^2; V(j,1)*V(j,2); (V(j,2))^2; V(j,1); V(j,2); 1]; %valori di f nei vertici
            arr_val_vert=[arr_val_vert;val];    
        end
        [~,ind_min]=min(arr_val_vert);
        for j=1:size(V,1)
            fprintf("f(%d,%d)= %+d\t",V(j,1),V(j,2),arr_val_vert(j));
            if j==ind_min
                fprintf(" <- MIN")
            end
            fprintf("\n");
        end
        p=V(ind_min,1:2)';
        fprintf("\nX_min=(%d,%d)\n",p(1),p(2));
        
        moltiplicatori=sistema_lkkt(grad_f,A_b,p');
        
        %cerco max
        %valuto gradiente = 0 se non esiste o è fuori
        %procedo con le restrizioni
        fprintf("cerco MAX ->\n grad_f=0:\n");
        fprintf("\t %+d x1  %+d x2 = %+d \n", grad_f(1,1), grad_f(1,2), -grad_f(1,3));
        fprintf("\t %+d x1  %+d x2 = %+d \n\n", grad_f(2,1), grad_f(2,2), -grad_f(2,3));
        
        %gestire il caso no sol e infinite sol (dovrebbe funzionare)
        if rank(grad_f(:,1:2)) < 2 %(rango non massimo perchè siamo in r^2)->non ci sono sol di gr(f)==0
            fprintf("det=0 -> procedo con le restrizioni\n"); %trovo la soluzione sia che ci siano infinite sol che 0
            p=calcola_m_su_restrizioni(1,f,V);
        else
            X_max= grad_d(:,1:2)\(-grad_f(:,3));
            if all(A*X_max <= b)
                fprintf("\t X_max = (%s, %s)",sym(X_max(1)), sym(X_max(2)));
            else
                fprintf("Punto che annulla il gradiente esterno -> procedo con restrizioni\n")
                p=calcola_m_su_restrizioni(1,f,V);
            end
        end 
        p=p';
        fprintf("\nX_max=(%s,%s)\n",sym(p(1)),sym(p(2)));
        sistema_lkkt(grad_f,A_b,p');
        %per max trovo ul punto che annulla il gradiente se è appartine al
        %poliedro è il min altirmenti trovo il min usando le restrizioni

        %trovo min su R con quadprog
        %trovo max su r con quadprog(non dovrebbe esistere e lo dimostro
        %con restrizione)
        
    elseif conv
        fprintf("Funzione CONVESSA -> uno dei vertici è MAX GLOBALE\n")
        for j = 1: size(V,1)
            val= f*[(V(j,1))^2; V(j,1)*V(j,2); (V(j,2))^2; V(j,1); V(j,2); 1]; %valori di f nei vertici
            arr_val_vert=[arr_val_vert;val];        
        end
        [~,ind_max]=max(arr_val_vert);
        for j=1:size(V,1)
            fprintf("f(%d,%d)= %+d\t",V(j,1),V(j,2),arr_val_vert(j));
            if j==ind_max
                fprintf(" <- MAX")
            end
            fprintf("\n");
        end
        p=V(ind_max,1:2)';
        fprintf("\nX_max=(%d,%d)\n",p(1),p(2));
        moltiplicatori=sistema_lkkt(grad_f,A_b,p');

        %cerco min
        %valuto gradiente = 0 se non esiste o è fuori
        %procedo con le restrizioni
        fprintf("cerco MIN ->\n grad_f=0:\n");
        fprintf("\t %+d x1  %+d x2 = %+d \n", grad_f(1,1), grad_f(1,2), -grad_f(1,3));
        fprintf("\t %+d x1  %+d x2 = %+d \n\n", grad_f(2,1), grad_f(2,2), -grad_f(2,3));
        
        %gestire il caso no sol e infinite sol (dovrebbe funzionare)
        if rank(grad_f(:,1:2)) < 2 %(rango non massimo perchè siamo in r^2)->non ci sono sol di gr(f)==0
            fprintf("det=0 -> procedo con le restrizioni\n"); %trovo la soluzione sia che ci siano infinite sol che 0
            p=calcola_m_su_restrizioni(0,f,V);
        else
            X_min= grad_f(:,1:2)\(-grad_f(:,3));
            if all(A*X_min <= b)
                fprintf("\t X_min = (%s, %s)",sym(X_min(1)), sym(X_min(2)));
            else
                fprintf("Punto che annulla il gradiente esterno -> procedo con restrizioni\n")
                p=calcola_m_su_restrizioni(0,f,V);
            end
        end 
        p=p';
        fprintf("\nX_min=(%s,%s)\n",sym(p(1)),sym(p(2)));
        sistema_lkkt(grad_f,A_b,p');
    
    else %funzione nè concava nè convessa
        %uso gradiente e restrizioni 
        %cerco max
        %valuto gradiente = 0 se esiste lo aggiungo ai punti da controllare
        %altrimenti procedo con sole restirzioni
        punti_max=[]; %punto,valore
        fprintf("cerco MAX \n grad_f=0:\n"); 
        fprintf("\t %+d x1  %+d x2 = %+d \n", grad_f(1,1), grad_f(1,2), -grad_f(1,3));
        fprintf("\t %+d x1  %+d x2 = %+d \n\n", grad_f(2,1), grad_f(2,2), -grad_f(2,3));
        
        %gestire il caso no sol e infinite sol (dovrebbe funzionare)
        
        if rank(grad_f(:,1:2)) < 2 %(rango non massimo perchè siamo in r^2)->non ci sono sol di gr(f)==0
            fprintf("det=0 -> procedo con le restrizioni\n"); %trovo la soluzione sia che ci siano infinite sol che 0
            p=calcola_m_su_restrizioni(1,f,V);
            punti_max=[ p , f(1)*p(1)*p(1)+f(2)*p(1)*p(2)+f(3)*p(2)*p(2)+f(4)*(p1)+f(5)*p(2)+f(6) ];
        else
            X_grad= grad_f(:,1:2)\(-grad_f(:,3)); %x che annulla il gradiente
            if all(A*X_grad <= b)
                p=X_grad';
                val=f(1)*p(1)*p(1)+f(2)*p(1)*p(2)+f(3)*p(2)*p(2)+f(4)*p(1)+f(5)*p(2)+f(6) ;
                fprintf("\tX = (%s, %s) -> f(X)= %s\n",sym(X_grad(1)), sym(X_grad(2)), sym(val) );
                punti_max=[ p , val];
                p=calcola_m_su_restrizioni(1,f,V);
                punti_max=[ punti_max; p, f(1)*p(1)*p(1)+f(2)*p(1)*p(2)+f(3)*p(2)*p(2)+f(4)*p(1)+f(5)*p(2)+f(6) ];
            else
                fprintf("Punto che annulla il gradiente esterno -> procedo con restrizioni\n")
                p=calcola_m_su_restrizioni(1,f,V);
                punti_max=[ punti_max; p , f(1)*p(1)*p(1)+f(2)*p(1)*p(2)+f(3)*p(2)*p(2)+f(4)*p(1)+f(5)*p(2)+f(6) ];
            end
        end 
        [val_max, ind_max] = max(punti_max(:,3));
        p=punti_max(ind_max,1:2);
        fprintf("\nX_max=(%s,%s)\n",sym(p(1)),sym(p(2)));
        sistema_lkkt(grad_f,A_b,p);
        
        %%MIN%%
        punti_min=[]; %punto,valore
        fprintf("cerco MIN \n grad_f=0:\n"); 
        fprintf("\t %+d x1  %+d x2 = %+d \n", grad_f(1,1), grad_f(1,2), -grad_f(1,3));
        fprintf("\t %+d x1  %+d x2 = %+d \n\n", grad_f(2,1), grad_f(2,2), -grad_f(2,3));
        
        %gestire il caso no sol e infinite sol (dovrebbe funzionare)
        
        if rank(grad_f(:,1:2)) < 2 %(rango non massimo perchè siamo in r^2)->non ci sono sol di gr(f)==0
            fprintf("det=0 -> procedo con le restrizioni\n"); %trovo la soluzione sia che ci siano infinite sol che 0
            p=calcola_m_su_restrizioni(0,f,V);
            punti_min=[ p , f(1)*p(1)*p(1)+f(2)*p(1)*p(2)+f(3)*p(2)*p(2)+f(4)*(p1)+f(5)*p(2)+f(6) ];
        else
            X_grad= grad_f(:,1:2)\(-grad_f(:,3)); %x che annulla il gradiente
            if all(A*X_grad <= b)
                p=X_grad';
                val=f(1)*p(1)*p(1)+f(2)*p(1)*p(2)+f(3)*p(2)*p(2)+f(4)*p(1)+f(5)*p(2)+f(6) ;
                fprintf("\t X = (%s, %s) -> f(X)= %s\n",sym(X_grad(1)), sym(X_grad(2)), sym(val) );
                punti_min=[ p , val];
                p=calcola_m_su_restrizioni(0,f,V);
                punti_min=[ punti_min; p, f(1)*p(1)*p(1)+f(2)*p(1)*p(2)+f(3)*p(2)*p(2)+f(4)*p(1)+f(5)*p(2)+f(6) ];
            else
                fprintf("Punto che annulla il gradiente esterno -> procedo con restrizioni\n")
                p=calcola_m_su_restrizioni(1,f,V);
                punti_min=[ punti_min; p , f(1)*p(1)*p(1)+f(2)*p(1)*p(2)+f(3)*p(2)*p(2)+f(4)*p(1)+f(5)*p(2)+f(6) ];
            end
        end 
        [val_min, ind_min] = min(punti_min(:,3));
        p=punti_min(ind_min,1:2);
        fprintf("\nX_min=(%s,%s)\n",sym(p(1)),sym(p(2)));
        sistema_lkkt(grad_f,A_b,p);

    end 

end