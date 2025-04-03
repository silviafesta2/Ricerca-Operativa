%funzione che dato il problema in formato primale standard
%stampa tutti i piani di taglio di Gomory
% c è vettore colonna

% NOTA IMPORTANTE
% puoi inserire la matrice di un primale senza i vincoli di positivita
% delle variabili PURCHE LE VARIABILI ABBIANO IL VINCOLO DI POSITIVITA


function Gomory(c,A,b)
% porto il problema in formato duale standard
fprintf("Passo 1: porto il problema in formato duale standard\n");
 [m,n] = size(A);
 Aeq = [A,eye(m)];
 f = [-c;zeros(m,1)];

 display(sym(Aeq),"Aeq");
 display(sym(f'),"c");
 
 fprintf("{ min c * x\n{ Ax = b\n{ x >= 0\n");
 [x,v] = linprog(f,[],[],Aeq,b,zeros(m + n,1),[]);
 display(sym(x),'X');

 % calcolo x frazionario
 x_frac = x - floor(x);
 indiciFrazionari = find(x_frac > 0);
 
 % se X e' a componenti intere allora sono all'ottimo anche per beni
 % indivisibili
 if all(x_frac == 0)
    fprintf("\tx e' a componenti intere quindi sono all'ottimo per beni indivisibili\n");
    return;
 end

 % se invece esiste una componente frazionaria allora posso fare dei piani
 % di taglio
 fprintf("Passo 2: calcolo i piani di taglio per le componenti frazionarie di x\n");
 indiciPositivi = find(x > 0);
 N = find(x == 0);

 A_tilde = calcolaGomory(Aeq,indiciPositivi);
 A_tilde = sym(A_tilde);

 % piani prima delle semplificazioni
 piani = [A_tilde,x_frac(indiciPositivi)];

 % piani dopo le semplificazioni
 piani_lavoro = piani;
 [~,denominatori] = rat(piani_lavoro);

 % devo prima trovarmi quanto valgono le variabili che ho aggiunto al
 % problema
 A_lavoro = [-A,b];

 % dove i denominatori sono nulli ci metto NaN
 % in modo che se ho degli interi non prendo il denominatore dell'intero
 denominatori(denominatori == 0) = NaN;
 
 % coefficienti per cui moltiplicare ogni riga
 coeff_moltiplicare = zeros(1,size(denominatori,1));

 % prima devo stabilire quali sono i coefficienti poi posso moltiplicare
 for i = 1:size(coeff_moltiplicare,2)
    coeff_moltiplicare(i) = lcm(sym(denominatori(i,:)));
    piani_lavoro(i,:) = piani_lavoro(i,:) * coeff_moltiplicare(i);
 end

 % adesso in piani lavoro abbiamo i piani semplificati
 % stampo comunque la rappresentazione simbolica
 piani_lavoro_simbolic = sym(piani_lavoro);


 % ottengo una matrice in cui le righe corrispondono ali indici frazionari
 % le colonne corrispondono alle variabili non di base
 % l'ultima colonna sono i termini noti

 fprintf("Piani di taglio di Gomory:\n\n")
    for i = 1:size(indiciPositivi,1)
        fprintf("\t%d: ",indiciPositivi(i));
        for j = 1:size(A_tilde,2)
            if j ~=1
                if (A_tilde(i,j) > 0)
                    fprintf(" +%s x%d ",sym(A_tilde(i,j)),N(j));
                else
                    fprintf(" %s x%d ",sym(A_tilde(i,j)),N(j));
                end
            else
                fprintf(" %s x%d ",sym(A_tilde(i,j)),N(j));
            end
        end
        fprintf (" >= %s\n",sym(x_frac(indiciPositivi(i))));
    end

    % stampo i piani semplificati
    fprintf("\nPiani di taglio di Gomory semplificati:\n\n")
    for i = 1:size(indiciPositivi,1)
        fprintf("\t%d: ",indiciPositivi(i));
        for j = 1:size(A_tilde,2)
            if j ~=1
                if (A_tilde(i,j) > 0)
                    fprintf(" +%s x%d ",piani_lavoro_simbolic(i,j),N(j));
                else
                    fprintf(" %s x%d ",piani_lavoro_simbolic(i,j),N(j));
                end
            else
                fprintf(" %s x%d ",piani_lavoro_simbolic(i,j),N(j));
            end
        end
        fprintf (" >= %s\n",piani_lavoro_simbolic(i,end));
    end

    % adesso devo essere in grado di passare questi piani al primale ho
    % aggiunto m variabili al problema, quindi se nei piani di taglio trovo
    % delle variabili del genere li devo portare al primale 
    
    fprintf("\nAdesso portiamo i piani in formato primale:\n\n");
    
    % devo prima scrivermi la definizione di ogni variabile slack

    for i = 1:m
        fprintf("\tx%d = ",i + n)
        for j = 1:(size(A_lavoro,2) - 1)
            if (j == 1)
                fprintf("%s x%d ",sym(A_lavoro(i,j)),j);
            else
                if (A_lavoro(i,j) < 0)
                    fprintf("%s x%d ",sym(A_lavoro(i,j)),j);
                else
                    fprintf("+%s x%d" ,sym(A_lavoro(i,j)),j);
                end
            end
        end
        fprintf("+ %s\n",sym(A_lavoro(i,end)))
    end

    
    % matrice di appoggio che contiene i coefficienti per ogni variabile
    % nei piani di gomory
    componenti_sommate = [zeros(size(piani,1),size(x,1) + 1)];

    for i = 1:size(piani_lavoro_simbolic,1)
    
        for j = 1:size(piani_lavoro_simbolic,2)

            % controllo se siamo a b
            if j == size(piani_lavoro_simbolic,2)
                componenti_sommate(i,end) = componenti_sommate(i,end) + piani_lavoro_simbolic(i,end);
                continue;
            end

            % per ogni componente di A_tilde, se questa e' una slack
            % variabile devo sommare in componenti_sommate la sua
            % definizione per il suo coefficiente

            % la mappa per ottenere da j all'indice della variabile e' N
            % sostanzialmente N(j) -> indice variabile
            if (N(j) > n)
                % e' una slack variable

                for ii = 1:n
                    componenti_sommate(i,ii) = componenti_sommate(i,ii) + piani_lavoro_simbolic(i,j) * A_lavoro(N(j)-n,ii);
                end

                % devo sommare b

                componenti_sommate(i,end) = componenti_sommate(i,end) - A_lavoro(N(j)-n,end) * piani_lavoro_simbolic(i,j);

            else
                % non e' una slack variable (non agisce neanche su b)
                componenti_sommate(i,j) = componenti_sommate(i,j) + piani_lavoro_simbolic(i,j);

            end

        end
    
    end

    divisori = zeros(size(componenti_sommate,1),1);

    for i = 1:size(componenti_sommate,1)
        divisori(i) = gcd(sym(componenti_sommate(i,:)));
    end


    % a questo punto abbiamo le equazioni dei piani
    % dentro la variabile
    componenti_sommate;

    % stampiamo le equazioni dei piani non semplificati
    fprintf("\nEquazioni dei piani al primale non semplificati\n\n")
    for i = 1:size(componenti_sommate,1)
        fprintf("\t%d: ",indiciPositivi(i));
        for j = 1:size(componenti_sommate,2) - 1

            if j == 1
            
                fprintf("%s x%d ",sym(componenti_sommate(i,j)),j)
                continue;

            end

            % se non c'e' questa componente non la stampo
            if componenti_sommate(i,j) == 0
                continue;
            end

            if componenti_sommate(i,j) < 0
                fprintf(" %s x%d ",sym(componenti_sommate(i,j)),j)
            else
                fprintf(" +%s x%d",sym(componenti_sommate(i,j)),j)
            end

        end

        % adesso stampiamo b
        fprintf(" >= %s\n",sym(componenti_sommate(i,end)))
    end

    % adesso pensiamo a semplificarli
    divisori;

    componenti_sommate_lavoro = componenti_sommate;

    for i = 1:size(componenti_sommate_lavoro,1)
        componenti_sommate_lavoro(i,:) = componenti_sommate_lavoro(i,:) / divisori(i) * -1;
    end


    % stampiamo le equazioni dei piani semplificati
    fprintf("\nEquazioni dei piani al primale semplificati\n\n")
    for i = 1:size(componenti_sommate_lavoro,1)

        % modifica barbara
        if isnan(componenti_sommate_lavoro(i,1))
            continue
        end

        fprintf("\t%d: ",indiciPositivi(i));
        for j = 1:size(componenti_sommate_lavoro,2) - 1


            % se non c'e' questa componente non la stampo
            if componenti_sommate_lavoro(i,j) == 0
                continue;
            end

            if componenti_sommate_lavoro(i,j) < 0
                fprintf(" %s x%d ",sym(componenti_sommate_lavoro(i,j)),j)
            else
                fprintf(" +%s x%d",sym(componenti_sommate_lavoro(i,j)),j)
            end

        end

        % adesso stampiamo b
        fprintf(" <= %s\n",sym(componenti_sommate_lavoro(i,end)))
    end

    fprintf("sol ottima con intlinprog\n");
    [x,v]=intlinprog(-c',1:(size(A,2)),A,b,[],[],zeros(size(A,2),1),[]);
    x
    % ritorniamo la riga da aggiungere al primale

end