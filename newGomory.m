% indicato il problema primale standard
% esclusi i vincoli di positivita delle variabili
function newGomory(c,A,b)

    format longg;

    [m,n] = size(A);

    % c puo' essere passato sia come vettore riga che come vettore colonna
    [cm,cn] = size(c);

    if (cm > cn)
            c = c';
    end

    % ho portato il vettore c a vettore di riga
    Aeq = [A,eye(m)];
    ceq = [c,zeros(1,m)];

    % dovrei stampare le seguenti informazioni:
    % 1. aeq
    % 2. ceq
    % 3. beq
    fprintf("Porto il problema in formato duale standard:\n")
    display(sym(Aeq),"Aeq");
    display(sym(b),"beq");
    display(sym(-ceq),"c");

    [xrc,vrc] = linprog(-ceq,[],[],Aeq,b,[zeros(1,m+n)],[]);
    % se la soluzione ha componenti intere allora devo dire che siamo
    % all'ottimo per componenti intere
    xrc_lavoro = xrc - floor(xrc);
    if (all(xrc_lavoro == 0))
        
        fprintf("Siamo all'ottimo per componenti intere\n");
        return;
    
    end


    % informazioni che dovrei stampare
    % 1. ottimo del rilassato continuo
    % 2. valore ottimo del rilassato continuo
    display(sym(xrc),"Ottimo del rilassato continuo")
    display(-sym(vrc),"Valore ottimo del rilassato continuo")

    % devo trovare la base della soluzione
    indiciNonNulli = find(xrc ~= 0)';

    % tuttavia la soluzione potrebbe essere degenere quindi potrebbe essere
    % necessario scegliere altri indici oltre a quelli diversi da 0
    % in totale ho bisogno di 2m - n indici in base
    indiciNulli = find(xrc == 0)';

    B = indiciNonNulli;

    % cardinalita' di indiciNulli
    cardNonNulli = size(indiciNonNulli,2);

    if cardNonNulli < (m)
        % ottengo il numero di componenti rimanenti per formare una base
        B = [B,indiciNulli(1,1:(2 * m - n - cardNonNulli))];
        B = sort(B);
    end

    N = [1:(m + n)];
    N = setdiff(N,B);
    
    % adesso abbiamo la base e l'insieme degli indici non di base
    B;
    N;

    % ottengo Ab, An, ed A tilde
    Ab = Aeq(:,B);
    An = Aeq(:,N);
    Atilde = inv(Ab) * An;
    
    for i = 1:size(Atilde,1)
        
        for j = 1:size(Atilde,2)
        
            if (abs(Atilde(i,j)) < 10e-6)
                Atilde(i,j) = 0;
            end
        end
    
    end

    display(sym(Atilde),"Matrice Atilde")
    variabili_lavoro = [-A,b];
    % ottengo A tilde frazionario
    Atilde_frac = Atilde - floor(Atilde);

    % sistemo per possibili errori di macchina
    for i = 1:size(Atilde_frac,1)
        for j = 1:size(Atilde_frac,2)
            [num,den] = rat(Atilde_frac(i,j));
            if (num/den < 10e-5)
                Atilde_frac(i,j) = 0;
            end
        end
    end

    display(sym(Atilde_frac),"Matrice Atilde frazionaria")
    % ottengo il vettore dei termini noti
    bb = xrc(B) - floor(xrc(B));




    % stampa i piani di Gomory
    fprintf("\nStampo i piani di Gomory:\n")

    for i = 1:m
    
        fprintf("%d: ",B(i))

        for j = 1:n
        
            if (j ~= 0) && (Atilde_frac(i,j) >= 0)
                % devo stampare anche il segno
                fprintf(" +")
            end
            % devo stampare il valore della variabile
            fprintf("%sx%d",sym(Atilde_frac(i,j)),N(j))
        
        end
        
        [num,den] = rat(bb(i));
        if (den > 1)
            fprintf(" >= %d/%d\n",num,den);
        else
            fprintf(" >= %d\n",num);
        end
    end


    % semplifico i piani di gomory


    % ottengo il minimo comune multiplo dei denominatori
    piani_lavoro = [Atilde_frac,bb];
    [~,denominatori] = rat(piani_lavoro);

    % ottengo le righe semplificate moltiplicando per il denominatore
    % comune
    for i = 1:size(B,2)
        % ottengo mcm della riga
        mcm = lcm(sym(denominatori(i,:)));
        piani_lavoro(i,:) = piani_lavoro(i,:) * mcm;

        % adesso semplifico ulteriormente
        [num,~] = rat(piani_lavoro(i,:));
        mcd = gcd(sym(num));

        piani_lavoro(i,:) = piani_lavoro(i,:) / mcd;
    end 

    % adesso possiedo i piani semplificati
    piani_lavoro;

    % stampo i piani di gomory semplificati
    fprintf("\nStampo i piani di Gomory semplificati:\n")

    for i = 1:m
    
        fprintf("%d: ",B(i))

        for j = 1:n
        
            if (j ~= 0) && (Atilde_frac(i,j) >= 0)
                % devo stampare anche il segno
                fprintf(" +")
            end
            % devo stampare il valore della variabile
            fprintf("%sx%d",sym(piani_lavoro(i,j)),N(j))
        
        end
        
        [num,den] = rat(piani_lavoro(i,n+1));
        if (den > 1)
            fprintf(" >= %d/%d\n",num,den);
        else
            fprintf(" >= %d\n",num);
        end
    
    end


    

    fprintf("\nPorto i piani al primale:\nVariabili ausiliarie:\n")
    
    % stampo le definizioni delle variabili
    for i = 1:m
        fprintf("x%d =",i + n);
        for j = 1:n
            fprintf(" %sx%d",sym(variabili_lavoro(i,j)),j)
        end
        if (variabili_lavoro(i,n+1) >= 0)
            fprintf(" +%s\n",sym(variabili_lavoro(i,n+1)))
        else
            fprintf(" %s\n",sym(variabili_lavoro(i,n+1)))
        end
    end

    % passiamo i piani al formato primale
    accumulatori = zeros(m,n + 1);


    % porto i piani di Gomory al primale
    for i = 1:size(B,2)
        for j = 1:size(N,2)
            % devo vedere se e' una variabile ausiliaria o no
            if N(j) > n
                % e' una variabile ausiliaria
                accumulatori(i,:) = accumulatori(i,:) + piani_lavoro(i,j) * variabili_lavoro(N(j) - n,:);
            else
                % non e' una variabile ausiliaria
                accumulatori(i,N(j)) = accumulatori(i,N(j)) + piani_lavoro(i,j);
            end

        end
    end

    accumulatori(:,n+1) = accumulatori(:,n+1) - piani_lavoro(:,n+1);
    accumulatori = -accumulatori;

    % adesso ho i coefficienti dei piani non semplificati
   
    accumulatori_lavoro = accumulatori;

    % semplificazione dei piani
    for i = 1:size(B,2)
        
        mcd = gcd(sym(accumulatori_lavoro(i,:)));
        accumulatori_lavoro(i,:) = accumulatori_lavoro(i,:)/mcd;
     
    end

    accumulatori_lavoro(:,n + 1) = -accumulatori_lavoro(:,n+1);

    % stampo i piani al primale
    fprintf("\nPiani semplificati:\n")
    for i = 1:m
    
        fprintf("%d: ",B(i));

        for j = 1:n
            if accumulatori_lavoro(i,j) >= 0
                fprintf(" +%sx%d",sym(accumulatori_lavoro(i,j)),j)
            else
                fprintf(" %sx%d",sym(accumulatori_lavoro(i,j)),j)
            end
        end

        fprintf(" <= %s\n",sym(accumulatori_lavoro(i,n+1)));
    
    end



end