% -1 -> variabile istanziata a 0
% 1 -> variabile istanziata a 1
% 0 -> variabile non istanziata

function [x,Vs,Vi,rule,xvi] = vsKnapsackBrancAndBound(node,c,A,b,Vi,Vs)

    rule = 0;
    xvi = 0;

    n = length(c);

    i = node.i;
    j = node.j;
    x = node.x;
    bres = b;
    x_lavoro = x;

    if (i ~= 0 && j ~= 0)
        % se non sono al root node stampo il problema e le variabili istanziate
        fprintf("P%d,%d:",i,j);
        
        for k = 1:n
            if (x(k) == -1)
                fprintf(" x%d = 0",k)
            elseif (x(k) == 1)
                fprintf(" x%d = 1",k)
            end
        end
        fprintf("\n");
    end

    % calcolo il vettore dei rendimenti
    r =  c ./ A;

    % ordino il vettore dei rendimenti ed ottengo una mappa dagli elementi
    % non ordinati a quelli ordinati
    r_lavoro = sort(r,'descend');
    m = [];
    for k = 1:n
        ind = find(r == r_lavoro(k));
        m = [m,ind];
    end



    % sistemo x istanziandolo
    for k = 1:n

        % variabile istanziata a 0
        if (x(k)) == -1
            x(k) = 0;
            continue;
        end

        % variabile istanziata ad 1
        if x(k) == 1
            bres = bres - A(k);

            if (bres <= 0)
                break
            end

            continue
        end

    end

    % controllo inammisibilita' ineliminabile
    if bres < 0
        fprintf("Scatta la prima regola: P%d%d = ∅\n",i,j)
        rule = 1;
        return;
    end

    % adesso saturo finche' posso
    for k = 1:n

        % ho riempito lo zaino
        if bres <= 0
            break
        end

        % se e' una variabile gia' istanziata devo continuare
        if (x_lavoro(m(k)) ~= 0)
            continue;
        end

        % devo controllare quanto ne posso inserire
        if A(m(k)) > bres
            % devo inserirlo frazionario
            x(m(k)) = bres/A(1,m(k));
            bres = 0;
            break
        else
            % lo posso inserire ad 1
            x(m(k)) = 1;
            bres = bres - A(1,m(k));
        end
    end

    % calcolo Vs
    Vs = floor(c * x);
    
    if (i == 0 && j == 0)
        % se sono al root node devo calcolare Vi
        
        x_lavoro = zeros(n,1);
        bres = b;
        for k = 1:n

            if bres <= 0
                break
            end
        
            % non c'e' spazio per questo bene
            if (A(m(k)) > bres)
                continue
            else
                x_lavoro(m(k)) = 1;
                bres = bres - A(m(k));
            end

        end

        xvi = x_lavoro;
        Vi = c * xvi;
        

    else
        % altrimenti devo stampare e verificare l'applicazione delle regole
        fprintf("Xrc = %s, Vi = %d, Vs = %d",char(sym(x)),Vi,Vs)

        if (Vs <= Vi)
            fprintf(", Scatta la regola 2: Vs(P%d%d) <= Vi(P)\n\n",i,j)
            rule = 2;
        else
            x_frac = x - floor(x);
            if all(x_frac == 0) % allora x e' a componenti intere
                rule = 3;
                Vi = Vs;
                fprintf(", Scatta la regola 3: Vs(P%d%d) > Vi(P) e soluzione ammissibile\n\n",i,j);
            end
        end

        if (rule == 0)
            fprintf(", Non scatta alcuna regola\n\n")
        end
    end

end