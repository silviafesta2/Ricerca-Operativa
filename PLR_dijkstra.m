function PLR_dijkstra(ind_CU,S,n)

    % prendo soltanto archi e costi
    archi = ind_CU(:,1:2);
    costi = ind_CU(:,4);

    % inizializzo N, p e pi
    N = 1:n;
    p = ones(1,n);
    p = -p;
    p(S) = 1;
    pi = Inf(1,n);
    pi(S) = 0;
    pi_lavoro = pi;

    % adesso entro in un ciclo in cui estraggo i nodi di etichetta minima
    while ~isempty(N)

        % estrazione etichetta
        [~,ind] = min(pi_lavoro);
        N = setdiff(N,ind);
        pi_lavoro(ind) = +Inf;

        % calcolo la stella uscente come tutti gli archi che hanno come
        % origine ind
        archi_lavoro = archi(archi(:,1) == ind,:);

        for k = 1:size(archi_lavoro,1)
            % mi calcolo il costo ridotto
            i = ind;
            j = archi_lavoro(k,2);

            % intanto ottengo il costo
            c = costi(ismember(archi,[i j],'rows'));

            % se il costo ridotto e' nullo lo metto in base
            if (c + pi(i) - pi(j) < 0)
                pi(j) = c + pi(i);
                p(j) = i;

                % aggiorno anche il vettore ausiliario
                pi_lavoro(j) = c + pi(i);
            end
        end
    
    end

    fprintf("π = %s, p = %s\n",mat2str(pi),mat2str(p))

end