function knapsackBranchAndBound(c, A, b)

    % porto c ed A a vettore di riga
    [n,~] = size(c);
    if (n > 1)
        c = c';
    end

    [n,~] = size(A);
    if (n > 1)
        A = A';
    end

    n = length(c);

    % calcolo il vettore dei rendimenti
    r =  c ./ A;

    % ordino il vettore dei rendimenti ed ottengo una mappa dagli elementi
    % non ordinati a quelli ordinati
    r_lavoro = sort(r,'descend');
    m = [];
    for i = 1:n
        ind = find(r == r_lavoro(i));
        m = [m,ind];
    end
    

    % Definisco il root node
    root.x = zeros(n,1);
    root.i = 0;
    root.j = 0;

    % calcolo Vi e Vs
    [x,Vs,Vi,~,xvi] = vsKnapsackBrancAndBound(root,c,A,b,0,0);

    % in xvi e' contenuta l'approssimazione della nostra soluzione ottima

    fprintf("Vs(P) = %d, X = %s\nVi(P) = %d, X = %s\n\n",Vs,char(sym(x)),Vi,mat2str(xvi))

    % inserisco i primi nodi all'interno della coda
    x_frac = x - floor(x);
    x_frac(x_frac ~= 0) = 1;
    x_frac_lavoro = -x_frac;


    % ogni nodo ha un id univoco calcolato come 2^i + j
    % la visita e' fatta per ampiezza -> id crescente

    node1.i = 1;
    node1.j = 1;
    node1.x = x_frac_lavoro;
    node1.id = 2^node1.i + node1.j;
    
    node2.i = 1;
    node2.j = 2;
    node2.x = x_frac;
    node2.id = 2^node2.i + node2.j;

    Q = [node1, node2];

    while ~isempty(Q)
        % Ordina la coda per id crescente
        [~, idx] = sort([Q.id]);
        Q = Q(idx);

        % Rimuovo l'indice con id minore
        node = Q(1);
        Q(1) = [];

        % ottengo la vs per il nodo corrente
        [x,Vs,Vi,rule] = vsKnapsackBrancAndBound(node,c,A,b,Vi,Vs);

        % quando scatta la regola 3 devo assegnare xvi ad x
        if (rule == 3)
            xvi = x;
        end
        
        % se non e' scattata alcuna regola di taglio devo aggiungere
        % 2 nodi figli
        if (rule ~= 0)
            continue;
        end

        % nodo di destra
        node1.i = node.i + 1;
        node1.j = node.j * 2;
        x_frac = x - floor(x);
        x_frac(x_frac ~= 0) = 1;
        node1.x = node.x + x_frac;
        node1.id = 2^node1.i + node1.j;

        x_frac_lavoro = -x_frac;

        % nodo di sinistra
        node2.i = node.i + 1;
        node2.j = node.j * 2 - 1;
        node2.x = node.x + x_frac_lavoro;
        node2.id = 2^node2.i + node2.j;

        Q = [node1,node2,Q];

    end

    fprintf("\nSoluzione ottima trovata: X = %s\n",mat2str(xvi))
        
end