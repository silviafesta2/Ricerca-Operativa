% funzione che calcola la valutazione inferiore per il TSP simmetrico
function assegnamentoTSP(C)


    % ottengo il vettore dei costi dell'assegnamento
    n = size(C,1) + 1;

    c = [zeros(n,1),[C;zeros(1,n-1)]];
    c = c + c';
    c(logical(eye(n))) = 10e10;
    c = c(:);

    A = matriceAssegnamento(n,n);
    b = ones(2 * n,1);

    [x,v] = linprog(c,[],[],A,b,zeros(1,n*n),ones(1,n*n));

    fprintf("Vi(P) = %d, X = %s\n\n",v, mat2str(x));

    for k = 1:(n*n)
        if (x(k) == 1)
            i = ceil(k/n);
            j = mod(k - 1,n) + 1;
            fprintf("%d -> %d : %d\n",sym(i),sym(j),c((i-1) * n + j));
        end
    end


    function A = matriceAssegnamento(n, m)
        % Verifica che le dimensioni siano positive
        if n <= 0 || m <= 0
            error('Le dimensioni devono essere numeri positivi.');
        end
        
        % Inizializza la matrice di trasporto A
        A = zeros((n + m), n * m);
        
        % Riempie la matrice A con vincoli di trasporto da luoghi di produzione a luoghi di raccolta
        for i = 1:n
            for j = 1:m
                row = i;
                col = (i - 1) * m + j;
                A(row, col) = 1;
                
                row = n + j;
                A(row, col) = 1;
            end
        end
    end

end