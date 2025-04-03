% passati due numeri n ed m, (n ingegneri ed m lavori ad esempio)
% questa funzione produce come risultato la relativa matrice
% per il problema dell'assegnamento o del traporto
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

