% argomenti: 
%       Aeq ->    matrice del problema duale standard per cui dobbiamo
%               trovare un piano di taglio
%       indiciPositivi ->    insieme degli indici di base passati come vettore riga
% risultato:    stampa la matrice Aeq tilde e le parti frazionarie, in piu
%               vengono stampati anche tutti i possibili piani di taglio
function gom = calcolaGomory(Aeq, indiciPositivi)
    % Verifica che indiciPositivi sia un sottoinsieme degli indici delle colonne di Aeq
    if any(indiciPositivi > size(Aeq, 2))
        error('indiciPositivi contiene indici fuori dai limiti delle colonne di Aeq.');
    end
    
    % Calcola N come l'insieme complementare di indiciPositivi rispetto agli indici delle colonne di Aeq
    colInd = 1:size(Aeq, 2);
    N = setdiff(colInd, indiciPositivi);
    
    % Calcola Ab come le colonne di Aeq corrispondenti agli indici in indiciPositivi
    Ab = Aeq(:, indiciPositivi);
    display(sym(Ab),"Ab");
    
    % Calcola bA come l'inversa di Ab
    bA = inv(sym(Ab));
    display(sym(bA),"Ab^(-1)");
    
    % Calcola An come le colonne di Aeq corrispondenti agli indici in N
    An = Aeq(:, N);
    display(sym(An),"An");

    % Calcola la matrice aTilde come prodotto tra bA e An
    aTilde = sym(bA) * sym(An);
    
    % Stampa la matrice aTilde
    disp('Matrice aTilde:');
    disp(sym(aTilde),"aTilde");
    
    % Calcola la matrice Aeq~ con le parti frazionarie di aTilde
    A_tilde_frac = sym(aTilde) - floor(sym(aTilde));
    
    % Stampa la matrice Aeq~
    disp('Matrice Aeq~ (parti frazionarie):');
    disp(sym(A_tilde_frac),"Parti Frazionarie");

    gom = A_tilde_frac;
end