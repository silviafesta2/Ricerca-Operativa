% questa funzione produce la matrice di incidenza di una rete
% argomenti:
%       A ->        matrice degli archi (i,j) passata come matrice [i1 j1; i2 j2; ...]
% risultato:
%       produce la matrice di incidenza E della rete
function incidence_matrix = matriceIncidenza(edges)
    num_nodes = max(max(edges)); % Determina il numero di nodi
    
    num_edges = size(edges, 1); % Determina il numero di archi
    
    incidence_matrix = zeros(num_nodes, num_edges); % Inizializza la matrice di incidenza
    
    for k = 1:num_edges
        i = edges(k, 1); % Nodo iniziale dell'arco
        j = edges(k, 2); % Nodo finale dell'arco
        
        incidence_matrix(i, k) = -1; % Decrementa il nodo iniziale
        incidence_matrix(j, k) = 1;  % Incrementa il nodo finale
    end
end
