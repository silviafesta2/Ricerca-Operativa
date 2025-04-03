function punti_ord = ordina_antiorario(p, r)
    % Calcola gli angoli rispetto a r
    angoli = atan2(p(:,2) - r(2), p(:,1) - r(1));
    
    % Calcola gli angoli in radianti e trasforma nell'intervallo [0, 2*pi)
    angoli = mod(angoli, 2*pi);
    
    % Ordina i punti in base agli angoli
    [~, ind] = sort(angoli);
    
    % Riordina i punti in base agli angoli ordinati
    punti_ord = p(ind, :);
end