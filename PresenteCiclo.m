function ciclo = PresenteCiclo(A)
    % in input ho una matrice di 0 o 1
    % devo capire se e' presente un ciclo
    %OUTPUT:1 c'e ciclo 0 non c'e ciclo
    G=graph(A+A');
    ciclo=hascycles(G);
end