function x = trova_ciclo(A)
    G = graph(A);
    x = allcycles(G);
    x = x{1};
end