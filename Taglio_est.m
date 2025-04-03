function ris= Taglio_est(b,ind,passi)
    N=3:2^(passi+2)-2;
    liv = ceil(log2(11+2))-1; %livello dell'indice
    liv_figli= passi+1-liv; % numero di livelli di figli del nodo
    B=[];
    if liv_figli~=0 %non dovrebbe succedere mai nel BB ma per evitare errori lo mettiamo
        B=[ind*2+1,ind*2+2];
        for i=2:liv_figli+1 %se liv figli minore di 2 questa parte non viene eseguita, non da errori
          F=sort([(B(end-i+1:end)*2)+1,(B(end-i+1:end)*2)+2]);
          B=[B,F]
        end
    end
    ris=b;
    ris(B)=0;
end
