% questa funzione svolge un passo del simplesso primale documentando tutti
% i passi da seguire:
% argomenti:
%       come linprog (ma c vettore colonna) (c,A,b) per risolvere {max c*x, Ax<=b
%       B vettore riga con gli indici della base
% USARE CON CAUTELA PERCHE' MANCANO I DOVUTI CONTROLLI SU ARGOMENTI SBAGLIATI
% attenzione nel codice y e yb li ho chiamati ytot e y non fare confusione
function B = passoSimplessoPrimale(c, A, b, B)
    fprintf("DATI:");
    display(sym(A),"A");
    display(sym(b),"b");
    display(sym(transpose(c)),"c^t");
    display(sym(B),"B");

    N = setdiff(1:size(A, 1), B);
    display(sym(N),"N");

    fprintf("Calcolo Ab, Ab^(-1), W e bB\n");

    Ab = A(B,:);
    
    % faccio un controllo che Ab sia non singolare
    if det(Ab) == 0
        fprintf("Ab non e' invertibile, fornire una base valida\n");
        return;
    end
    
    % calcolo le matrici che mi serviranno
    bA = Ab^(-1);
    W = -bA;
    bb = b(B);
    
    display(sym(Ab),"Ab");
    display(sym(bA),"Ab^(-1)");
    display(sym(W),"W=-Ab^(-1)");
    display(sym(bb),"bB");

    % trovo l'indice uscente
    fprintf("Calcolo x,v,y e TEST OTTIMALITA'");
    x = bA * bb;
    y = transpose(c) * bA;

    ytot = zeros(1, size(A, 1));
    for i = 1:length(B)
        ytot(B(i)) = y(i);
    end

    display(sym(x),"x = Ab^(-1) * bb");
    display(sym(c' * x),'v = c^t * x');
    display(sym(y),"yb = c * Ab^(-1)");
    display(sym(ytot),"y=(yb,yn)");

    if all(y >= 0)
        fprintf('y AMMISSIBILE (yi >= 0 %s i) -> B ottima e X ottima \n', char(0x2200)); %il carattere è la codifica unicode del simbolo PER OGNI
        return;
    else
        NegInd = find(ytot < 0); %array degli indici negativi
        NegIndStr = sprintf('%d,', NegInd);
        NegIndStr(end) = []; %stringa con tutti gli indici negativi

        fprintf('y NON AMMISSIBILE (y%s <0) -> B NON OTTIMA\n\n', NegIndStr);
        fprintf("Calcolo INDICE USCENTE: h\n");
        h = find(y < 0, 1); % Trova la PRIMA componente negativa di y
        hnorm = B(h);
        fprintf("h = min{i%sB t.c. yi<0} = min{%s} = %d\n",char(0x2208), NegIndStr, hnorm);
        fprintf("h = %d\n\n", hnorm);
    end
    
    fprintf("Calcolo Ai * W^h %si%sN\n",char(0x2200),char(0x2208));
    wh = W(:, h); % prendo la colonna corretta di W
    
    r = [];
    rapporti = [];
    k = 0;
    
    display(sym(wh), sprintf("W^h= W^%d",hnorm));

    for i = 1:length(N)
        val = A(N(i),:) * wh;
        fprintf("A%d * w^%d = %s\n",N(i),hnorm,sym(val));
        r = [r,val];
    end

    if all(r <= 0)
        fprintf("PRODOTTI SCALARI TUTTI NEGATIVI o nulli:\nAi * W^h < 0 %si%sN ->v(P)=+%s\n",char(0x2200),char(0x2208),char(0x221E));
        return;
    end

    fprintf("\nCalcolo ri %si%sN t.c. Ai*W^h > 0\n", char(0x2200), char(0x2208));

    for i = 1:length(r)

        if (abs(r(i)) < 1e-10)
            fprintf("Prodotto scalare nullo (vincolo parallelo)\n");
            rapporti = [rapporti;-1];
            continue;
        end

        if (r(i) < 0)
            fprintf("r%d = Prodotto scalare negativo\n",N(i));
            rapporti = [rapporti;-1];
            continue;
        end

        rapporto = (b(N(i)) - A(N(i),:) * x) / (r(i));
        fprintf("r%d = %s\n",N(i),sym(rapporto));
        rapporti = [rapporti;rapporto];
    end

    fprintf("\nCalcolo INDICE ENTRANTE: k\n");
    indiciNonNegativi = find(rapporti >= 0);
    elementiNonNegativi = rapporti(indiciNonNegativi);
    [~, indiceRelativo] = min(elementiNonNegativi);
    indice = indiciNonNegativi(indiceRelativo);

    StrRapporti=sprintf('%s,', sym(elementiNonNegativi));
    StrRapporti(end)=[];
    fprintf("k t.c. rk= min{ri}=min{%s}=%s=r%d\nK=%d\n",StrRapporti,sym(rapporti(indice)),N(indice),N(indice))
    % creo la nuova base
    B = setdiff(B,hnorm);
    B = [B, N(indice)];
    B = sort(B);
    fprintf('\nCalcolo Nuova BASE:\n');
    display(B);

end
