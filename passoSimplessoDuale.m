% questa funzione svolge un passo del simplesso duale documentando tutti
% i passi da seguire:
% argomenti:
%       come linprog (ma c vettore colonna) (c,A,b) per risolvere 
%       {min b^t*x, Ax=c, x>=0
%       B vettore riga con gli indici della base
% USARE CON CAUTELA PERCHE' MANCANO I DOVUTI CONTROLLI SU ARGOMENTI SBAGLIATI
% attenzione nel codice y e yb li ho chiamati ytot e y non fare confusione
function B = passoSimplessoDuale(c, A, b, B)
    fprintf("DATI:");
    display(sym(A),"A");
    display(sym(b),"b");
    display(sym(transpose(c)),"c^t");
    display(sym(B),"B");

    N = setdiff(1:size(A, 1), B);
    display(sym(N),"N");

    fprintf("Calcolo Ab, Ab^(-1), W e bB");

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
    bn = b(N);
    An = A(N,:);
    
    display(sym(Ab),"Ab");
    display(sym(bA),"Ab^(-1)");
    display(sym(W),"W=-Ab^(-1)");
    display(sym(bb),"bB");

    % trovo l'indice entrante
    fprintf("Calcolo y,v,x e TEST OTTIMALITA'");
    x = bA * bb;
    y = transpose(c) * bA;

    ytot = zeros(1, size(A, 1));
    for i = 1:length(B)
        ytot(B(i)) = y(i);
    end

    display(sym(y),"yb = c * Ab^(-1)");
    display(sym(ytot),"y=(yb,yn)");
    display(sym(b' * ytot'),'v = b^t * y');
    display(sym(x),"x = Ab^(-1) * bb");

    fprintf("\nAn*x <= b?\n");
    display(sym(An*x), "An*x");
    display(sym(bn), "bn");
    Anx=An*x;

    IndN=zeros(length(N),1);
    

    for i= 1:length(N)
        if ( abs(Anx(i) - bn(i)) < 1e-10 ) %serve una certa tolleranza per l'=per colpa delle approssimazioni
          fprintf("= ");
          continue
      end
        
      if ( Anx(i)<(bn(i)) )
          fprintf("%s ",char(0x2713));
      end
      if ( Anx(i)>(bn(i)) )
          IndN(i)=-1;
          fprintf("X ");
      end
    end
    fprintf("\n");

    if all(IndN >= 0)
        fprintf('x AMMISSIBILE -> B ottima e X ottima \n');
        return;
    else

        Ind = zeros(size(A, 1),1);
        for i = 1:length(N)
            Ind(N(i)) = IndN(i);
        end


        NegInd = find(Ind < 0); %array degli indici negativi
        NegIndStr = sprintf('%d,', NegInd);
        NegIndStr(end) = []; %stringa con tutti gli indici negativi

        fprintf('X NON AMMISSIBILE -> B NON OTTIMA\n\n');
        fprintf("Calcolo INDICE ENTRANTE: k\n");
        k = find(IndN < 0, 1); % Trova la PRIMA componente negativa di IndN
        knorm = N(k);
        fprintf("k = min{i%sN t.c. Ai*x>bi} = min{%s} = %d\n",char(0x2208), NegIndStr, knorm);
        fprintf("k = %d\n\n", knorm);
    end
    
    fprintf("Calcolo Ak * W^i %si%sB",char(0x2200),char(0x2208));
    Ak = A(knorm, :); % prendo la riga
    
    r = [];
    rapporti = [];

    
    display(sym(Ak), sprintf("Ak= A^%d",knorm));

    for i = 1:length(B)
        val = Ak * W(:,i);
        fprintf("A%d * w^%d = %s\n",knorm, B(i),sym(val));
        r = [r,val];
    end

    if all(r >= 0)
        fprintf("PRODOTTI SCALARI TUTTI POSITIVI o nulli:\nAk * W^i > 0 %si%sB -> v(P)=-%s\n",char(0x2200),char(0x2208),char(0x221E));
        return;
    end

    fprintf("\nCalcolo ri %si%sB t.c. Ak * W^i < 0\n", char(0x2200), char(0x2208));

    for i = 1:length(r)
        if (abs(r(i)) < 1e-10) %dovrebbe essere == 0 ma bisogna aumentare un pò la tolleranza per gestire piccoli arrotondamenti sul calcolo di W
            fprintf("Prodotto scalare nullo (vincolo parallelo)\n");
            rapporti = [rapporti;-1];
            continue;
        end

        if (r(i) > 0)
            fprintf("r%d = Prodotto scalare positivo\n",N(i));
            rapporti = [rapporti;-1];
            continue;
        end

        rapporto = -y(i) / (r(i));
        fprintf("r%d = %s\n",B(i),sym(rapporto));
        rapporti = [rapporti;rapporto];
    end

    fprintf("\nCalcolo INDICE USCENTE: h\n");
    indiciNonNegativi = find(rapporti >= 0);
    elementiNonNegativi = rapporti(indiciNonNegativi);
    [~, indiceRelativo] = min(elementiNonNegativi);
    indice = indiciNonNegativi(indiceRelativo);

    StrRapporti=sprintf('%s,', sym(elementiNonNegativi));
    StrRapporti(end)=[];
    fprintf("h t.c. rh= min{ri}=min{%s}=%s=r%d\nh=%d\n",StrRapporti,sym(rapporti(indice)),B(indice),B(indice));
    % creo la nuova base
    B = setdiff(B,B(indice));
    B = [B, knorm];
    B = sort(B);

    fprintf('\nCalcolo Nuova BASE:\n');
    display(sym(B),"B");
end
