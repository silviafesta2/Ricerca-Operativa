% questa funzione svolge un passo del simplesso duale documentando tutti
% i passi da seguire:
% RICHIEDE LO SCRIPT PresenteCiclo.m e taglio.m
% NON FUNZIONA NEL CASO DI COSTI NULLI
% argomenti:
%       C=matrice dei costi simm (n-1*n-1) (-1 nelle caselle sotto diagonale)
%       k= k-albero
%       n1= nodo da cui partire per nodo +vicino
%       V= vettore 3*2 su ciascuna riga ci deve essere una coppia di indici
%          questi indici rappresentano le 3 var da istanziare (vengono istanziate in ordine dalla prima alla 3 riga)
% 
% USARE CON CAUTELA PERCHE' MANCANO I DOVUTI CONTROLLI SU ARGOMENTI
% SBAGLIATI(in particolar modo fare attenzione all'ordine degli indici delle var)

function B = Passi3_BB_TSPSim(C, k, n1, V)
    passi=3; %per ora funziona solo così
    %per estendere: 
    % funzione taglio estesa su p passi -> Taglio_est.m
    % funzione che genera la matrice mappa M per un numero qualunque di var
    % istanziate ->da fare
    % aggiungere p come argomento della funzione
    %togliere l'assegnamento passi=3 e mettere passi = p;

    %variabili per costruire tabella finale:
    Vect_Vs=[]; %contiene tutte le vs trovate
    Vect_Cicli=[]; %contiene tutti i cicli pseudo ottimi che ci hanno fornito le Vs
    Vect_X=[]; %contiene tutte le soluzioni pseudo ottime 

    fprintf("\nDATI:");
    display(sym(C),"C");
    fprintf("%d-ALBERO\nNODO: %d\nVAR: X%d%d, X%d%d, X%d%d\n\n",k,n1,V(1,1),V(1,2),V(2,1),V(2,2),V(3,1),V(3,2));

    fprintf("---------------------------------------------------------------\n")
    fprintf("Vi(P): %d-ALBERO",k);
    n= size(C, 1)+1; %numero di nodi

    Vs_Corrente=-1;
    Ciclo_pseudo_ottimo=zeros(n+1,1);

    A=zeros(n-1,n-1);%matrice da stampare a video con i tondini delle locazioni scelte
    CEst=[zeros(n-1,1),C;zeros(1,n)];

    
    CdiffEst = CEst([1:k-1 k+1:end], [1:k-1 k+1:end]);

    Lavoro=CdiffEst;
    Lavoro(Lavoro <= 0) = NaN;
    valori=zeros(n,1);
    VectArchi=zeros(n,2);


    i=1;
    while i<=n-2
        [v, LinInd] = min(Lavoro(:)); %trova il min in lavoro come valore e indice linearee
        [r, c] = ind2sub(size(Lavoro), LinInd); %mi da le coordinate del min
        A(r,c)=1;
        valori(i)=v;

        if r>=k
            r1=r+1;
        else 
            r1=r;
        end

        if c>=k
            c1=c+1;
        else 
            c1=c;
        end

        VectArchi(i,:)=[r1,c1];
        if(PresenteCiclo(A))
            A(r,c)=0;
            Lavoro(r,c)=NaN;
            i=i-1;
        end
        Lavoro(r,c)=NaN;
        i=i+1;
    end


    A=[A(1:k-1, :);zeros(1,n-1);A(k:end,:)]; %aggiungo riga k
    A=[A(:,1:k-1),zeros(n,1),A(:,k:end)]; %aggiungo colonna k



    Lavoro2=zeros(n,n);
    Lavoro2(k,:)=CEst(k,:);
    Lavoro2(:,k)=CEst(:,k);
    Lavoro2(Lavoro2 <= 0) = NaN;

    i=1;
    while i<=2
        [v, LinInd] = min(Lavoro2(:)); %trova il min in lavoro come valore e indice linearee
        [r, c] = ind2sub(size(Lavoro2), LinInd); %mi da le coordinate del min
        A(r,c)=1;
        Lavoro2(r,c)=NaN;
        valori(n-2+i)=v;
        VectArchi(n-2+i,:)=[r,c];
        i=i+1;
    end


    A = A(1:n-1,2:n);
    Matrice_Scelte=A.*C;
    %display(Matrice_Scelte);

    fprintf('\nArchi %d-albero con costo:\n',k);
    for h = 1:size(VectArchi, 1)
        fprintf('\t%d -- %d : %d\n', VectArchi(h, 1), VectArchi(h, 2), valori(h));
    end
    
    StrVal1=sprintf('%d+', valori(1:n-2));
    StrVal2=sprintf('%d+', valori(n-1:end));
    StrVal2(end)=[];
    fprintf('\nVi(P)=%s%s = %d\n',StrVal1,StrVal2,sum(valori)); %se sono in grado attivare la stampa a colori e mettere in rosso la prima%s
    
    fprintf("---------------------------------------------------------------\n");
    fprintf("\nCalcolo Vs(P): NODO + VICINO partendo da n=%d\n", n1);
    Lavoro3=CEst;
    Lavoro3(Lavoro3 <= 0) = NaN;
    valori=zeros(n,1);
    sequenza=zeros(n+1,1);
    A=zeros(n);

    i=1;
    ind=n1;
    while i<=n-1      

        [v1, LinInd1] = min(Lavoro3(:,ind)); 
        [v2, LinInd2] = min(Lavoro3(ind,:)); 

        if v1 <= v2 || isnan(v2)
            v = v1;
            r=LinInd1;
            c=ind;

        elseif v2 < v1 || isnan(v1)
            v = v2;
            r=ind;
            c=LinInd2;
        end

        A(r,c)=1;
        valori(i)=v;
        sequenza(i)=ind;

        Lavoro3(ind,:)=NaN;
        Lavoro3(:,ind)=NaN;

        if v1 <= v2 || isnan(v2)
            ind=r;
        elseif v2 < v1 || isnan(v1)
            ind=c;
        end

        i=i+1;
    end
    
    last=sort([ind,n1]);
    A(last(1),last(2))=1;
    A_copia=A;
    A=A(1:n-1,2:end);
    Matrice_Scelte= A.*C;
    %display(Matrice_Scelte);

    sequenza(n)=ind;
    sequenza(n+1)=n1;
    Ciclo_pseudo_ottimo=sequenza;

    valori(n)=CEst(last(1),last(2));
    Vs_Corrente=sum(valori);
    

    StrSequenza=sprintf('%d-',sequenza);
    StrSequenza(end)=[];

    StrValori=sprintf('%d+',valori);
    StrValori(end)=[];

    %inizializzazione vettori recap
    Vect_Vs=Vs_Corrente;
    Vect_Cicli=StrSequenza;

    A_copia=A_copia+1; %gli 1 diventano 2 e gli diventano 1
    A_copia=triu(A_copia); %metto tutti 0 sotto la diagonale (sulla diagonale ci saranno tutti 1 perche in a sulla diagonale c'erano tutti 0)
    A_copia=A_copia-eye(size(A_copia)); %tolgo anche la diagonale
    A_copia=A_copia';
    A_riga=A_copia(:); %lo rendo un vettore riga (serve trasporre per mantenere l'ordine corretto)

    A_riga=A_riga(A_riga~=0);
    x_soluzione=(A_riga-1)';
    Vect_X=x_soluzione;


    fprintf("SEQUENZA: %s\n\nVs(P)=%s=%d\n\n",StrSequenza,StrValori,sum(valori));


    %preparo variabili per esplorazione albero delle soluzioni
    %Vs_Corrente -> vedi riga
    %Ciclo_pseudo_ottimo -> vedi riga
    %Sol_Pseudo_ottima ->vedi riga
    X=zeros(1,passi);
    X(:,:)=-1;
    DaVisitare=zeros(1, 2^(passi+2)-2);
    DaVisitare(:,:)=-1;

    M=[1,1;1,2;2,1;2,2;2,3;2,4;3,1;3,2;3,3;3,4;3,5;3,6;3,7;3,8;-1,-1]; %mappa (se si vuole fare completo meglio creare una funzione che restituisce i e j data la posizione in DaVisitare
    
    i=1;
    j=1;
    ind=1;
    StrIstanziate=[];
    
    while ~all(DaVisitare(1:2^(passi+1)-2)==0) 
        if DaVisitare(ind)==0
            ind=ind+1;
            i=M(ind,1);
            j=M(ind,2);
            continue
        end
        

        X(1:i)=(dec2bin(j-1,i) - '0'); 
        X(i+1:passi)=-1;

        %sistemare meglio questa stringa, evitare strcat se possibile
        %RiempimentoStringa=[V,X']';
        %RiempimentoStringa=RiempimentoStringa(RiempimentoStringa())
        %->prendi tutte le colonne tranne quelle che terminano con -1
        %StrIstanziate=sprintf('X%d%d = %d,  ',RiempimentoStringa)

        StrIstanziate=[];
        for h=1:passi
            if X(h) ~= -1
            StrIstanziate=strcat(StrIstanziate,sprintf('  X%d%d = %d,',V(h,1),V(h,2),X(h)) );
            end
        end
        StrIstanziate(end)=[];
        fprintf("---------------------------------------------------------------\n");
        fprintf("\nVi(P%d%d)  -> %s \n\n",i,j,StrIstanziate);
        

        valori=zeros(n,1);
        VectArchi=zeros(n,2);
        Lavoro4=CEst;
        A=zeros(n);

        
        ind1=1;
        num1=0; %numero di valori ROSSI già assegnati (vedi appunti)
        num2=0; %numero di valori BLU già assegnati (vedi appunti)
        %preparo la matrice Lavoro e la matrice A per fare il k albero
        libero=1;%indice prima pos libera in valori e vectArchi
        while ind1<=passi
            if X(ind1) == 0
               Lavoro4(V(ind1,1),V(ind1,2))=NaN;

            elseif X(ind1) == 1
               Lavoro4(V(ind1,1),V(ind1,2))=NaN;
               A(V(ind1,1),V(ind1,2))=1;
               valori(libero)=CEst(V(ind1,1),V(ind1,2));
               VectArchi(libero,:)=[V(ind1,1),V(ind1,2)];
               libero=libero+1;

               if V(ind1,1)==k || V(ind1,2)==k %ho già preso un blu
                    num2=num2+1;
               else %ho già preso un rosso
                   num1=num1+1;
               end

            end
            ind1=ind1+1;
        end

        
        Lavoro=Lavoro4([1:k-1 k+1:end], [1:k-1 k+1:end]);
        L_riga_k=Lavoro4(k,:);
        L_col_k=Lavoro4(:,k);
        Lavoro(Lavoro <= 0) = NaN;
        A_trunc=A([1:k-1 k+1:end], [1:k-1 k+1:end]);
        A_riga_k=A(k,:);
        A_col_k=A(:,k);

        %controllo Inammissibilità ineliminabile
        % 1. NO INAMMISSIBILITA' INELIMINABILE
        % 1. VIOLATO VINCOLO DI GRADO AL NODO ?: grado > 2 -> TAGLIO
        % 1. VIOLATO VINCOLO DI GRADO AL NODO ?: grado < 2 -> TAGLIO
        % 1. SI E' CREATO UN CICLO

        Inam=-1;        % =-1 no inammissibilità, =1 gr <2, =3 gr>2, 4 CICLO
        primo_nodo_violato=-1;
        Lavoro4(Lavoro4<=0)=NaN;
        if PresenteCiclo(A)
            Inam=4;
            fprintf("1. E' presente un CICLO -> TAGLIO\n\n")
        end
        for nodo=1:n
            %controllo 1
            riga_i=Lavoro4(nodo,:);
            colonna_i=Lavoro4(:,nodo);
            quanti = sum(~isnan(riga_i));
            quanti = quanti + sum(~isnan(colonna_i)); 
            % ho sommato tutti quelli sceglibili perchè non sono NaN
            %ora sommo anche quelli che sono NaN ma li ho già scelti.
            riga_i=A(nodo,:);
            colonna_i=A(:,nodo);
            quanti= quanti+sum(colonna_i==1); 
            quanti = quanti + sum(riga_i==1);

            if quanti<2
                primo_nodo_violato=nodo;
                Inam=1;
                fprintf("1. VIOLATO VINCOLO DI GRADO AL NODO %d:\n\t grado_del_nodo_%d < 2 -> TAGLIO\n",primo_nodo_violato,primo_nodo_violato)
                break;
            end
            %controllo 2
            riga_i=A(nodo,:);
            colonna_i=A(:,nodo);
            quanti = sum(colonna_i==1);
            quanti = quanti + sum(riga_i==1);
            if quanti>2
                primo_nodo_violato=nodo;
                Inam=3;
                fprintf("1. VIOLATO VINCOLO DI GRADO AL NODO %d:\n\t grado_del_nodo_%d > 2 -> TAGLIO\n",primo_nodo_violato,primo_nodo_violato)
                break;
            end
        end
       
        if Inam==-1
           % fprintf("1. NO INAMMISSIBILITA' INELIMINABILE\n\n");
        else
            DaVisitare=taglio(DaVisitare,ind);
            DaVisitare(ind)=0;
            ind=ind+1;
            i=M(ind,1);
            j=M(ind,2);
            continue;
        end

        %calcola K albero con le matrici preparate

        i1=1+num1;
        while i1<=n-2
            [v, LinInd] = min(Lavoro(:)); %trova il min in lavoro come valore e indice linearee
            [r, c] = ind2sub(size(Lavoro), LinInd); %mi da le coordinate del min
            A_trunc(r,c)=1;
            valori(i1+num2)=v;
            if r>=k
                r1=r+1;
            else 
                r1=r;
            end

            if c>=k
                c1=c+1;
            else 
                c1=c;
            end

            VectArchi(i1+num2,:)=[r1,c1];
            if(PresenteCiclo(A_trunc))
                A_trunc(r,c)=0;
                Lavoro(r,c)=NaN;
                i1=i1-1;
            end
            Lavoro(r,c)=NaN;
            i1=i1+1;
        end

        A=[A_trunc(1:k-1, :);zeros(1,n-1);A_trunc(k:end,:)]; %aggiungo riga k (senza elemento k-esimo) della vecchia A
        A=[A(:,1:k-1),zeros(n,1),A(:,k:end)]; %aggiungo colonna k della vecchia A
        A(k,:)=A_riga_k;
        A(:,k)=A_col_k;

        Lavoro2=zeros(n,n);
        Lavoro2(k,:)=L_riga_k;
        Lavoro2(:,k)=L_col_k;
        Lavoro2(Lavoro2 <= 0) = NaN;

        i1=1+num2;
        while i1<=2
            [v, LinInd] = min(Lavoro2(:)); %trova il min in lavoro come valore e indice linearee
            [r, c] = ind2sub(size(Lavoro2), LinInd); %mi da le coordinate del min
            A(r,c)=1;
            Lavoro2(r,c)=NaN;
            valori(n-2+i1)=v;
            VectArchi(n-2+i1,:)=[r,c];
            i1=i1+1;
        end
        
        %per l'ammissibilità controllo che siano rispettati tutti i vincoli
        %di grado
        Ammissibile=true;
        StrAm="AMMISSIBILE";
        for nodo=1:n
            %controllo 1
            riga_i=A(nodo,:);
            colonna_i=A(:,nodo);
            quanti = sum(riga_i==1);
            quanti = quanti + sum(colonna_i==1);

            if quanti~=2
                Ammissibile=false;
                StrAm = "NON AMMISSIBILE";
                break;
            end

        end
        
        A_copia=A;
        A_copia2=A;
        A = A(1:n-1,2:n);
        Matrice_Scelte=A.*C;
        Vi=sum(valori);
        %display(Matrice_Scelte);


        %fprintf('Archi %d-albero con costo:\n',k);
        %for h = 1:size(VectArchi, 1)
        %    fprintf('\t%d -- %d : %d\n', VectArchi(h, 1), VectArchi(h, 2), valori(h));
        %end
    
        StrVal1=sprintf('%d+', valori(1:n-2));
        StrVal2=sprintf('%d+', valori(n-1:end));
        StrVal2(end)=[];
        fprintf('\nVi(P%d%d)=%s%s = %d\t|\tSOLUZIONE %s\n\n',i,j,StrVal1,StrVal2,Vi,StrAm); %se sono in grado attivare la stampa a colori e mettere in rosso la prima%s
    
        %Controllo 2 per il taglio
        %fprintf("2. Vi(P%d%d)>=Vs(P) ?\t ->\t %d >= %d?",i,j,Vi,Vs_Corrente);
        
        if Vi>=Vs_Corrente
            fprintf("2. Vi(P%d%d)>=Vs(P) ?\t ->\t %d >= %d?",i,j,Vi,Vs_Corrente);
            fprintf("\t->\t %s TAGLIO\n\n",char(0x2714));
            DaVisitare=taglio(DaVisitare,ind);
            DaVisitare(ind)=0;
            ind=ind+1;
            i=M(ind,1);
            j=M(ind,2);
            continue;
        else
            %fprintf("\t->\t %s NON scatta la regola 2\n\n",char(0x2716));
        end

        %Controllo 3 per il taglio e aggiorno
        
        if Vi<Vs_Corrente
            simb1 = char(0x2714); % Simbolo del visto
        else
            simb1 = char(0x2716); % Simbolo della croce
        end

        if Ammissibile
            simb2 = char(0x2714); % Simbolo del visto
        else
            simb2 = char(0x2716); % Simbolo della croce
        end 

        if Vi<Vs_Corrente && Ammissibile
            StrRis3 = ' TAGLIO';
        else
            StrRis3 = ' NON TAGLIO';
        end 

        %fprintf("   | Vi(P%d%d)<Vs(P) ?\t ->\t %d < %d? ->\t%s\n",i,j,Vi,Vs_Corrente,simb1);
        %fprintf("3. | && \n");
        %fprintf("   | SOL AMMISIBILE ?\t ->\t %s\t->\t%s\n",StrAm, simb2);

        if Vi<Vs_Corrente && Ammissibile
            fprintf("\n\t%s\t TAGLIO e AGGIORNO: \t\tNUOVA Vs(P)=%d\n\n",char(0x21D2),Vi);
            Vs_Corrente=Vi; %aggiorno

            %aggiorno valori recap
            Vect_Vs=[Vect_Vs;Vs_Corrente];

            A_copia=A_copia+1; %gli 1 diventano 2 e gli diventano 1
            A_copia=triu(A_copia); %metto tutti 0 sotto la diagonale (sulla diagonale ci saranno tutti 1 perche in a sulla diagonale c'erano tutti 0)
            A_copia=A_copia-eye(size(A_copia)); %tolgo anche la diagonale
            A_copia=A_copia';
            A_riga=A_copia(:); %lo rendo un vettore riga (serve trasporre per mantenere l'ordine corretto)
            A_riga=A_riga(A_riga~=0);
            x_soluzione=(A_riga-1)';

            Vect_X=[Vect_X;x_soluzione];
            
            
            nodi_ciclo=zeros(n+1,1);
            nodi_ciclo(1)=1;
            A_sym=A_copia2+A_copia2';
            nodo_corrente=1;
            for ind3=1:n
                nodi_ciclo(ind3+1)=find(A_sym(:,nodo_corrente)==1,1);
                A_sym(nodi_ciclo(ind3+1),nodo_corrente)=0;
                A_sym(nodo_corrente,nodi_ciclo(ind3+1))=0;
                nodo_corrente=nodi_ciclo(ind3+1);
            end 
            StrCiclo=sprintf('%d-', nodi_ciclo);
            StrCiclo(end)=[];
            
            Vect_Cicli=[Vect_Cicli;StrCiclo];

            DaVisitare=taglio(DaVisitare,ind); %taglio
            DaVisitare(ind)=0;
            ind=ind+1;
            i=M(ind,1);
            j=M(ind,2);
            continue;
        else
            fprintf("\n\t%s\t NON TAGLIO\n\n",char(0x21D2));
        end


        %aggiornamento var ciclo
        DaVisitare(ind)=0;
        ind=ind+1;
        i=M(ind,1);
        j=M(ind,2);
    end
    
    %stampa valori finali
    dim=size(Vect_Vs);

    fprintf("\n\n---------------------------------------------------------------\n");
    fprintf("TABELLA RIASSUNTIVA (più si scende più i ris sono aggiornati)\n");
    fprintf("VS_CORRENTE  \t CICLO_PESUDO_OTTIMO \t  SOLUZIONE_PSEUDO_OTTIMA\n");

    for ind6=1:dim
        x= sprintf('%d, ',Vect_X(ind6,:));
        x(end)=[];
        x(end)=[];
        fprintf("\t%d\t\t\t\t%s\t\t\t\t[%s]\n",Vect_Vs(ind6), Vect_Cicli(ind6,:),x);
    end

    ottima=false;
    if all(DaVisitare==0)
        ottima=true;
        fprintf("ALBERO COMPLETAMENTE POTATO -> la soluzione e' OTTIMA %s\n", char(0x2714));
    else
        fprintf("RIMANGONO RAMI NON POTATI -> la soluzione POTREBBE NON essere ottima %s\n",char(0x2716));
    B=x;
end
