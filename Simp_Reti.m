function Simp_Reti(ind_cu,num_nodi,b)
    %ARGOMENTI

    %ind_cu: 5 componenti
    % 1: nodo di PARTENZA dell'arco
    % 2: nodo di ARRIVO dell'arco
    % 3: il loro tipo: 1=T, 2=U, 3=L
    % 4: COSTO dell'arco
    % 5: CAPACITA' dell'arco

    %num_nodi = numero di nodi della rete

    %b= vettore bilanci ai nodi (PER INTERO)

    if sum(b)~=0
        fprintf("rete non bilanciata, BILANCI ERRATI\n");
        return;
    end
    ind_cu=sortrows(ind_cu,1:2);
    
    indici_x=ind_cu(:,1:3);
    cu=ind_cu(:,4:5);

    L = indici_x(indici_x(:,3) == 3, :);
    L=L(:,1:2);
    %display(L);
    
    T = indici_x(indici_x(:,3) == 1, :);
    T=T(:,1:2);
    %display(T);

    U = indici_x(indici_x(:,3) == 2, :);
    U=U(:,1:2);
    %display(U);

    fprintf("\nDATI:\n");
    
    %display T
    fprintf("\tT={")
    Str=sprintf('(%d,%d) ',T');
    Str(end)=[];
    fprintf("%s",Str);
    fprintf("}\n");

    %display U
    fprintf("\tU={")
    Str=sprintf('(%d,%d) ',U');
    Str(end)=[];
    fprintf("%s",Str);
    fprintf("}\n");
    
    %display L
    fprintf("\tL={")
    Str=sprintf('(%d,%d) ',L');
    Str(end)=[];
    fprintf("%s",Str);
    fprintf("}\n");

    [num_archi,~]=size(indici_x);
    %calcolo matrice di incidenza della rete
    Et= zeros(num_archi,num_nodi);
    for i=1:num_archi;
        Et(i,indici_x(i,1))=-1;
        Et(i,indici_x(i,2))=1;
    end
    %display(Et);

    %indici_x=[];
    %for i=1:size(Et)
    %    inizio=find(Et(i,:)==-1);
    %    fine=find(Et(i,:)==1);
    %    indici_x(i,:)=[inizio,fine];
    %end
    %display(indici_x);

    c=cu(:,1);
    u=cu(:,2);
    
    fprintf("\nCALCOLO x e %s:\n",char(960));
    x=zeros(1,num_archi);
    x_bin=zeros(1,num_archi);

    b_lavoro=b;
    
    for i=1:num_archi
        if indici_x(i,3)== 2 %arco appartente a U
            x(i)=u(i);
            x_bin(i)=1;
            b_lavoro(indici_x(i,1))=b_lavoro(indici_x(i,1))+u(i); %aggiorno bilancio nodo partenza
            b_lavoro(indici_x(i,2))=b_lavoro(indici_x(i,2))-u(i); %aggiorno bilancio nodo arrivo
        end
        if indici_x(i,3)== 3 %arco appartenente a L
            x_bin(i)=1;
        end 
    end

    %display(x);
    %display(x_bin);
    
    
    lavoro_archi=[indici_x,cu,(1:num_archi)']; 
    %display(lavoro_archi);
    lavoro_archi= lavoro_archi(lavoro_archi(:,3) == 1, :); %prendo solo gli archi in T

    %5componenti nodo partenza, nodo arrivo, costo arco, capacità, numero della relativa componente di x
    lavoro_archi= lavoro_archi(:, [1:2, 4:end]);

    %display(lavoro_archi);
    while ~all(x_bin==1)
        n_ue=zeros(num_nodi,2); %numero di archi entranti e uscenti dal nodo

        for i=1:num_nodi
            n_ue(i,:)=[sum(lavoro_archi(:,1)==i),sum(lavoro_archi(:,2)==i)];

            if n_ue(i,1)==1 && n_ue(i,2)==0  %c'e' solo un arco uscente
                ind=find(lavoro_archi(:,1)==i); %trovo dentro la tabella lavoro l'indice di riga dell'arco uscente
                nodo_arrivo=lavoro_archi(ind,2);
                
                x(lavoro_archi(ind,5))=-b_lavoro(i);
                x_bin(lavoro_archi(ind,5))=1;
                b_lavoro(nodo_arrivo)=b_lavoro(nodo_arrivo)+b_lavoro(i);
                b_lavoro(i)=0;
                %tolgo l'arco processato
                lavoro_archi=lavoro_archi([1:ind-1 ind+1:end],:);
            end

            if n_ue(i,1)==0 && n_ue(i,2)==1  %c'e' solo un arco entrante
                ind=find(lavoro_archi(:,2)==i); %trovo dentro la tabella lavoro l'indice di riga dell'arco uscente
                nodo_partenza=lavoro_archi(ind,1); %trovo nodo di partenza
                
                x(lavoro_archi(ind,5))=b_lavoro(i);
                x_bin(lavoro_archi(ind,5))=1; %x assegnata
                b_lavoro(nodo_partenza)=b_lavoro(nodo_partenza)+b_lavoro(i);
                b_lavoro(i)=0;
                lavoro_archi=lavoro_archi([1:ind-1 ind+1:end],:);
            end
        end        
    end
    

    lavoro_archi=[indici_x,cu,(1:num_archi)'];
    lavoro_archi= lavoro_archi(lavoro_archi(:,3) == 1, :); %prendo solo gli archi in T
    lavoro_archi= lavoro_archi(:, [1:2, 4:end]);
    ind_T=lavoro_archi(:,end);
    Eb=Et(ind_T,2:end);
    cb=c(ind_T);
    pi=[0,((Eb^(-1))*cb)'];

    %calcolo valore funzione obiettivo
    costo=c'*x';
    %display x
    fprintf("\tX=[");
    str=sprintf('%d,',x);
    str(end)=']';
    fprintf("%s\t->\tcosto=%d\n",str,costo);

    %display pi
    fprintf("\t%s=[",char(960));
    str=sprintf('%d,',pi);
    str(end)=']';
    fprintf("%s\n",str);    
    
    fprintf("\nTEST OTTIMALITA' CON BELLMAN\n")
    lavoro_archi=[indici_x,cu,(1:num_archi)'];
    lavoro_archi= lavoro_archi(lavoro_archi(:,3) ~= 1, :); %prendo solo gli archi in L e U, SONO ORDINATI perchè provengono da lavoro archi
    num_costi_rid=size(L(:,1))+size(U(:,1));
    num_costi_rid(end)=[];
    c_rid=zeros(num_costi_rid,1);
    c_rid=lavoro_archi(:,4)+pi(lavoro_archi(:,1))'-pi(lavoro_archi(:,2))';

    violano_Bellman=[];
    for i=1:num_costi_rid
        tipo='L';
        str='>= 0';
        simbolo=char(0x2714); %simbolo visto


        if lavoro_archi(i,3)==2
            tipo='U';
            str ='<= 0';
            if c_rid(i) > 0
                simbolo=char(0x2716); %croce
                violano_Bellman=[violano_Bellman;lavoro_archi(i,[1:2 6])];
            end

        elseif c_rid(i) < 0
            simbolo=char(0x2716);
            violano_Bellman=[violano_Bellman;lavoro_archi(i,[1:2 6])];
        end

        fprintf('\t%s c%d%d = %d%+d%+d\t=\t%d\t %s %s\n',tipo,lavoro_archi(i,1),lavoro_archi(i,2),lavoro_archi(i,4),pi(lavoro_archi(i,1)),-pi(lavoro_archi(i,2)),c_rid(i),str,simbolo);
    end
    
    display(violano_Bellman);
    if(isempty(violano_Bellman))
        fprintf("Siamo all'OTTIMO\n")
        return;
    end

    arco_entrante=violano_Bellman(1,:);

    fprintf("\n\tArco ENTRANTE -> \t(%d,%d)\n",violano_Bellman(1,1:2));
    
    lavoro_archi=[indici_x,cu,(1:num_archi)'];
    riga_arco_entrante=lavoro_archi(arco_entrante(3),:);
    provenienza_arco_entrante=riga_arco_entrante(3);
    lavoro_archi= lavoro_archi(lavoro_archi(:,3) == 1, :); %prendo solo gli archi in T e l'arco entrante
    lavoro_archi=[lavoro_archi;riga_arco_entrante];

    A=zeros(num_nodi);
    A_sym=A;
    for i=1:num_nodi %nell'albero ci sono num_nodi-1 archi + 1 che è quello entrante
        A(lavoro_archi(i,1),lavoro_archi(i,2))=1;
        A_sym(lavoro_archi(i,1),lavoro_archi(i,2))=1;
        A_sym(lavoro_archi(i,2),lavoro_archi(i,1))=1;
    end
    
    C=trova_ciclo(A_sym);
    Presunti_Archi_ciclo=[C',circshift(C, -1)'];
    [num_archi_ciclo,~]=size(Presunti_Archi_ciclo);
    C_piu=[];
    C_meno=[];
    C1=[];
    C2=[];
    for i=1:num_archi_ciclo
        if A(Presunti_Archi_ciclo(i,1),Presunti_Archi_ciclo(i,2))==1
            C1=[C1;Presunti_Archi_ciclo(i,:)];
        else
            arco_vero=[Presunti_Archi_ciclo(i,2),Presunti_Archi_ciclo(i,1)];
            C2=[C2;arco_vero];
        end
    end

    if ~isempty(C1)
        C1=sortrows(C1, 1:2);
    end
    if ~isempty(C2)
        C2=sortrows(C2, 1:2);
    end

    if ( (~isempty(C1) && any(all(C1 == arco_entrante(1:2), 2)) && provenienza_arco_entrante==3) ) || (~isempty(C2) &&( any(all(C2 == arco_entrante(1:2), 2)) && provenienza_arco_entrante==2) )
        C_piu = C1;
        C_meno= C2;
    else
        C_piu = C2;   
        C_meno = C1;
    end

    fprintf("\nTROVO C+ e C-\n")
    fprintf("\tC+ = { ");
    fprintf("(%d,%d) ",C_piu');
    fprintf("}\n");

    fprintf("\tC- = { ");
    fprintf("(%d,%d) ",C_meno');
    fprintf("}\n");
    
    %display(C_piu);
    %display(C_meno);

    fprintf("\nCALCOLO %s+, %s- e %s\n",char(0x03B8),char(0x03B8),char(0x03B8));
    %calcolo theta + e l'arco che lo genera
    indici_righe_selezionate=[];
    indici_piu=[];
    if ~isempty(C_piu)
        indici_righe_selezionate = ismember(lavoro_archi(:, 1:2), C_piu, 'rows');
        lavoro_archi_1 = lavoro_archi(indici_righe_selezionate, :);
        indici_piu=lavoro_archi_1(:,6);
    end
    %vettore di possibili valori per theta+ e indici relativi
    Val_ind=[u(indici_piu)-x(indici_piu)',indici_piu];
    if ~isempty(Val_ind)
        Val_ind=sortrows(Val_ind, 2);%bland
    end

    valore_minimo=+Inf;
    ind=0;
    riga_theta_piu=[+Inf,0];
    if ~isempty(C_piu)
        valore_minimo=min(Val_ind(:,1));
        ind=find(Val_ind(:, 1) == valore_minimo, 1);
        riga_theta_piu=Val_ind(ind,:);
    end

    if ~isempty(C_piu)
        fprintf("\t%s+ = min{",char(0x03B8));
        str=sprintf('%d,',Val_ind(:,1));
        str(end)=[];
        fprintf("%s}=%d -> arco relativo:(%d,%d)\n",str,riga_theta_piu(1),indici_x(riga_theta_piu(2),1),indici_x(riga_theta_piu(2),2))
    else
        fprintf("\t%s+ = +∞\n",char(0x03B8));
    end
    %calcolo theta - e l'arco che lo genera
    indici_righe_selezionate = ismember(lavoro_archi(:, 1:2), C_meno, 'rows');
    lavoro_archi_2 = lavoro_archi(indici_righe_selezionate, :);
    indici_meno=lavoro_archi_2(:,6);
    Val_ind=[x(indici_meno)',indici_meno];
    if ~isempty(Val_ind)
    Val_ind=sortrows(Val_ind, 2);%BLAND
    end
    valore_minimo=min(Val_ind(:,1));
    ind=find(Val_ind(:, 1) == valore_minimo, 1);
    riga_theta_meno=Val_ind(ind,:);
    
    fprintf("\t%s+ = min{",char(0x03B8));
    str=sprintf('%d,',Val_ind(:,1));
    str(end)=[];
    fprintf("%s}=%d -> arco relativo:(%d,%d)\n",str,riga_theta_meno(1),indici_x(riga_theta_meno(2),1),indici_x(riga_theta_meno(2),2))
    
    %calcolo theta
    mat_theta=sortrows([riga_theta_meno;riga_theta_piu],2); %BLAND
    valore_minimo=min(mat_theta(:,1));
    ind=find(mat_theta(:, 1) == valore_minimo, 1);
    riga_theta=mat_theta(ind,:);

    fprintf("\t%s = min{%s+,%s-} = min{%d,%d}=%d  -> arco relativo: (%d,%d)\n\n",char(0x03B8),char(0x03B8),char(0x03B8),riga_theta_piu(1),riga_theta_meno(1),riga_theta(1),indici_x(riga_theta(2),1),indici_x(riga_theta(2),2))

    fprintf("\tARCO USCENTE ->\t(%d,%d)\n",indici_x(riga_theta(2),1),indici_x(riga_theta(2),2));

    %update
    fprintf("\nUPDATE FLUSSO,TRIPARTIZIONE, COSTO\n");
    %aggiorno x
    x(indici_meno)=x(indici_meno)-riga_theta(1);
    x(indici_piu)=x(indici_piu)+riga_theta(1);
    %aggiorno v
    costo=c'*x';
    %aggiorno tripartizione
    T=sortrows([T;riga_arco_entrante(1:2)],1:2);
    T(ismember(T, [indici_x(riga_theta(2),1),indici_x(riga_theta(2),2)], 'rows'), :) = [];
    %cancello l'arco entrante ovunque esso si trovi (una delle due
    %istruzioni non avrà nessun effetto)
    L(ismember(L,riga_arco_entrante(1:2),'rows'),:)=[];
    U(ismember(U,riga_arco_entrante(1:2),'rows'),:)=[];
    %inserisco l'arco uscente nell'insieme giusto a seconda se proviene da
    %C+ o C-


    if ismember(riga_theta(2),indici_piu) %l'indice dell'arco uscente appareneva c+
        U=[U;indici_x(riga_theta(2), 1:2)]; %lo mando in U
    else
        L=[L;indici_x(riga_theta(2), 1:2)]; %sennò lo mando in L
    end

    fprintf("\tFLUSSO:\tX=[");
    str=sprintf('%d,',x);
    str(end)=']';
    fprintf("%s\n\tCOSTO: v=%d\n",str,costo);
    fprintf("\tTRIPARTIZIONE:\n") 

    %display T
    fprintf("\t\tT={")
    Str=sprintf('(%d,%d) ',T');
    Str(end)=[];
    fprintf("%s",Str);
    fprintf("}\n");

    %display U
    fprintf("\t\tU={")
    Str=sprintf('(%d,%d) ',U');
    Str(end)=[];
    fprintf("%s",Str);
    fprintf("}\n");
    
    %display L
    fprintf("\t\tL={")
    Str=sprintf('(%d,%d) ',L');
    Str(end)=[];
    fprintf("%s",Str);
    fprintf("}\n");
    fprintf("\n\n\n\n\n\n\n\n\n\n\n\n\n")
    %vettore d
end