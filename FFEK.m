function FEEK(b,ind_cu,num_nodi);

%ARGOMENTI
    
    %ind_cu: 5 componenti

    % 1: nodo di PARTENZA dell'arco
    % 2: nodo di ARRIVO dell'arco
    % 3: il loro tipo: 1=T, 2=U, 3=L   ->si puo mettere 0 (non servono)
    % 4: COSTO dell'arco               ->si puo mettere 0 (non servono)
    % 5: CAPACITA' dell'arco

    %num_nodi = numero di nodi della rete

    %b= vettore bilanci ai nodi (PER INTERO) ->si puo mettere 0 (non servono)
    %se sorgente e destinazione non sono questi modificare questi parametri
    S=1;
    T=6;

    ind_cu=sortrows(ind_cu,1:2);
    num_archi=(size(ind_cu,1));
    val=0;
    no_caum=0;

    lavoro=[ind_cu(:,[1,2,5]),(1:num_archi)']; %nodo 1, nodo 2,capacità indice in X
    lavoro1=lavoro;
    
    x=zeros(1,num_archi);

    R=zeros(num_archi);

    for i= 1:num_archi
        R(lavoro(i,1),lavoro(i,2))=lavoro(i,3);
    end

    while no_caum==0

        v=[S];
        dv=[]; %nodo da visitare,nodo di provenienza,visitato o no (1 da visitare 0 visitato)

        for i=1:num_archi
            if lavoro(i,1)==S
                dv=[dv;[lavoro(i,2),S,1]];
            end
        end

        i=1;
        while ~all(dv(:,3)==0)
            v=[v;dv(i,1)]; %sposto il nodo fra i visitati
            dv(i,3)=0;  %segno che l'ho visitato
            for k=1:num_archi
                if lavoro(k,1)==v(i+1) && lavoro(k,3)>0 && ~ismember(lavoro(k,2), dv(:,1)) %c'e un arco verso il nodo, ha ancora capacità a disposizione, e non ho ancora messo quel nodo in dv
                    dv=[dv;[lavoro(k,2),v(i+1),1]];
                end
            end
            i=i+1;
        end

        nodo =T;
        Caum=[];
        if ismember(T,dv(:,1))
            Caum=[T];
            while nodo ~= S
                nuova_riga = find(dv(:, 1) == nodo);
                nodo_prec = dv(nuova_riga, 2);
                nodo=nodo_prec;
                Caum =[nodo, Caum];
            end

            vect_u=[]; %vettore delle capacità
            for i=1:size(Caum,2)-1
                vect_u=[vect_u,R(Caum(i),Caum(i+1))];
            end
            delta=min(vect_u);
            val=val+delta;
            %aggiorno la capacità sia in R che in lavoro
            for i=1:size(Caum,2)-1
                R(Caum(i),Caum(i+1)) = R(Caum(i),Caum(i+1))-delta;
                R(Caum(i+1),Caum(i)) = R(Caum(i+1),Caum(i))+delta;
                riga = find(all(lavoro(:,1:2) == [Caum(i),Caum(i+1)], 2)); %trova la riga di lavoro da aggiornare
                lavoro(riga,3)=lavoro(riga,3)-delta;
            end
            %aggiorno X
            for i=1:size(Caum,2)-1
                riga = find(all(lavoro(:,1:2) == [Caum(i),Caum(i+1)], 2));
                x(riga)=x(riga)+delta;
            end
            fprintf("\tCaum= ")
            str=sprintf('%d-',Caum);
            str(end)=[];
            fprintf("%s\n",str);
            fprintf("\tdelta=%d, v=%d\n", delta,val);
            fprintf("\tX=[")
            str=sprintf('%d,',x);
            str(end)=']';
            fprintf("%s\n\n",str);
            %display(Caum);
            %display(vect_u);
            %display(x);
        else
            no_caum=1;

        end

    end
    
    fprintf("Fine percorsi su grafo normale, procedo con grafo residuo\n")
    parte2=[lavoro1(:,2), lavoro1(:,1), x'];
    lavoro=[lavoro(:,1:3);parte2];
    num_archi=num_archi*2;
    no_caum=0;
    
    
    caum=1;
    while no_caum==0

        v=[S];
        dv=[]; %nodo da visitare,nodo di provenienza,visitato o no (1 da visitare 0 visitato)

        for i=1:num_archi
            if lavoro(i,1)==S
                dv=[dv;[lavoro(i,2),S,1]];
            end
        end

        i=1;
        while ~all(dv(:,3)==0)
            v=[v;dv(i,1)]; %sposto il nodo fra i visitati
            dv(i,3)=0;  %segno che l'ho visitato
            for k=1:num_archi
                if lavoro(k,1)==v(i+1) && lavoro(k,3)>0 && ~ismember(lavoro(k,2), dv(:,1)) %c'e un arco verso il nodo, ha ancora capacità a disposizione, e non ho ancora messo quel nodo in dv
                    dv=[dv;[lavoro(k,2),v(i+1),1]];
                end
            end
            i=i+1;
        end

        nodo =T;
        Caum=[];
        if ismember(T,dv(:,1))
            Caum=[T];
            while nodo ~= S
                nuova_riga = find(dv(:, 1) == nodo);
                nodo_prec = dv(nuova_riga, 2);
                nodo=nodo_prec;
                Caum =[nodo, Caum];
            end

            vect_u=[]; %vettore delle capacità
            for i=1:size(Caum,2)-1
                vect_u=[vect_u,R(Caum(i),Caum(i+1))];
            end
            delta=min(vect_u);
            val=val+delta;
            %aggiorno la capacità sia in R che in lavoro
            for i=1:size(Caum,2)-1
                R(Caum(i),Caum(i+1)) = R(Caum(i),Caum(i+1))-delta;
                R(Caum(i+1),Caum(i)) = R(Caum(i+1),Caum(i))+delta;
                riga = find(all(lavoro(:,1:2) == [Caum(i),Caum(i+1)], 2));
                lavoro(riga,3)=lavoro(riga,3)-delta;
            end
            %aggiorno X
            for i=1:size(Caum,2)-1
                riga = find(all(lavoro(:,1:2) == [Caum(i),Caum(i+1)], 2));
                x(riga)=x(riga)+delta;
            end
            fprintf("\tCaum= ")
            str=sprintf('%d-',Caum);
            str(end)=[];
            fprintf("%s\n",str);
            fprintf("\tdelta=%d, v=%d\n", delta,val);
            fprintf("\tX=[")
            str=sprintf('%d,',x);
            str(end)=']';
            fprintf("%s\n\n",str);
            %display(Caum);
            %display(vect_u);
            %display(x);
        else
            no_caum=1;
            fprintf("finiti i Caum anche sul grafo residuo -> x ottima\n")
        end

    end
   


end
