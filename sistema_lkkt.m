function [ris] = sistema_lkkt(grad_f,A_b,p) %deve restituire il vettore dei lamda
    p=p';
    A=A_b(:,1:2);
    b=A_b(:,3);
    num_moltiplicatori = size(A,1);
    %vettore dei moltiplicatori
    lamda=zeros(1,size(A,1));
    lamda_bin=zeros(1,size(A,1));
    
    %gradiente di f VALUTATO IN P
    grad_f_p = grad_f(:,1:2)*p+grad_f(:,3);

    %ogni COLONNA contiene il gradiente di una g(x) (VALUTATO O NO è UGUALE
    %TANTO è COSTANTE)
    grad_g=(A');

    
    %vettore contente le valutazioni di g in p (ogni riga ne contiene una)
    g_p=(A*p)-b;
    
    lamda_nulli=[];
    lamda_non_nulli=[];

    for j=1:num_moltiplicatori
        if g_p(j)~=0
            lamda_bin(j)=1;
            lamda_nulli=[lamda_nulli,j];
        end
    end   
    lamda_non_nulli=setdiff(1:num_moltiplicatori,lamda_nulli);
    matrice_sistema=grad_g(:, (lamda_bin == 0)); %prendo solo le colonne relative a lamda non nulli
    
    %gestire il caso in cui il sistema non ha soluzioni (GESTITO)
    fprintf("SISTEMA LKKT:\n"); %se i lamda non nulli sono più di due può creare problemi e succede
    
    fprintf('\tλ%d=0 \n',lamda_nulli);
    for j=1:2
        fprintf('\t')
        for k=1:size(lamda_non_nulli,2)
            fprintf('%s%s λ%d  ', c_segno(sym(matrice_sistema(j,k))), abs(sym(matrice_sistema(j,k))), lamda_non_nulli(k));
        end
        fprintf('= %s%s\n', c_segno(sym(-grad_f_p(j))) ,abs(sym(-grad_f_p(j))) );
    end

    

    if rank(matrice_sistema)~=min(size(matrice_sistema))  %sistema senza soluzioni
        fprintf("\n PUNTO NON STAZIONARIO (non esistono sol di LLKT con x=p)\n")
    else    
        sol=matrice_sistema\(-grad_f_p);
        lamda(lamda_non_nulli)=sol;
        display(sym(lamda), "moltiplicatori (lamda)");
    end
    ris=sym(lamda);
end