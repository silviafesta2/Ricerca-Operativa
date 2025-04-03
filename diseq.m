function  [a1,a2,b]=diseq(c)
%c(1:2) punto 1 
%c(3:4) punto 2
%c(5) puo assumere solo valori 0, 1 ,-1 
    %=0 sifgnifica che (0,0) appartiene al semipiano, 1 = non appartiene,
    % se la retta passa per (0,0) e non è verticale o orizzontale 1-> il punto (0,1) deve appartenere al
    % semipiano, 0-> il punto (0,1) non deve appartenere al semipiano
    %nel caso di  rette parallele ai due assi usare +1 e -1 per indicare il sempiano positivo e
    %quello negativo

    Xa=c(1);
    Ya=c(2);
    Xb=c(3);
    Yb=c(4);
    z=c(5); 
    
    if size(c,1) ~= 1 || size(c,2)~=5
        printf("dimensione di c errata\n");
        return;
    end
    if c(5)~=0 && c(5)~=1 && c(5)~=-1
        printf("valore di c(5) errato")
    end

    %coeff in forma aX+bY=c
    Coeff=[Yb-Ya, Xa-Xb, Xa*(Yb-Ya)-Ya*(Xb-Xa)];

    %controllo se aX+bY<=c è il sempiano corretto altrimenti inverto tutti
    %i segni
    if (z==0 && 0>Coeff(3)) || (z~=0 && 0<Coeff(3))
        Coeff=-Coeff;
    end

    %assi
    if (Coeff(3)==0 && Coeff(1)==0) && (Coeff(2)*z > 0) %asse X2
        Coeff=-Coeff;
    elseif (Coeff(3)==0 && Coeff(2)==0)&& (Coeff(1)*z > 0) %asse X1
        Coeff=-Coeff;
    elseif Coeff(3)==0 && ( (Coeff(2)>0 && z==1) || (Coeff(2)<0 && z==0) ) %passa per l'origine ma non è orizzontale o verticale prova effettuata con il punto (0,1)
        Coeff=-Coeff;
    end
    
    %rette orizzontali e verticali
    if (Coeff(1)==0 && Coeff(2)*z > 0) || (Coeff(2)==0 && Coeff(1)*z > 0)
        Coeff=-Coeff;
    end

    a1=Coeff(1);
    a2=Coeff(2);
    b=Coeff(3);

    [~,d1] = numden(sym(abs(a1)));
    [~,d2] =numden(sym(abs(a2)));
    [~,d3] =numden(sym(abs(b)));

   
    mcm=lcm(lcm(d1,d2), d3);

    a1 = a1*mcm;
    a2 = a2*mcm;
    b = b*mcm;

    Coeff(1) = a1;
    Coeff(2) = a2;
    Coeff(3) = b;


    MCD=gcd(gcd(abs(a1),abs(a2)),abs(b));

    a1 = a1/MCD;
    a2 = a2/MCD;
    b = b/MCD;

    Coeff(1) = a1;
    Coeff(2) = a2;
    Coeff(3) = b;

    fprintf("\t%+d X1  %+d X2 <= %+d\n",Coeff);
    %gestire il caso di coeff non interi
end

