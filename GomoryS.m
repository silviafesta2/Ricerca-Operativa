function GomoryS(c,A,b)
    Aeq=[A,eye(size(A,1))];
    display(Aeq);
    f=[-c',zeros(1,(size(A,1)))];
    display(f);
    beq=b;
    display(beq);
    Xrc=linprog(f,[],[],Aeq,beq,zeros(size(Aeq,2),1));
    display(sym(Xrc));
    for 
end