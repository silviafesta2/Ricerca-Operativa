function ris = taglio(b, ind)
    B=[];
    switch ind
        case 1
            B=[3,4,7,8,9,10,15,16,17,18,19,20,21,22];
        case 2
            B=[5,6,11,12,13,14,23,24,25,26,27,28,29,30];
        case 3
            B=[7,8,15,16,17,18];
        case 4
            B=[9,10,19,20,21,22];
        case 5
            B=[11,12,23,24,25,26];
        case 6
            B=[13,14,27,28,29,30];
        case 7
            B=[15,16];
        case 8
            B=[17,18];
        case 9
            B=[19,20];
        case 10
            B=[21,22];
        case 11
            B=[23,24];
        case 12
            B=[25,26];
        case 13
            B=[27,28];
        case 14
            B=[29,30];
        otherwise
            B=[];
    end
    ris = b;
    ris(B)=0;
end
