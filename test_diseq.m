function test_diseq()
    % Initialize random seed for reproducibility
    rng(42);
    
    all_tests_passed = true;

    for i = 1:200
        % Generate random points
        Xa = randi([-100, 100]);
        Ya = randi([-100, 100]);
        Xb = randi([-100, 100]);
        Yb = randi([-100, 100]);
        
        % Ensure points are not the same
        while Xa == Xb && Ya == Yb
            Xa = randi([-100, 100]);
            Ya = randi([-100, 100]);
        end
        
        % Random value for c(5), can be 0, 1, or -1
        z = randi([-1, 1]);
        
        % Combine into input vector
        c = [Xa, Ya, Xb, Yb, z];
        
        % Call the function
        [a1, a2, b] = diseq(c);
        
        % Check the inequality
        if z == 0
            % (0,0) should satisfy the inequality
            test_result = (a1 * 0 + a2 * 0 <= b);
        elseif z == 1
            % (0,0) should not satisfy the inequality
            test_result = ~(a1 * 0 + a2 * 0 <= b) && (a1 * 0 + a2 * 1 <= b);
        elseif z == -1
            % Special case where (0,1) should be checked
            test_result = (a1 * 0 + a2 * 1 <= b) == (a1 == 0 && a2 > 0);
        end
        
        if ~test_result
            all_tests_passed = false;
            break;
        end
    end
    
    if all_tests_passed
        fprintf('OK\n');
    else
        fprintf('ERROR\n');
    end
end

function [a1, a2, b] = diseq(c)
    Xa = c(1);
    Ya = c(2);
    Xb = c(3);
    Yb = c(4);
    z = c(5);

    if size(c, 1) ~= 1 || size(c, 2) ~= 5
        fprintf("dimensione di c errata\n");
        return;
    end
    
    if z ~= 0 && z ~= 1 && z ~= -1
        fprintf("valore di c(5) errato\n");
        return;
    end

    % Coefficients in form aX + bY = c
    Coeff = [Yb - Ya, Xa - Xb, Xa * (Yb - Ya) - Ya * (Xb - Xa)];

    % Check if aX + bY <= c is the correct half-plane, otherwise invert signs
    if (z == 0 && 0 > Coeff(3)) || (z ~= 0 && 0 < Coeff(3))
        Coeff = -Coeff;
    end

    if (Coeff(3) == 0 && Coeff(1) == 0) && (Coeff(2) * z > 0)
        Coeff = -Coeff;
    elseif Coeff(2) == 0 && Coeff(1) == 0 && (Coeff(3) * z > 0)
        Coeff = -Coeff;
    elseif Coeff(3) == 0 && ((Coeff(2) > 0 && z == 1) || (Coeff(2) < 0 && z == 0)) % test with point (0,1)
        Coeff = -Coeff;
    end

    a1 = Coeff(1);
    a2 = Coeff(2);
    b = Coeff(3);

    fprintf("\t%+d X1  %+d X2 <= %+d\n", Coeff);
end

% Run the test

