% file:        simplex.m
% description: Solves a linear program using the revised 
%              simplex algorithm
% input:       aN, matrix (A, in the canonical form)
%              b,  transposed array
%              cN, array (c, in the canonical form) 

function simplex (aN, b, cN)
    m = size(aN)(1);
    n = size(aN)(2);

    B = B_inv = eye(m);
    xN = 1:n;
    xB = n+1:n+m;
    cB = zeros(1, m);

    o = optimality(cB, B_inv, aN, cN);
    
    while o > 0
        aK = aN(:, o);
        mc = minimalCoeficient(B_inv, b, aK);

        % Switching columns
        [aN(:, o), B(:, mc)] = deal(B(:, mc), aN(:, o));

        % Switching cost variables
        [cN(o), cB(mc)] = deal(cB(mc), cN(o));

        % Switching variables
        [xN(o), xB(mc)] = deal(xB(mc), xN(o));

        B_inv = inv(B);

        o = optimality(cB, B_inv, aN, cN);
    end

    variables = xB;
    resp = (B_inv*b).';
    u = cB*B_inv*b;

    fprintf('Maximum value is %f.\n', u);

    for i = 1:m
        fprintf('x_%d = %f\n', xB(i), resp(i));
    end

    fprintf('Other variables are 0.\n');

endfunction

% Solves the optimality criteria. The returned index will be 
% matched in the matrix that holds the non-basic variables.
% If it's -1, the linear program is solved.
function o = optimality(cB, B_inv, aN, cN)
    criteria = cB*B_inv*aN - cN;
    [m, index] = min(criteria);
    if m >= 0
        o = -1;
    else
        o = index;
    end 
endfunction

% Solves the minimal coeficient criteira. The returned 
% index will be matched in the matrix that holds the basic 
% variables.
function mc = minimalCoeficient(B_inv, b, aK)
    criteria = (B_inv*b)./(B_inv*aK);
    m = min(criteria(criteria > 0));
    if m == zeros(1, 0)
        mc = -1;
    else 
        mc = find(criteria == m);   
    end
endfunction
