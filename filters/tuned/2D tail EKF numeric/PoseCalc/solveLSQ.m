function[x, norm_err,delta_x_plus_x] = solveLSQ(A,b,x0)
    x = (A'*A)^-1*A'*b;
    % show the sum of the squares of the pixel error 
    norm_err = sum(b.^2);
    delta_x_plus_x = x + x0;
end
