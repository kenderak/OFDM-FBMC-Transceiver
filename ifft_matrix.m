%% IFFT algorithm
function Xk = ifft_matrix(xn,W)
    Xk = fixp(W*xn);
end

