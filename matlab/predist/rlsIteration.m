function [h_out, P_out] = rlsIteration(phi_n, y_n, h, P, lambda)
    
    alpha = y_n-h.'*phi_n;

    % update gain vector
    Px = P*phi_n;
    k = Px/(lambda+phi_n'*Px);

    % matrix update
    P_out = (1/lambda)*(P-k*phi_n'*P);

    % update kernel vector (check if conj(alpha) needed)
    h_out = h + k*alpha;
    
end

