function [h_i, P_i, h_q, P_q] = volterraRlsOfflineIq(x, y, M, K, options)
    arguments
        x (:, 1) double
        y (:, 1) double
        M double {mustBePositive}
        K double {mustBePositive}
        options.lambda double = 1
        options.f_yTransform function_handle = @(y) y
    end
    lambda = options.lambda;

    currentFile = mfilename('fullpath');
    [pathstr, ~, ~] = fileparts(currentFile);
    addpath(fullfile(pathstr, '../methods/direct'));

    N = length(x);
    if length(y) ~= N
        error('Expected x and y to be the same length. Instead length(x) = %d, length(y) = %d.', N, length(y));
    end

    M_iq = 2 * M;
    szKernel = sum(arrayfun(@(k) nchoosek(M_iq+k-1, k), 1:K));

    % Parameter Estimation Setup
    h_i = [1; zeros(szKernel-1, 1)];
    P_i = eye(szKernel);
    h_q = [zeros(M, 1); 1; zeros(szKernel-M-1, 1)];
    P_q = eye(szKernel);

    % Companion Vector (phi) setup
    phi_n = zeros(szKernel, 1);
    x_slice_n = zeros(M, 1);
%     m_table = makeTablePM(M_iq, K);
    figure
    % main loop
    for ii=1:N
        x_n = x(ii);
        y_n = y(ii);

        % generate companion vector
        x_slice_n = [x_n; x_slice_n(1:end-1)];
        x_slice_iq_n = [real(x_slice_n); imag(x_slice_n)];
        phi_n = voltVecGen(x_slice_iq_n, 1:K); %applyTableKernelPM(phi_n, y_slice_iq_n, M_iq, m_table);

        % update parameter estimates
        [h_i, P_i] = rlsIteration(phi_n, real(y_n), h_i, P_i, lambda);
        [h_q, P_q] = rlsIteration(phi_n, imag(y_n), h_q, P_q, lambda);

        if mod(ii-1, 5) == 0
            
            h_ind = 0:szKernel-1;
            stem(h_ind, h_i);
            hold on 
            stem(h_ind, h_q)
            hold off
            title(sprintf('ii=%d', ii))
            drawnow
        end
    end

end