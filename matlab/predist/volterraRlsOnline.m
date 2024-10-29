function [y, h, P] = volterraRlsOnline(x, M, K, f_sys_forward, sys_state, options)
    arguments
        x (:, 1) double
        M double {mustBePositive}
        K double {mustBePositive}
        f_sys_forward function_handle
        sys_state structure
        options.lambda double = 1
        options.f_yTransform function_handle = @(y) y
    end
    lambda = options.lambda;
    f_yTransform = options.f_yTransform;

    addpath(fullfile(pathstr, '../methods/pm'));

    N = length(x);
    y = zeros(N, 1);

    szKernel = sum(arrayfun(@(k) nchoosek(M+k-1, k), 1:K));

    % Parameter Estimation Setup
    h = [1; zeros(szKernel-1, 1)];
    P = eye(szKernel);

    % Companion Vector (phi) setup
    phi_n = zeros(szKernel, 1);
    y_slice_n = zeros(M, 1);
    m_table = makeTablePM(M, K);

    % main loop
    for ii=1:N
        x_n = x(ii);
        
        % perform forward simulation
        [y_n, sys_state]  = f_sys_forward(x_n, M, K, h, sys_state);

        % generate companion vector
        y_slice_n = [y_n; y_slice_n(1:end-1)];
        y_slice_transformed_n = f_yTransform(y_slice_n);
        phi_n = applyTableKernelPM(phi_n, y_slice_transformed_n, M, m_table);

        % update parameter estimates
        [h, P] = rlsIteration(phi_n, x_n, h, P, lambda);

        % store output
        y(ii) = y_n;
    end

end