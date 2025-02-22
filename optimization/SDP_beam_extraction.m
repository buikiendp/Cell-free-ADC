function [F_star_comm, F_star_sensing] = SDP_beam_extraction(F_star, H_comm)
    % Rank-1 Approx.
    [U, M, N] = size(H_comm);
    S = size(F_star, 3);
    F_star_comm = zeros(U, M, N);
    H_comm_stacked = reshape(H_comm, U, []);
    F_comm_sum = 0;
    for u=1:U
        Q_u = F_star(:, :, u).';
        h_u = H_comm_stacked(u, :).';
        f_u = (h_u' * Q_u * h_u)^(-1/2) * Q_u * h_u;
        F_star_u = reshape(f_u, M, N);
        F_star_comm(u, :, :) = F_star_u;
        F_u_hat = f_u * f_u';
        F_comm_sum = F_comm_sum + F_u_hat;
    end
    
    F_star_sensing = zeros(max(S-U, 1), M, N);
    
    %fprintf( 'Matrix A: %d %d\n', M, N);
    
    if S > U
        F_star_sum = sum(F_star, 3).';
        F_sens_sum = F_star_sum-F_comm_sum;        
        
        [rows, cols] = size(F_sens_sum); % Lấy kích thước ma trận
        fprintf( 'Matrix A: %d %d\n', rows, cols);
        %{
        for i = 1:rows
            for j = 1:cols
                if isnan(real(F_sens_sum(i, j)))
                    F_sens_sum(i, j) = 0;
                end
                if isnan(imag(F_sens_sum(i, j)))
                    F_sens_sum(i, j) = complex(real(F_sens_sum(i, j)), 0);
                end
            end
        end
        %}
        %if any(isnan(F_sens_sum), 'all')
        %    fprintf('F_sens_sum contains NaN values3.');
            %[row, col] = find(isnan(F_sens_sum));
            %fprintf('NaN found at rows: %s\n', mat2str(row));
            %fprintf('NaN found at cols: %s\n', mat2str(col));
        %end
        %F_sens_sum
        F_sens_sum(isnan(real(F_sens_sum))) = 0;
        F_sens_sum(isnan(imag(F_sens_sum))) = 0;

        %if any(isnan(F_sens_sum), 'all')
        %    error('F_sens_sum contains NaN values2.');
        %end
        
        %if any(isinf(F_sens_sum), 'all')
        %    error('F_sens_sum contains Inf values.');
        %end
        
        [eigenvec, D] = eig(F_sens_sum);
        [lambda, order] = sort(diag(D), 'descend');
        for s=U+1:S
            F_star_sensing(s-U, :, :) = reshape(sqrt(lambda(s-U))*eigenvec(:, order(s-U)), 1, M, N);
        end
     end

end