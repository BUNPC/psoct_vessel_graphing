function [V,D] = eig2(x)
% Batched eigendecomposition of 2x2 matrices
%
% FORMAT [V,D] = eig2(X)
% X - [N, 2, 2] array - input matrices
% V - [N, 2, 2] array - columns are eigenvectors
% D - [N, 2]    array - eigenvalues
%
% WARNING: for stability I take the absolute value before the square root.
% This function will therefore give wrong results if you expect complex 
% eigenvalues (the returned values will never be complex).

% YB 2022/03/18

    x00 = x(:, 1, 1);
    x11 = x(:, 2, 2);
    x01 = x(:, 1, 2);
    x10 = x(:, 2, 1);
    
    % eigenvalues
    t = x00 + x11;
    d = x00 .* x11 - x01 .* x10;
    D1 = (t - sqrt(abs(t.*t - 4*d))) / 2;
    D2 = (t + sqrt(abs(t.*t - 4*d))) / 2;
    
    D = cat(2, D1, D2);
    
    % sort eigenvalues (not needed?)
%     mask = D1 < D2;
%     D(mask, 1) = D2(mask);
%     D(mask, 2) = D1(mask);
%     D1 = D(:, 1);
%     D2 = D(:, 2);
    
    if (nargout > 1)
        % eigenvectors
        v1 = cat(2, -x01, x00 - D1);
        mask = all(v1 == 0, 2);
        v1 = (1 - mask) .* v1 + mask .* cat(2, x11 - D1, -x10);
        v1 = v1 ./ sqrt(sum(v1.*v1, 2));

        v2 = cat(2, -x01, x00 - D2);
        mask = all(v2 == 0, 2);
        v2 = (1 - mask) .* v2 + mask .* cat(2, x11 - D2, -x10);
        v2 = v2 ./ sqrt(sum(v2.*v2, 2));

        V = cat(3, v1, v2);
    end
end


% test
% ----
%
% x = randn([2, 2]) + 10 * eye(2);
% x = reshape(x, [1, 2, 2]);
% 
% [V, D] = eig(reshape(x, [2, 2]))
% 
% [V, D] = eig2(x);
% V = reshape(V, [2, 2])
% D