function [M,D] = ainv_mc80_right(A,order,droptol,droptol_type)

% basic right looking AINV with absolute drop tolerance

n = length(A);
M = speye(n,n);
D = sparse(n,n);

sorder = zeros(n, 1);

for i = 1:n
    if order(i) > 0, sorder(i) = 1;
    else, sorder(i) = 2;
    end
end

i = 1;
while i <= n
    s = sorder(i);
    i_pvt_idx = i:i - 1 + s;

    j = i;
    while j <= n
        j_pvt_idx = j:j - 1 + s;
        D(j_pvt_idx, j_pvt_idx) = (A(i_pvt_idx,:) * M(:,j_pvt_idx));
        if nnz(D(j_pvt_idx, j_pvt_idx)) ~= 0 && j ~= i
            exact = M(:, j_pvt_idx) - M(:, i_pvt_idx) * inv(D(i_pvt_idx, i_pvt_idx)) * D(j_pvt_idx, j_pvt_idx);

            for jj = 1:s
                % remove elements by absolute drop tolerance
                [nzidx,~,nzv] = find(exact(:,jj));
                M(nzidx, jj+j-1) = nzv.*(abs(nzv) >= droptol);
            end
        end
        j = j + s;
    end
    i = i + s;
end

D
M