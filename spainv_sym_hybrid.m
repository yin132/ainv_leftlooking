function [M,D,p] = spainv_sym_hybrid(A,order,droptol,droptol_type)
%
% This subroutine gives an approximate inverse factorization of a real
% symmetric matrix such that M'*A*M = D, in left-looking style. The input 
% matrix A has been preprocessed by HSL's MC80, and we are not doing any 
% pivoting in this subrouine. Instead, we simply follow the partition of A 
% generated by mc80, passed in here through the input `order'. 
% order(dd) > 0 means the dd-th row/column of A is partitioned as a single 
% row/column (1x1 pivot should be used), whereas order(dd) < 0 means the 
% dd-th and (dd+1)-st rows/colums are partitioned together (2x2 pivot)
%


n = length(A);

normA = norm(A,1);
normAFro = norm(A,'fro');
% initial unit upper triangular M and block diagonal D (stored as a tri-
% diagonal and hence can be constructed explicitly from the main diagonal
% and subdiagonal, due to the real symmetry of A)
M = speye(n,n);
D = sparse(n,n);
invD = sparse(n,n);

% alpha = (1+sqrt(17))/8;
alpha = 1;
p = 1:n;
accept_small_pvt = false;
c = 1;

droptol_max = 5e-1;

pvttypesum = zeros(3,1);
small_sv_D = realmax;
large_sv_D = realmin;
progress = zeros(10,1);

dd = 1;
while dd <= n
    if order(dd) > 0,   s = 1;
    else,               s = 2;
    end
    pvt_idx = dd:dd-1+s;
    prj_pvt_cols = ((A*M(:,pvt_idx))'*M(:,c:dd-1))';
    
    M_pvtcols = M(:,pvt_idx) - M(:,c:dd-1)*(invD(c:dd-1,c:dd-1)*prj_pvt_cols);

    % dropping small elements in the pivot columns, one column at a time
    for jj = 1 : s
        % assemble all nonzero elements in the current pivot column
        [nzidx,~,tmpv] = find(M_pvtcols(1:dd-2+jj,jj));
        if ~isempty(nzidx)
            if strcmpi(droptol_type,'relative')
                % sort the nodnzero elements in the current pivot column,
                % and drop the smallest ones whose sum is no greater than
                % droptol times the total sum of all nonzero elements;
                % (note that we could sum up the elements themselves or
                % of their squares)
                [ordtmpv,ordix] = sort(conj(tmpv).*tmpv);   %sort(abs(tmpv));  %
                cumsum_ordtmpv = cumsum(ordtmpv);
                dropidx = min([find(cumsum_ordtmpv >= cumsum_ordtmpv(end)*droptol,1) ...
                    find(ordtmpv >= droptol_max,1)])-1;
                if dropidx > 0
                    tmpv(ordix(1:dropidx)) = 0;
                end
                %M(p(nzidx),p(dd-1+jj)) = tmpv;
                M_pvtcols(nzidx,jj) = tmpv;
            else
                % absolute drop tolerance means to simply drop all nonzero
                % elements smaller than droptol
                filter = min([droptol droptol_max]);
                %M(p(nzidx),p(dd-1+jj)) = tmpv.*(abs(tmpv) >= filter);
                M_pvtcols(nzidx,jj) = tmpv.*(abs(tmpv) >= filter);
            end
        end
    end
    
    nxt_pvt = M_pvtcols'*(A*M_pvtcols);
    
    %check is next entries to D are acceptable, if not, then switch to
    %right looking to look for new pivot
    if s == 1
        req_pivoting = (abs(nxt_pvt) <= 1e-6);
    else
        req_pivoting = norm(nxt_pvt(:,1),"inf") <= 10e-4 || norm(nxt_pvt(:,2),"inf") <= 10e-4;
    end
    
    if req_pivoting && (~accept_small_pvt)
        accept_small_pvt = false;
        %A-orthogonalize next 10? columns against the already processed columns
%         search_cols = min(10, n - dd);
        search_cols = n - dd;
        k = dd;
        temp = sparse(n,search_cols + 1);
        while k <= dd + search_cols && k <= n
            if order(k) > 0,    ks = 1;
            else,               ks = 2;
            end
            k_pvt_idx = k:k-1+ks;
            k_prj_pvt_cols = ((A*M(:,k_pvt_idx))'*M(:,c:dd-1))';
            k_M_pvtcols = M(:,k_pvt_idx) - M(:,c:dd-1)*(invD(c:dd-1,c:dd-1)*k_prj_pvt_cols);

            for jj = 1:ks
                % remove elements by absolute drop tolerance
                [nzidx,~,nzv] = find(k_M_pvtcols(:,jj));
                temp(nzidx, k-dd+jj) = nzv.*(abs(nzv) >= droptol);
%                 M(nzidx, k+jj) = nzv.*(abs(nzv) >= droptol);
            end

            m = k - dd + 1;
            k = k + ks;
        end
        c = dd;
        
        %pick out a pivot and perform swaps on A and M and order
        A1 = ((A*temp(:,1))'*temp)';
        A1_skip1 = A1;  A1_skip1(1) = 0;
        % use the threshold beta to find the second candidate pivot row/column
        if dd < n
            omega_1 = norm(A1_skip1(1:search_cols+1),'inf');
            r = find(abs(A1_skip1) >= omega_1,1);
        else
            omega_1 = 0;
        end

        % the leading 1x1 pivot is greater than or equal to alpha times the
        % selected off-diagonal element corresponding to the second candidate
        % coordinate pivot row/column (which itself is greather than or equal 
        % to beta times the largest off-diagonal element)
        if abs(A1(1)) >= alpha*omega_1
            % original is good enough, do nothing
            pvttypesum(1) = pvttypesum(1)+1;
            accept_small_pvt = true;
        % otherwise, the second candidate coordinate pivot row/column must be 
        % formed for further comparison and pivot selection
        else
            Ar = ((A*temp(:,r))'*temp)';
            Ar_skipr = Ar;      Ar_skipr(r) = 0;
            omega_r = norm(Ar_skipr,'inf');
            % the leading 1x1 pivot (i.e., (1,1) element) is relatively large
            if abs(A1(1))*omega_r >= alpha*omega_1^2
                % do nothing
                accept_small_pvt = true;
                pvttypesum(1) = pvttypesum(1)+1;
            % the second 1x1 pivot (i.e., (r,r) element) is relatively large
            elseif abs(Ar(r)) >= alpha*omega_r
                tmp = p(dd);     p(dd) = p(dd+r-1);    p(dd+r-1) = tmp;
                order(dd+r-1) = order(dd); order(dd) = 1;
                tmp_col = A(:,dd);
                A(:,dd) = A(:,dd+r-1); 
                A(:,dd+r-1) = tmp_col;
                tmp_rw = A(dd,:);
                A(dd,:) = A(dd+r-1,:); 
                A(dd+r-1,:) = tmp_rw;
                D(dd,dd) = Ar(r);
                pvttypesum(2) = pvttypesum(2)+1;
            % both the leading and the second 1x1 pivots are relatively small
            % compared to the (1,r) (r,1) elements; use 2x2 pivot instead
            else
                tmp = p(dd+1);   p(dd+1) = p(dd+r-1);  p(dd+r-1) = tmp;
                tmp_col = A(:,dd+1); A(:,dd+1) = A(:,dd+r-1);  A(:,dd+r-1) = tmp_col;
                tmp_rw = A(dd+1,:);
                A(dd+1,:) = A(dd+r-1,:); 
                A(dd+r-1,:) = tmp_rw;
                order(dd+r-1) = order(dd+1); order(dd) = -1; order(dd+1) = -1;
                D(dd,dd) = A1(1);   D(dd+1, dd+1) = Ar(r);
                D(dd,dd+1) = A1(r); D(dd+1, dd) = A1(r);
                pvttype = 3;    pvttypesum(3) = pvttypesum(3)+1;
            end
        end
        continue;
    end

    %return to left looking

    M(:,pvt_idx) = M_pvtcols;
    D(pvt_idx,pvt_idx) = nxt_pvt;

    if s == 1 && abs(D(dd,dd)) < 64*sqrt(eps)*normA%) || (s == 2 && cond(full(D(dd:dd+1,dd:dd+1)))>1e16)
        D(dd,dd) = sign(0.5+sign(real(D(dd,dd))))*1e-2*normA;
    elseif s == 2 && cond(full(D(dd:dd+1,dd:dd+1))) > 1/eps
        D(dd+1,dd+1) = -sign(0.5+sign(real(D(dd,dd))))*1e-2*normA;
        %disp('Step %d, singular block D appears.\n',dd);
    end
    invD(pvt_idx,pvt_idx) = D(pvt_idx,pvt_idx)\speye(s,s);
    if nnz(isnan(D(dd:dd-1+s,dd:dd-1+s)))>0 || nnz(isinf(D(dd:dd-1+s,dd:dd-1+s))) > 0
        fprintf('What is wrong with D?\n');
    end
    svdPvt = svd(full(D(dd:dd-1+s,dd:dd-1+s)));
    if min(svdPvt) < small_sv_D,    small_sv_D = min(svdPvt);   end
    if max(svdPvt) > large_sv_D,    large_sv_D = max(svdPvt);   end
    condD = large_sv_D/small_sv_D;

    dd = dd + s;
%     if nnz(isnan(M)) > 0
%         fprintf('What is wrong with M?\n');
%     end
    tens_percent_done = floor(dd/n*10);
    if tens_percent_done > 0 && tens_percent_done < 10 && progress(tens_percent_done) == 0
        progress(tens_percent_done) = 1;
        nnzM = nnz(M);
        %D = spdiags([DSubd DMain [0; conj(DSubd(1:end-1))]],-1:1,n,n);
        %R = M(:,p(1:dd-1))'*A*M(:,p(1:dd-1))-D(1:(dd-1),1:(dd-1));
        R = M(:,1:dd-1)'*A*M(:,1:dd-1)-D(1:dd-1,1:dd-1);
        fprintf(' %d0%%: nnz(M) = %d (%.3f x nnz(A)), r-res = %.2d, cond(Dk) = %.2d, pvt type = [%d %d %d], norm(M) = %.2d.\n',...
            tens_percent_done,nnzM,nnzM/nnz(A),norm(R,'fro')/normAFro, condD,[pvttypesum(1) pvttypesum(2) pvttypesum(3)],norm(M,1));
    end
%     if tens_percent_done >= 5
%         fprintf('Step %d, nnz(M) = %d, norm(M) = %.3d, cond(D) = %.3d.\n',dd,nnz(M),norm(M,'fro'),condD);
%     end
end
nnzM = nnz(M);
% if nnz(sort(p)-(1:n)) > 0
%     error('Wrong permutation vector p\n');
% end
%D = spdiags([DSubd DMain [0; conj(DSubd(1:end-1))]],-1:1,n,n);

M = M(p,:);

R = M'*A*M-D;
fprintf('100%%: nnz(M) = %d (%.3f x nnz(A)), r-res = %.2d, cond(Dk) = %.2d, pvt type = [%d %d %d], norm(M) = %.2d\n',...
    nnzM,nnzM/nnz(A),norm(R,'fro')/normAFro,condD,[pvttypesum(1) pvttypesum(2) pvttypesum(3)],norm(M,1));
fprintf('Drop tolerance = %d\n',droptol);
% M = M(p,p); % Unit upper triangular

D;
end