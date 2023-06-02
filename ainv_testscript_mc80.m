function ainv_testscript_mc80(problem,droptol,droptype,shift,shiftA_flag,linsolve_tol)

linsolve_tol = str2num(linsolve_tol);
droptol = str2num(droptol);
shift = str2num(shift);
fprintf(strcat('\nProblem\t',problem,'\n'));
A = mmread(problem);
B = A - shift*speye(size(A));
if strcmpi(shiftA_flag,'true')
    A = B;
    fprintf('Coefficient matrix is A-s*I with s = %d. AINV also constructed for A-s*I.\n\n',shift);
else
    fprintf('Coefficient matrix is A. AINV constructed for A-s*I with s = %d.\n\n',shift);
end
    
% fprintf('Performing LDLT factorization of A.\n');
% [L,D,~,~] = ldl(A,0.5);
% fprintf('LDLT factorization done. Computing the inertia of A.\n');
% [np,nn,nz] = inertia_blkdiag(D);
% fprintf('The inertia of A is [%d %d %d]\n',np,nn,nz);
% nnzl = nnz(L);
% %%% nnzinvl = nnzl;  %nnz(invL);
% 
% fprintf('nnz(A) = %d, nnz(L) = %d\n\n',nnz(A),nnzl);
nnzl = nnz(A);
fprintf('Computing a favorable diagonally strong permutation.\n');
fprintf('Using hsl_mc80... \n');
n = length(A);

control.unmatched_last = false;
[order,scale,info] = hsl_mc80_order(B,'md',control);
disp(info);

S_prpc = spdiags(scale,0,n,n);
P_prpc = speye(size(A));
P_prpc = P_prpc(:,abs(order));

fprintf('Diagonally strong permutation computed.\n\n');
SB = P_prpc'*(S_prpc*B*S_prpc)*P_prpc;
SA = P_prpc'*(S_prpc*A*S_prpc)*P_prpc;

n = length(A);
s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);
b = randn(n,1);     b = b/norm(b);      %b = P_prpc'*(S_prpc*b);

fprintf('Solving the original linear system by unpreconditioned SQMR...\n');
[~,flag,relres,steps] = sqmr(A,b,linsolve_tol,1);
if flag == 0
    fprintf('Unpreconditioned SQMR converged in %d steps.\n',steps);
else
    fprintf('Unpreconditioned SQMR did not converge within %d steps ...\n',n);
    fprintf('Minimal residual %d reached at step %d.\n',relres,steps);
end
fprintf('Solving the MC64 preprocessed linear system by unpreconditioned SQMR...\n');
[~,flag,relres,steps] = sqmr(SA,P_prpc'*(S_prpc*b),linsolve_tol,1);
if flag == 0
    fprintf('Unpreconditioned SQMR converged in %d steps.\n\n',steps);
else
    fprintf('Unpreconditioned SQMR did not converge within %d steps ...\n',n);
    fprintf('Minimal residual %d reached at step %d.\n\n',relres,steps);
end

lindata.A = A;      lindata.b = b;      lindata.SB = SB;

lindata.P_prpc = P_prpc;	lindata.S_prpc = S_prpc;    lindata.nnzl = nnzl;
lindata.order = order;

mincost = realmax;
    
optdata = matrixtest_update_ainv_mc80(lindata,droptol,droptype,linsolve_tol);
if optdata.cost < mincost
    mincost = optdata.cost;
    minoptdata = optdata;
end
fprintf('\n');

fprintf('Optimal preconditioner droptol = %.2d\n\tdensity factors = [%.3d %.3d], iter = %d, cost = %.2d\n\n\n',...
    minoptdata.droptol,minoptdata.density(1),minoptdata.density(2),minoptdata.iter,minoptdata.cost);

end