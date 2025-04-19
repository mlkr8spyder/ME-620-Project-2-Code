function x_0 = gauss_seidel(A, b, x_0, iter, tol)
% Gauss-Seidel Solver
% A: coefficient Matrix
% B: right hand side vector
% x_0: initial solution guess
% iter: max number of iterations
% tol: tolerance of convergence

    for i = 1:iter
        x_old = x_0;
        for j = 1:size(A,1)
            x_0(j) = (b(j) - sum(A(j,:)'.*x_0) + A(j,j)*x_0(j)) / A(j,j);
        end
        if norm(x_0 - x_old, "inf") < tol
            fprintf('%d iterations \n', i)
            return
        end
    end
end

