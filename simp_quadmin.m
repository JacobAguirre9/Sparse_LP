function [x, mu, z] = simp_quadmin(q, gamma)
% minimize q'*x + (gamma/2)*norm(x)^2 s.t. e'*x<=1, x>=0
% mu and z are the KKT multipliers of the solution
[sortq,ind] = sort(q);
cumsumq = cumsum(sortq);
% In this first loop we check for a solution with e'*x=1.
foundj = 0;
n = length(q);
for j = 1 : n
    %  suppose sortq(j) <= -mu <= sortq(j+1) 
    % testz = [0,..,0,sortq(j+1)+mu,...,sortq(n)+mu]
    % testx = [-mu-sortq(1),...,-mu-sortq(j),0,...,0] / gamma
    %e'* testx = (-mu*j - cumsumq(j)) / gamma
    % Since that has to equal 1, we get mu = -(gamma+cumsumq(j))/j
    mu = -(gamma + cumsumq(j)) / j;
    if sortq(j) <= -mu && (j == n || sortq(j+1) >= -mu)
        foundj = j;
        break
    end
end
assert(foundj > 0)
if mu >= 0
    xf = (-mu * ones(foundj,1) - sortq(1:foundj)) / gamma;
    xf = xf / sum(xf);  % mathematically unnecessary but fixes some roundoff issues
    x = zeros(n,1);
    x(ind(1:foundj)) = xf;
    z = q + mu;
    z(ind(1:foundj)) = 0;
    return
end

x = max(zeros(n,1), -q / gamma);
mu = 0;
z = q + gamma * x;
end



%
function [x, mu, z] = simp_quadmin(q, gamma)
    % x = simp_quadmin(q, gamma)
    %
    % Solves the quadratic minimization problem
    %
    % minimize q'*x + (gamma/2)*norm(x)^2 s.t. e'*x<=1, x>=0
    % 
    % where q is a vector of length n and gamma is a positive scalar.
    %
    % This function is based on the algorithm in
    %
    % P. K. Andersen and L. Vandenberghe, "CVXOPT: A Python Package for Convex
    % Optimization," 2014.
    %
    % It is not intended to be efficient for large problems, but it is useful for
    % small problems and for testing purposes.
    %
    % The function returns the optimal x, the optimal value of mu, and the
    % optimal value of z.
    
    % The algorithm is based on the following observation. Suppose x is a
    % solution to the problem. Then there exist mu and z such that
    %
    % 1. z = q + gamma*x
    % 2. e'*x = 1
    % 3. z >= 0
    % 4. z'*x = 0
    % 5. mu >= 0
    % 6. mu*(e'*x - 1) = 0
    % 7. z >= 0
    % 8. z'*x = 0
    %
    % We can then find x by solving the following problem:
    %
    % minimize q'*x + (gamma/2)*norm(x)^2 s.t. e'*x<=1, x>=0
    %
    % If mu = 0, then we have a solution. If mu > 0, then we have a solution if
    % and only if there exists a j such that
    %
    % 1. sort(q)(j) <= -mu <= sort(q)(j+1)
    %
    % In this case, the solution is given by
    %
    % x = [-mu-sort(q)(1),...,-mu-sort(q)(j),0,...,0] / gamma
    %
    % If mu < 0, then the solution is given by
    %
    % x = max(0,-q/gamma)
    
    % This function is based on the algorithm in
    %
    % P. K. Andersen and L. Vandenberghe, "CVXOPT: A Python Package for Convex
    % Optimization," 2014.
    
    % The algorithm is based on the following
    % minimize q'*x + (gamma/2)*norm(x)^2 s.t. e'*x<=1, x>=0
    % mu and z are the KKT multipliers of the solution
    [sortq,ind] = sort(q);
    cumsumq = cumsum(sortq);
    % In this first loop we check for a solution with e'*x=1.
    foundj = 0;
    n = length(q);
    for j = 1 : n
        %  suppose sortq(j) <= -mu <= sortq(j+1) 
        % testz = [0,..,0,sortq(j+1)+mu,...,sortq(n)+mu]
        % testx = [-mu-sortq(1),...,-mu-sortq(j),0,...,0] / gamma
        %e'* testx = (-mu*j - cumsumq(j)) / gamma
        % Since that has to equal 1, we get mu = -(gamma+cumsumq(j))/j
        mu = -(gamma + cumsumq(j)) / j;
        if sortq(j) <= -mu && (j == n || sortq(j+1) >= -mu)
            foundj = j;
            break
        end
    end
    assert(foundj > 0)
    if mu >= 0
        xf = (-mu * ones(foundj,1) - sortq(1:foundj)) / gamma;
        xf = xf / sum(xf);  % mathematically unnecessary but fixes some roundoff issues
        x = zeros(n,1);
        x(ind(1:foundj)) = xf;
        z = q + mu;
        z(ind(1:foundj)) = 0;
        return
    end
    
    x = max(zeros(n,1), -q / gamma);
    mu = 0;
    z = q + gamma * x;
    end
    
    
