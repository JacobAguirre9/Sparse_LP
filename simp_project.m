function x = simp_project(p)
% Computes: x = proj_{Gamma}(p) where Gamma = {x: e'*x<=1, x>=0}
[x,~,~] = simp_quadmin(-p,1);
end
