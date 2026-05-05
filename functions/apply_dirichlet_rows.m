function [Amod, bmod] = apply_dirichlet_rows(A, b, topOn, topVal, botOn, botVal)
% Strongly enforce Dirichlet BCs by replacing the top/bottom rows with
% identity rows and RHS equal to the boundary value.
Amod = A; bmod = b;
nz = size(A,1);

if topOn
    Amod(1,:) = 0; Amod(1,1) = 1;
    if isnan(topVal), topVal = 0; end
    bmod(1)   = topVal;
end

if botOn
    Amod(end,:) = 0; Amod(end,end) = 1;
    if isnan(botVal), botVal = 0; end
    bmod(end)  = botVal;
end
end
