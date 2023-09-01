function PHI = phi_calc(p, i, K, numericSTT)
PHI = 0;
if p == 1
    PHI = numericSTT{1}(i,K(1));
elseif p == 2
    PHI = numericSTT{2}(i,K(1),K(2));
end
end