H = 0.642;
FA = 0.3;
k = sqrt(2)*norminv(FA/2);

syms dp
eqn = .642 == normcdf((k+dp)/sqrt(2)) + normcdf((k-dp)/sqrt(2))

solx = solve(eqn,dp)



H = hitRate_first_adj_all
FA = FARate_first_adj_some
k = sqrt(2).*norminv(FA/2)

syms dp
dprime = zeros(size(hitRate_first_adj_all))
for i = 1:size(hitRate_first_adj_all,1)
    for j = 1:size(hitRate_first_adj_all,2)
    H = hitRate_first_adj_all(i,j);
    k = sqrt(2).*norminv(FA(i,j)/2);
    eqn = H == normcdf((k+dp)/sqrt(2)) + normcdf((k-dp)/sqrt(2));
    
    dprime(i,j) = abs(solve(eqn,dp));
    
    end
end