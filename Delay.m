function X = Delay(Tc, n, q)

X = zeros(n*q, 1);

%% My way    
    
% for(i = 1:n)
%     for(j = 1:q)
%         X((i-1)*q+j, k) = Tc(i, j);
%     end
% end
   
%% Mezic way

for(i = 1:q)
    for(j = 1:n)
        X((i-1)*n+j, 1) = Tc(j, i);
    end
end

return