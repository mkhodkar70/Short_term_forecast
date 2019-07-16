function [X, Y] = Hankel_DMD_Matrices(Tc, n, m, q)

H = zeros(n*q, m-q+1);

%% My way    
    
% for(k = 1:m-q+1)
%     for(i = 1:n)
%         for(j = 1:q)
%             H((i-1)*q+j, k) = Tc(i, j+k-1);
%         end
%     end
% end
   
%% Mezic way

for(k = 1:m-q+1)
    for(i = 1:q)
        for(j = 1:n)
            H((i-1)*n+j, k) = Tc(j, i+k-1);
        end
    end
end


X = H(:, 1:end-1);
Y = H(:, 2:end);

return