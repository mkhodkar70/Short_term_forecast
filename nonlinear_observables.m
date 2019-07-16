function u_nl = nonlinear_observables(u_l, J)

n = size(u_l, 1);
m = size(u_l, 2);

u_nl = zeros(J*n, m);

k = 0;
for(i = 1:n)
   for(j = i:i+J)
      k = k + 1;
      if(j > n)
        j_p = j - n;
      else
        j_p = j;
      end
      u_nl(k, :) = u_l(i, :).*u_l(j_p, :);
   end
end

return

