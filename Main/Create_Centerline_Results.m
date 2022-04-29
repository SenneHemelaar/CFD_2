function [results] = Create_Centerline_Results(results)

results.p_x05  = results.p(:,round(length(results.p)/2));
results.p_y05  = results.p(round(length(results.p)/2),:);

results.u_x05  = results.vec_u(:,round(length(results.vec_u)/2));
results.u_y05  = results.vec_u(round(length(results.vec_u)/2),:);

results.v_x05  = results.vec_v(:,round(length(results.vec_v)/2));
results.v_y05  = results.vec_v(round(length(results.vec_v)/2),:);

results.xi_x05 = results.xi(:,round(length(results.xi)/2));
results.xi_y05 = results.xi(round(length(results.xi)/2),:);

[~, neworder]  = sort(lower(fieldnames(results)));
results        = orderfields(results, neworder);
end