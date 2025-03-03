function [cost,condition,min_half_lambda,ang_locs_tagmode,hadamard_full,e] = fun_FindOptimalEncod(vesloc,modmat_phase,modmat_lookup,niters,cost_temp,offset)
condition_array = zeros(niters,1);
min_half_lambda_array = zeros(niters,1);

for i = 1:niters
    [cost,condition,min_half_lambda,ang_locs_tagmode,hadamard_full,e] = fun_CalEncodCost(vesloc,modmat_phase,modmat_lookup,offset);
    condition_array(i,1) = condition;
    min_half_lambda_array(i,1) = min_half_lambda;
    if cost < cost_temp
        cost_temp = cost;
        ang_locs_tagmode_temp = ang_locs_tagmode;
        condition_temp = condition;
        min_half_lambda_temp = min_half_lambda;

        hadamard_full_temp = hadamard_full;
        e_temp = e;
    end
end

cost = cost_temp;
condition = condition_temp;
min_half_lambda = min_half_lambda_temp;

ang_locs_tagmode = ang_locs_tagmode_temp;
[ang_locs_tagmode,ind] = unique(ang_locs_tagmode,'rows');

hadamard_full = hadamard_full_temp(ind,:);
e = e_temp;

end