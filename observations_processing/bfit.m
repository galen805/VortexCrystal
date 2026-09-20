%% calculate RMSE
rmse_map = nan(length(b_cpc_vec),length(b_pc_vec));

for ii = 1:length(b_cpc_vec)

    b_cpc_test = b_cpc_vec(ii);

    for jj = 1:length(b_pc_vec)

        b_pc_test = b_pc_vec(jj);

        ptest = [b_cpc_test,b_pc_test];

        rmse_map(ii,jj) = signed_velocity_rmse( ...
            ptest,xdat,ydat,u_obs,v_obs, ...
            centers_x0,centers_y0, ...
            Vm_cpc,Rm_cpc,Vm_pc,Rm_pc);
    end
end


%% Best-fit values

[min_rmse,idx] = min(rmse_map(:));
[i_best,j_best] = ind2sub(size(rmse_map),idx);

b_cpc_best = b_cpc_vec(i_best);
b_pc_best  = b_pc_vec(j_best);




%% Local function
function rmse = signed_velocity_rmse(p,x,y,u_obs,v_obs, ...
    xc,yc, ...
    Vm_cpc,Rm_cpc,Vm_pc,Rm_pc)

Vcyclone = @(Vm,Rm,r,b) Vm./Rm .* r .* ...
    exp((1./b) .* (1 - (r./Rm).^b));

u_calc = zeros(size(x));
v_calc = zeros(size(x));

% CPCs: centers 2:numK
for k = 2:length(xc)

    dx = x - xc(k);
    dy = y - yc(k);
    r  = hypot(dx,dy);
    rs = max(r,eps);

    V = Vcyclone(Vm_cpc,Rm_cpc,r,p(1));

    % Counterclockwise cyclonic velocity
    u_calc = u_calc - V .* dy ./ rs;
    v_calc = v_calc + V .* dx ./ rs;
end

% Polar cyclone: center 1
dx = x - xc(1);
dy = y - yc(1);
r  = hypot(dx,dy);
rs = max(r,eps);

V = Vcyclone(Vm_pc,Rm_pc,r,p(2));

u_calc = u_calc - V .* dy ./ rs;
v_calc = v_calc + V .* dx ./ rs;

% Vector velocity RMSE
rmse = sqrt(mean((u_calc-u_obs).^2 + (v_calc-v_obs).^2));
end