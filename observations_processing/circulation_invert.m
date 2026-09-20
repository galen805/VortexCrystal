vort_mean = cell(numK,1);
ut_inv    = cell(numK,1);
psi_prof  = cell(numK,1);
rc_keep   = cell(numK,1);

for k = 1:numK

    % center
    i0 = centers(k,1);
    j0 = centers(k,2);
    x0 = Xg(i0,j0);
    y0 = Yg(i0,j0);

    % polar radius
    r = hypot(Xg - x0, Yg - y0);

    % ring mean vorticity
    vort_bar = NaN(size(rc));
    for j = 1:numel(rc)
        mask = rb(j) <= r & r < rb(j+1) & isfinite(vort_specfilt_calc);
        vort_bar(j) = mean(vort_specfilt_calc(mask),'omitnan');
    end

    % invert for u_theta via circulation
    % takes the place of symmetrization procedure!
    vort_use = vort_bar;
    vort_use(~isfinite(vort_use)) = 0;
    Gamma = cumtrapz(rc, vort_use .* 2*pi.*rc);
    ut_bar = Gamma ./ (2*pi.*max(rc, dr/2));

    % trim at first zero crossing of inverted velocity
    ut_trim   = ut_bar;
    vort_trim = vort_bar;
    r_trim    = rb(1:end-1); % or rc?

    s = sign(ut_trim);
    iref = find(s ~= 0 & isfinite(ut_trim), 1, 'first');

    if ~isempty(iref) && iref < numel(rc)
        s0 = s(iref);
        icross = find(sign(ut_trim(iref+1:end)) == -s0 | ut_trim(iref+1:end) == 0, 1, 'first');
        if ~isempty(icross)
            icut = iref + icross - 1;
            ut_trim(icut+1:end)   = NaN;
            vort_trim(icut+1:end) = NaN;
            r_trim(icut+1:end)    = NaN;
        end
    end

    keep = isfinite(ut_trim) & isfinite(vort_trim) & isfinite(r_trim);

    r_k    = r_trim(keep);
    zeta_k = vort_trim(keep);
    ut_k   = ut_trim(keep);

    % streamfunction from u_theta: u_theta = dpsi/dr
    psi_k = cumtrapz(r_k, ut_k);
    psi_k = psi_k - psi_k(1);

    vort_mean{k} = zeta_k*sign(zeta_k(1));
    ut_inv{k}    = ut_k*sign(ut_k(10));
    rc_keep{k}   = r_k;
end

