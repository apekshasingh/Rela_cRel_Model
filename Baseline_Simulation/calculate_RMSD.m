function RMSD = calculate_RMSD(t, expt, pts, rr_mdl, cr_mdl, norm)

Rela_expt = interp1(t, expt(:, 1), pts);
cRel_expt = interp1(t, expt(:, 2), pts);

RMSD = zeros(1, 2);
scale = zeros(1, 2);

RMSD(1) = rms(Rela_expt-rr_mdl', 'omitnan');
scale(1) = rms(Rela_expt,'omitnan');

RMSD(2) = rms(cRel_expt-cr_mdl', 'omitnan');
scale(2) = rms(cRel_expt,'omitnan');

if norm
    RMSD = RMSD./scale;
end
end