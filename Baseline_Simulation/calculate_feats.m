function feats = calculate_feats(t, expt, pts)

if isnan(t)
    Rela_expt = expt(:, 1);
    cRel_expt = expt(:, 2);
else
    Rela_expt = interp1(t, expt(:, 1), pts);
    cRel_expt = interp1(t, expt(:, 2), pts);
end

[pk_r, pkt_r] = max(Rela_expt(pts<=4));
[pk_c, pkt_c] = max(cRel_expt(pts<=4));

pkt_r = pts(pkt_r);
pkt_c = pts(pkt_c);

ta_r = trapz(pts(~isnan(Rela_expt)), Rela_expt(~isnan(Rela_expt)));
ta_c = trapz(pts(~isnan(cRel_expt)), cRel_expt(~isnan(cRel_expt)));

hm_r = find((Rela_expt>=pk_r/2), 1);
hm_c = find((cRel_expt>=pk_c/2), 1);

hm_r = pts(hm_r);
hm_c = pts(hm_c);

feats = [pk_r, pk_c, pkt_r, pkt_c, ta_r, ta_c, hm_r, hm_c];
end