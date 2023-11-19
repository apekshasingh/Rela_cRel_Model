function [feats_r, feats_c] = get_traj_feats(traj_rela, traj_crel)

%traj = timepoints by samples

% maximum amplitude (peak amp and time)
pk_r = max(traj_rela, [], 1);
pk_c = max(traj_crel, [], 1);

% total activity
ta_r = trapz(traj_rela);
ta_c = trapz(traj_crel);

% time up to half max
%[~, hm_r] = max((traj_rela>=(pk_r/2)), [], 1);
%[~, hm_c] = max((traj_crel>=(pk_c/2)), [], 1);

% time to peak
[~, hm_r] = max((traj_rela>=(pk_r/1)), [], 1);
[~, hm_c] = max((traj_crel>=(pk_c/1)), [], 1);

% number of peaks
npk_r = zeros(1, size(traj_rela, 2));
npk_c = zeros(1, size(traj_crel, 2));

for i = 1:size(traj_rela, 2)
    npk_r(i) = numel(findpeaks(traj_rela(:, i)));
    npk_c(i) = numel(findpeaks(traj_crel(:, i)));
end

feats_r = [pk_r; ta_r; log(1./hm_r); npk_r];
feats_c = [pk_c; ta_c; log(1./hm_c); npk_c];
end