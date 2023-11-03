fpath = mfilename('fullpath');
search = strfind(fpath, '/');
if isempty(search)
    search = strfind(fpath, '\');
end
spath = fpath(1:search(end-1));

for i = 1:10
    params = readmatrix(strcat("params_", num2str(i), ".txt"));
    traj_rela = nan(721, 2, 10000);
    traj_crel = nan(721, 2, 10000);
    for j = 1:10
        traj_rela(:, 1, (((j-1)*1000)+1):(j*1000)) = readmatrix(strcat("traj_rela_TNF_", num2str(i), "_", num2str(j), ".txt"));
        traj_rela(:, 2, (((j-1)*1000)+1):(j*1000)) = readmatrix(strcat("traj_rela_LPS_", num2str(i), "_", num2str(j), ".txt"));
        traj_crel(:, 1, (((j-1)*1000)+1):(j*1000)) = readmatrix(strcat("traj_crel_TNF_", num2str(i), "_", num2str(j), ".txt"));
        traj_crel(:, 2, (((j-1)*1000)+1):(j*1000)) = readmatrix(strcat("traj_crel_LPS_", num2str(i), "_", num2str(j), ".txt"));
    end
    save(strcat(spath, "simulation_results_", num2str(i), ".mat"), "params", "traj_rela", "traj_crel")
end