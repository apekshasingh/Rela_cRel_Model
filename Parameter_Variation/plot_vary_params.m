%% Get current directory
fpath = mfilename('fullpath');
search = strfind(fpath, '/');
if isempty(search)
    search = strfind(fpath, '\');
end
fpath = fpath(1:search(end));
cd(fpath)

%% Load Simulation Results & Compute Trajectory Features

TNF_feats_r = [];
TNF_feats_c = [];
LPS_feats_r = [];
LPS_feats_c = [];
params_collect = [];
traj_rela_collect = [];
traj_crel_collect = [];

for i = 1:10
    load(strcat('simulation_results_', num2str(i), '.mat'))
    params_collect = cat(2, params_collect, params);
    traj_rela_collect = cat(3, traj_rela_collect, traj_rela);
    traj_crel_collect = cat(3, traj_crel_collect, traj_crel);
    [feats_r, feats_c] = get_traj_feats(squeeze(traj_rela(:, 1, :)),  squeeze(traj_crel(:, 1, :)));
    TNF_feats_r = cat(2, TNF_feats_r, feats_r);
    TNF_feats_c = cat(2, TNF_feats_c, feats_c);
    [feats_r, feats_c] = get_traj_feats(squeeze(traj_rela(:, 2, :)),  squeeze(traj_crel(:, 2, :)));
    LPS_feats_r = cat(2, LPS_feats_r, feats_r);
    LPS_feats_c = cat(2, LPS_feats_c, feats_c);
end
%% Navigate to model directory
model_path = strcat(fpath(1:search(end-1)), "Model_Script/");
cd(model_path)
%% Run baseline simulation (WT)
names = {'IkBaRelA','IkBeRelA','IKKIkBaRelA','IKKIkBeRelA','RelA','IkBaRelAn','IkBeRelAn','RelAn','IkBacRel','IkBecRel','IKKIkBacRel','IKKIkBecRel','cRel','IkBacReln','IkBecReln','cReln','IKK'};
options = struct;
length_t = 12;
options.SIM_TIME = length_t*60;
traj_rela = [];
traj_crel = [];
% TNF SIMULATION
dose_scale = 1/5200; % Dose scaling from Shannon et al (2007)
doses = 10;
output_tnf = [];
[~,x,simdata] = nfkbSimulate({'TNF',doses*dose_scale},names, [], {},options);
output_tnf = cat(3,output_tnf,x);
[traj_rela(:,1), traj_crel(:,1)] = convert_to_ratio(output_tnf);
% LPS SIMULATION
doses = 10;
dose_scale = 1/24000; % LPS molecular weight estimated between 10KDa and 100KDa
options.STEADY_STATE = simdata.STEADY_STATE;
output_lps = [];
[~,x] = nfkbSimulate({'LPS',doses*dose_scale},names, [], {},options);
output_lps = cat(3,output_lps,x);
[traj_rela(:,2), traj_crel(:,2)] = convert_to_ratio(output_lps);  
% Append to simulations
cd(fpath)
params = [8e-3; 0.00003576; 0.064; 0.08; 190; 46.1429];
params_collect = [params_collect, params];
traj_rela_collect = cat(3, traj_rela_collect, traj_rela);
traj_crel_collect = cat(3, traj_crel_collect, traj_crel);
[feats_r, feats_c] = get_traj_feats(traj_rela(:, 1),  traj_crel(:, 1));
TNF_feats_r = [TNF_feats_r, feats_r];
TNF_feats_c = [TNF_feats_c, feats_c];
[feats_r, feats_c] = get_traj_feats(traj_rela(:, 2),  traj_crel(:, 2));
LPS_feats_r = [LPS_feats_r, feats_r];
LPS_feats_c = [LPS_feats_c, feats_c];
%traj_feats = normalize([TNF_feats_r; TNF_feats_c; LPS_feats_r; LPS_feats_c]');
traj_feats = normalize([TNF_feats_r-TNF_feats_c; LPS_feats_r-LPS_feats_c]');
clear("params", "feats_r", "feats_c", "TNF_feats_r", "TNF_feats_c", "LPS_feats_r", "LPS_feats_c", "traj_rela", "traj_crel")
%% Convert Parameters to Kd
params_collect(1, :) = params_collect(1, :)/200;
params_collect(2, :) = params_collect(2, :)/55.93008739;
params_collect(3, :) = params_collect(3, :)/125.093633;
params_collect(4, :) = params_collect(4, :)/55.93008739;
params_collect(5, :) = 38./params_collect(5, :);
params_collect(6, :) = 38./params_collect(6, :);
param_names = ["IkBa-Rela Kd", "IkBe-cRel Kd", "IkBa-cRel Kd", "IkBe-Rela Kd", "IKK-IkBa Kd", "IKK-IkBe Kd"]; 
%% Define Clusters based on Features
rng(1)
[idx, C, sumd, D] = kmeans(traj_feats, 5);
rela_col = sscanf('155098','%2x%2x%2x',[1 3])/255;
crel_col = sscanf('C4005B','%2x%2x%2x',[1 3])/255;
sim_time = linspace(0,options.SIM_TIME/60,options.SIM_TIME+1);
%find representatives for each cluster & plot
clus_rep = zeros(1, 5);
for clus = 1:5
   [~, id] = min(D(:, clus));
   clus_rep(clus) = id;
   sumd(clus) = sumd(clus)/sum(idx==clus);
   figure()
   subplot(1, 2, 1)
   plot(sim_time, traj_rela_collect(:, 1, id), 'Color', rela_col, 'LineWidth', 2)
   hold on
   plot(sim_time, traj_crel_collect(:, 1, id), 'Color', crel_col, 'LineWidth', 2)
   ylim([0, 2])
   xlim([0, options.SIM_TIME/60])
   title(strcat("Cluster ", num2str(clus), ":TNF"))
   ylabel("Nuclear/Total")
   xlabel("Time (hr)")
   subplot(1, 2, 2)
   plot(sim_time, traj_rela_collect(:, 2, id), 'Color', rela_col, 'LineWidth', 2)
   hold on
   plot(sim_time, traj_crel_collect(:, 2, id), 'Color', crel_col, 'LineWidth', 2)
   text(1, 1.7, strcat(param_names', " = ", num2str(params_collect(:, id))))
   ylim([0, 2])
   xlim([0, options.SIM_TIME/60])
   title(strcat("Cluster ", num2str(clus), ":LPS"))
end
hold off

% plot baseline
figure()
subplot(1, 2, 1)
plot(sim_time, traj_rela_collect(:, 1, end), 'Color', rela_col, 'LineWidth', 2)
hold on
plot(sim_time, traj_crel_collect(:, 1, end), 'Color', crel_col, 'LineWidth', 2)
ylim([0, 2])
xlim([0, options.SIM_TIME/60])
title("Baseline: TNF")
ylabel("Nuclear/Total")
xlabel("Time (hr)")
subplot(1, 2, 2)
plot(sim_time, traj_rela_collect(:, 2, end), 'Color', rela_col, 'LineWidth', 2)
hold on
plot(sim_time, traj_crel_collect(:, 2, end), 'Color', crel_col, 'LineWidth', 2)
text(1, 1.7, strcat(param_names', " = ", num2str(params_collect(:, end))))
ylim([0, 2])
xlim([0, options.SIM_TIME/60])
title("Baseline: LPS")

figure()
histogram(categorical(idx), 'FaceColor', [.7 .7 .7])
title("Num Param Sets per Cluster")
figure()
bar(sumd, 'FaceColor', [.7 .7 .7])
title("Avg Dist from Cluster Mean")

figure()
[coef, score, ~, ~, explained] = pca(traj_feats);
for i= 1:5
    scatter(score(idx==i, 1), score(idx==i, 2), 1, "filled");
    hold on
end
scatter(score(end, 1), score(end, 2), 100, "k", 'Marker', "x");
xlabel(strcat("PC 1 (", num2str(explained(1)), "%)"))
ylabel(strcat("PC 2 (", num2str(explained(2)), "%)"))

figure()
for i = 1:6
    figure()
    scatter(score(:, 1), score(:, 2), 1, params_collect(i, :), 'filled')
    xlabel(strcat("PC 1 (", num2str(explained(1)), "%)"))
    ylabel(strcat("PC 2 (", num2str(explained(2)), "%)"))
    title(param_names(i))
end
