%% Navigate to Model Functions
fpath = mfilename('fullpath');
search = strfind(fpath, '/');
if isempty(search)
    search = strfind(fpath, '\');
end
fpath = fpath(1:search(end));
model_path = strcat(fpath(1:search(end-1)), "Model_Script/");
%% Run Simulations
parpool(28)
num_samples = 10000;
iter = 10;
for i = 1:iter
    traj_rela = zeros(721, 2, num_samples);
    traj_crel = zeros(721, 2, num_samples);
    params = zeros(6, num_samples);
    cd(model_path)
    parfor n=1:num_samples
        [traj_rela(:, :, n) , traj_crel(:, :, n), params(:, n)] = simulate_param_vary()
    end
    cd(fpath)
    save(strcat('simulation_results_', num2str(i), '.mat'), 'traj_rela', 'traj_crel', 'params')
end
%% Parameter Variation Simulation Function
function [traj_rela, traj_crel, params] = simulate_param_vary()

params = rand(1, 6);

%IkB-NFkB Kd
%sample values between 3.2e-7 and 2.9e-3:
% IkBaRela dissociation
params(1) = exp((params(1)* (log(2.9e-3) - log(3.2e-7))) + log(3.2e-7)) * 200;
% IkBecRel dissociation
params(2) = exp((params(2)* (log(2.9e-3) - log(3.2e-7))) + log(3.2e-7)) * 55.93008739;
% IkBacRel dissociation
params(3) = exp((params(3)* (log(2.9e-3) - log(3.2e-7))) + log(3.2e-7)) * 125.093633;
% IkBeRela dissociation
params(4) = exp((params(4)* (log(2.9e-3) - log(3.2e-7))) + log(3.2e-7)) * 55.93008739;

%IKK-IkB Kd
%sample values between 0.1 and 1.6
% IKK+IkBa association
params(5) = 38/exp((params(5)* (log(1.6) - log(0.1))) + log(0.1)) ;
% IKK+IkBe association
params(6) = 38/exp((params(6)* (log(1.6) - log(0.1))) + log(0.1)) ;

p_mod = [
    19 1 params(1)
    20 1 params(1)
    41 1 params(2)         
    42 1 params(2)         
    58 1 params(3)         
    60 1 params(3)         
    59 1 params(4)         
    61 1 params(4)         
    21 1 params(5)         
    22 1 params(5)
    62 1 params(5)
    43 1 params(6)         
    44 1 params(6)
    63 1 params(6)];    

%% Run simulations with parameter modifications
names = {'IkBaRelA','IkBeRelA','IKKIkBaRelA','IKKIkBeRelA','RelA','IkBaRelAn','IkBeRelAn','RelAn','IkBacRel','IkBecRel','IKKIkBacRel','IKKIkBecRel','cRel','IkBacReln','IkBecReln','cReln','IKK'};
options = struct;
length_t = 12;
options.SIM_TIME = length_t*60;

traj_rela = [];
traj_crel = [];
%% TNF SIMULATION
dose_scale = 1/5200; % Dose scaling from Shannon et al (2007)
doses = 10;
output_tnf = [];
[~,x,simdata] = nfkbSimulate({'TNF',doses*dose_scale},names, p_mod, {},options);
output_tnf = cat(3,output_tnf,x);
[traj_rela(:,1), traj_crel(:,1)] = convert_to_ratio(output_tnf);
%% LPS SIMULATION
doses = 10;
dose_scale = 1/24000; % LPS molecular weight estimated between 10KDa and 100KDa
options.STEADY_STATE = simdata.STEADY_STATE;
output_lps = [];
[~,x] = nfkbSimulate({'LPS',doses*dose_scale},names, p_mod, {},options);
output_lps = cat(3,output_lps,x);
[traj_rela(:,2), traj_crel(:,2)] = convert_to_ratio(output_lps);   
end