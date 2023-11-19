%% Get current directory
fpath = mfilename('fullpath');
search = strfind(fpath, '/');
if isempty(search)
    search = strfind(fpath, '\');
end
fpath = fpath(1:search(end));
cd(fpath)
% Specify model directory
model_path = strcat(fpath(1:search(end-1)), "Model_Script/");
%% Specify Stimulations
run_TNF = true;
run_lps = true;
run_cpg= false;
run_pam= false;
run_poly=false;
run_fla=false;
run_r84=false;
%% Specify Parameter Modification
alpha_mut = false;
epsilon_KO = false;
modify_IkBeRelA = false;
modify_IKKIkBe = true;
%% Specify Plotting Options
plot_expt = false;
plot_sim = true;
%% Plotting experimental data and model simulations
length_t = 12;

untran_r = 0.130;
untran_c = 0.035;

names = {'IkBaRelA','IkBeRelA','IKKIkBaRelA','IKKIkBeRelA','RelA','IkBaRelAn','IkBeRelAn','RelAn','IkBacRel','IkBecRel','IKKIkBacRel','IKKIkBecRel','cRel','IkBacReln','IkBecReln','cReln'};
options = struct;
options.DEBUG = 1;
options.SIM_TIME = length_t*60;
time = linspace(0,options.SIM_TIME/60,options.SIM_TIME+1);

%set colors
rela_color = ['#1E87BB';'#31A6C3';'#4DC6C7';'#68D6D0';'#7BE0D6'];
crel_color = ['#CC204F';'#D34144';'#DB6138';'#E2822D';'#EAA221'];
for i = 1:5
    r_color(i,:) = [sscanf(rela_color(i,2:end),'%2x%2x%2x',[1 3])/255,0.6];
    c_color(i,:) = [sscanf(crel_color(i,2:end),'%2x%2x%2x',[1 3])/255,0.6];
end

%% TNF SIMULATION
cd(fpath)
if run_TNF
    
    if alpha_mut %reduced IkBa induced synthesis by Rela
        p_mod = [
            6 1 1.5e-5];
        expt = readmatrix('alphamut_TNF_traj.csv');
        expt = expt(1:2:end, :);
    elseif epsilon_KO %set translation to zero (hence no IkBe in system)
        p_mod = [
            30 1 0];
        expt = readmatrix('epsilonKO_TNF_traj.csv');
        expt = expt(1:2:end, :);
    elseif modify_IkBeRelA %original dissociation 0.08 (KD = 0.0014 to KD = 6.3937e-07)
        p_mod = [
            59 1 0.00003576
            61 1 0.00003576];
    elseif modify_IKKIkBe %original dissociation 38 (KD = 0.8235 to KD = 0.20)
        p_mod = [
            45 1 9.228
            46 1 9.228
            65 1 9.228];
    else
        p_mod = [];
        expt = readmatrix('WT_TNF_traj.csv');
        expt = expt(1:2:end, :);
    end
    
    dose_scale = 1/5200; % Dose scaling from Shannon et al (2007)
    doses = [10];
    
    cd(model_path)
    output_tnf = [];
    for i = 1:length(doses)
        if isempty(output_tnf)
            [t,x,simdata] = nfkbSimulate({'TNF',doses(i)*dose_scale},names, p_mod, {},options);
        else
            options.STEADY_STATE = simdata.STEADY_STATE;
            [~,x] = nfkbSimulate({'TNF',doses(i)*dose_scale},names, p_mod, {},options);
        end
        output_tnf = cat(3,output_tnf,x);
    end
    cd(fpath)

    figure
    hold on
    %Real Data
    if plot_expt
        for cell = 1:(size(expt, 2)/3)
            plot(expt(:,(3*(cell-1))+1),expt(:,(3*(cell-1))+2),'Color', r_color((cell-(5*floor((cell-1)/5))),:),'LineWidth',3)
            plot(expt(:,(3*(cell-1))+1),expt(:,(3*(cell-1))+3),'Color', c_color((cell-(5*floor((cell-1)/5))),:),'LineWidth',3)
        end
    end

    %RelA Calculations
    rela = [];
    rela_ratio = [];
    rela_cyto = sum(output_tnf(:,1:5),2);
    rela_nuc = sum(output_tnf(:,6:8),2);
    rela(:,1) = rela_nuc;
    rela_total = (3.5*rela_cyto+rela_nuc)/4.5;
    rela_ratio(:,1) = rela_nuc/(rela_total(1)+untran_r);
    rela_ratio(:,1) = rela_ratio(:,1) - rela_ratio(1,1);

    %cRel Calculations
    crel = [];
    crel_ratio = [];
    crel_cyto = sum(output_tnf(:,9:13),2);
    crel_nuc = sum(output_tnf(:,14:16),2);
    crel(:,1) = crel_nuc;
    crel_total = (3.5*crel_cyto+crel_nuc)/4.5;
    crel_ratio(:,1) = crel_nuc/(crel_total(1)+untran_c);
    crel_ratio(:,1) = crel_ratio(:,1) - crel_ratio(1,1);

    %legend('RelA ratio (sim)','cRel ratio (sim)', 'RelA ratio', 'cRel ratio')
    title('TNF Simulation')
    xlabel('Time(hr)')
    ylabel('Nuclear/Total Cell Intensity')
    axis([0 options.SIM_TIME/60 0 2])
    
    if plot_sim
        plot(time, rela_ratio(:,1),'Color',sscanf('155098','%2x%2x%2x',[1 3])/255, 'LineWidth',4)
        plot(time, crel_ratio(:,1),'Color',sscanf('C4005B','%2x%2x%2x',[1 3])/255,'LineWidth',4)
    end

    hold off

end
%% LPS Simulation
cd(fpath)
if run_lps
    
    if alpha_mut %reduced IkBa induced synthesis by Rela
        p_mod = [
            6 1 1.5e-5];
        expt = readmatrix('alphamut_LPS_traj.csv');
        expt = expt(1:2:end, :);
    elseif epsilon_KO %set translation to zero (hence no IkBe in system)
        p_mod = [
            30 1 0];
        expt = readmatrix('epsilonKO_LPS_traj.csv');
        expt = expt(1:2:end, :);
    elseif modify_IkBeRelA
        p_mod = [
            59 1 0.00003576
            61 1 0.00003576];
    elseif modify_IKKIkBe
        p_mod = [
            45 1 9.228
            46 1 9.228
            65 1 9.228];
    else
        p_mod = [];
        expt = readmatrix('WT_LPS_traj.csv');
        expt = expt(1:2:end, :);
    end

    %LPS Simulation
    doses = [10];
    dose_scale = 1/24000; % LPS molecular weight estimated between 10KDa and 100KDa
    
    cd(model_path)
    output_lps = [];
    for i = 1:length(doses)
        if isempty(output_lps)
            [t,x,simdata] = nfkbSimulate({'LPS',doses(i)*dose_scale},names, p_mod, {},options);
        else
            options.STEADY_STATE = simdata.STEADY_STATE;
            [~,x] = nfkbSimulate({'LPS',doses(i)*dose_scale},names, p_mod, {},options);
        end
        output_lps = cat(3,output_lps,x);
    end
    cd(fpath)

    %RelA Calculations
    rela_cyto = sum(output_lps(:,1:5),2);
    rela_nuc = sum(output_lps(:,6:8),2);
    rela(:,2) = rela_nuc;
    rela_total = (3.5*rela_cyto+rela_nuc)/4.5;
    rela_ratio(:,2) = rela_nuc/(rela_total(1)+untran_r);
    rela_ratio(:,2) = rela_ratio(:,2) - rela_ratio(1,2);

    %cRel Calculations
    crel_cyto = sum(output_lps(:,9:13),2);
    crel_nuc = sum(output_lps(:,14:16),2);
    crel(:,2) = crel_nuc;
    crel_total = (3.5*crel_cyto+crel_nuc)/4.5;
    crel_ratio(:,2) = crel_nuc/(crel_total(1)+untran_c);
    crel_ratio(:,2) = crel_ratio(:,2) - crel_ratio(1,2);

    figure
    hold on
    %Real Data
    if plot_expt
        for cell = 1:(size(expt, 2)/3)
            plot(expt(:,(3*(cell-1))+1),expt(:,(3*(cell-1))+2),'Color', r_color((cell-(5*floor((cell-1)/5))),:),'LineWidth',3)
            plot(expt(:,(3*(cell-1))+1),expt(:,(3*(cell-1))+3),'Color', c_color((cell-(5*floor((cell-1)/5))),:),'LineWidth',3)
        end
    end
    
    title('LPS Simulation')
    xlabel('Time(hr)')
    ylabel('Nuclear/Total Cell Intensity')
    axis([0 options.SIM_TIME/60 0 2])

    if plot_sim
        plot(time, rela_ratio(:,2),'Color',sscanf('155098','%2x%2x%2x',[1 3])/255, 'LineWidth',4)
        plot(time, crel_ratio(:,2),'Color',sscanf('C4005B','%2x%2x%2x',[1 3])/255,'LineWidth',4)
    end
    
    hold off

end
%% CpG Simulation
cd(fpath)
if run_cpg

    if alpha_mut %reduced IkBa induced synthesis by Rela
        p_mod = [
            6 1 1.5e-5];
    elseif epsilon_KO %set translation to zero (hence no IkBe in system)
        p_mod = [
            30 1 0];
    elseif modify_IkBeRelA
        p_mod = [
            59 1 0.00003576
            61 1 0.00003576];
    elseif modify_IKKIkBe
        p_mod = [
            45 1 9.228
            46 1 9.228
            65 1 9.228];
    else
        p_mod = [];
    end

    doses = [25];
    dose_scale = 0.14185; % Convert to uM (from ug/ml)
    output_cpg = [];
    cd(model_path)
    for i = 1:length(doses)
        if isempty(output_cpg)
            [t,x,simdata] = nfkbSimulate({'CpG',doses(i)*dose_scale},names, p_mod, {},options);
        else
            options.STEADY_STATE = simdata.STEADY_STATE;
            [~,x] = nfkbSimulate({'CpG',doses(i)*dose_scale},names, p_mod, {},options);
        end
        output_cpg = cat(3,output_cpg,x);
    end
    cd(fpath)

    figure
    hold on

    %RelA Calculations
    rela_cyto = sum(output_cpg(:,1:5),2);
    rela_nuc = sum(output_cpg(:,6:8),2);
    rela(:,3) = rela_nuc;
    rela_total = (3.5*rela_cyto+rela_nuc)/4.5;
    rela_ratio(:,3) = rela_nuc/(rela_total(1)+untran_r);
    rela_ratio(:,3) = rela_ratio(:,3) - rela_ratio(1,3);

    %cRel Calculations
    crel_cyto = sum(output_cpg(:,9:13),2);
    crel_nuc = sum(output_cpg(:,14:16),2);
    crel(:,3) = crel_nuc;
    crel_total = (3.5*crel_cyto+crel_nuc)/4.5;
    crel_ratio(:,3) = crel_nuc/(crel_total(1)+untran_c);
    crel_ratio(:,3) = crel_ratio(:,3) - crel_ratio(1,3);
    
    title('CpG Simulation')
    xlabel('Time(hr)')
    ylabel('Nuclear/Total Cell Intensity')
    axis([0 options.SIM_TIME/60 0 2])

    plot(time, rela_ratio(:,3),'Color',sscanf('155098','%2x%2x%2x',[1 3])/255, 'LineWidth',4)
    plot(time, crel_ratio(:,3),'Color',sscanf('C4005B','%2x%2x%2x',[1 3])/255, 'LineWidth',4)
    
    hold off

end
%% Poly(I:C) Simulation
cd(fpath)
if run_poly
   
    if alpha_mut %reduced IkBa induced synthesis by Rela
        p_mod = [
            6 1 1.5e-5];
    elseif epsilon_KO %set translation to zero (hence no IkBe in system)
        p_mod = [
            30 1 0];
    elseif modify_IkBeRelA
        p_mod = [
            59 1 0.00003576
            61 1 0.00003576];
    elseif modify_IKKIkBe
        p_mod = [
            45 1 9.228
            46 1 9.228
            65 1 9.228];
    else
        p_mod = [];
    end

    doses = 1000*[20];
    dose_scale = 1/5e6; % Convert to uM. PolyI:C molecular weight: 1000KDa(+)
    output_poly = [];
    cd(model_path)
    for i = 1:length(doses)
        if isempty(output_poly)
            [t,x,simdata] = nfkbSimulate({'polyIC',doses(i)*dose_scale},names, p_mod, {},options);
        else
            options.STEADY_STATE = simdata.STEADY_STATE;
            [~,x] = nfkbSimulate({'polyIC',doses(i)*dose_scale},names, p_mod, {},options);
        end
        output_poly = cat(3,output_poly,x);
    end
    cd(fpath)

    figure
    hold on

    %RelA Calculations
    rela_cyto = sum(output_poly(:,1:5),2);
    rela_nuc = sum(output_poly(:,6:8),2);
    rela(:,4) = rela_nuc;
    rela_total = (3.5*rela_cyto+rela_nuc)/4.5;
    rela_ratio(:,4) = rela_nuc/(rela_total(1)+untran_r);
    rela_ratio(:,4) = rela_ratio(:,4) - rela_ratio(1,4);

    %cRel Calculations
    crel_cyto = sum(output_poly(:,9:13),2);
    crel_nuc = sum(output_poly(:,14:16),2);
    crel(:,4) = crel_nuc;
    crel_total = (3.5*crel_cyto+crel_nuc)/4.5;
    crel_ratio(:,4) = crel_nuc/(crel_total(1)+untran_c);
    crel_ratio(:,4) = crel_ratio(:,4) - crel_ratio(1,4);
    
    title('Poly(I:C) Simulation')
    xlabel('Time(hr)')
    ylabel('Nuclear/Total Cell Intensity')
    axis([0 options.SIM_TIME/60 0 2])
    
    plot(time, rela_ratio(:,4),'Color',sscanf('155098','%2x%2x%2x',[1 3])/255, 'LineWidth',4)
    plot(time, crel_ratio(:,4),'Color',sscanf('C4005B','%2x%2x%2x',[1 3])/255,'LineWidth',4)
    
    hold off

end
%% FLAGELLIN Simulation
cd(fpath)
if run_fla
   
   if alpha_mut %reduced IkBa induced synthesis by Rela
        p_mod = [
            6 1 1.5e-5];
    elseif epsilon_KO %set translation to zero (hence no IkBe in system)
        p_mod = [
            30 1 0];
    elseif modify_IkBeRelA
        p_mod = [
            59 1 0.00003576
            61 1 0.00003576];
    elseif modify_IKKIkBe
        p_mod = [
            45 1 9.228
            46 1 9.228
            65 1 9.228];
    else
        p_mod = [];
    end

    doses = [250]; %250 ng/mL; Flagellin = 50kDa
    dose_scale = 1/50000; % convert ng/mL to uM

    output_fla = [];
    cd(model_path)
    for i = 1:length(doses)
        if isempty(output_fla)
            [t,x,simdata] = nfkbSimulate({'FLA',doses(i)*dose_scale},names, p_mod, {},options);
        else
            options.STEADY_STATE = simdata.STEADY_STATE;
            [~,x] = nfkbSimulate({'FLA',doses(i)*dose_scale},names, p_mod, {},options);
        end
        output_fla = cat(3,output_fla,x);
    end
    cd(fpath)

    figure
    hold on

    %RelA Calculations
    rela_cyto = sum(output_fla(:,1:5),2);
    rela_nuc = sum(output_fla(:,6:8),2);
    rela(:,5) = rela_nuc;
    rela_total = (3.5*rela_cyto+rela_nuc)/4.5;
    rela_ratio(:,5) = rela_nuc/(rela_total(1)+untran_r);
    rela_ratio(:,5) = rela_ratio(:,5) - rela_ratio(1,5);

    %cRel Calculations
    crel_cyto = sum(output_fla(:,9:13),2);
    crel_nuc = sum(output_fla(:,14:16),2);
    crel(:,5) = crel_nuc;
    crel_total = (3.5*crel_cyto+crel_nuc)/4.5;
    crel_ratio(:,5) = crel_nuc/(crel_total(1)+untran_c);
    crel_ratio(:,5) = crel_ratio(:,5) - crel_ratio(1,5);
    
    title('Flagellin Simulation')
    xlabel('Time(hr)')
    ylabel('Nuclear/Total Cell Intensity')
    axis([0 options.SIM_TIME/60 0 2])

    plot(time, rela_ratio(:,5),'Color',sscanf('155098','%2x%2x%2x',[1 3])/255, 'LineWidth',4)
    plot(time, crel_ratio(:,5),'Color',sscanf('C4005B','%2x%2x%2x',[1 3])/255,'LineWidth',4)
    
    hold off

end
%% R848 Simulation
cd(fpath)
if run_r84

   if alpha_mut %reduced IkBa induced synthesis by Rela
        p_mod = [
            6 1 1.5e-5];
   elseif epsilon_KO %set translation to zero (hence no IkBe in system)
        p_mod = [
            30 1 0];
   elseif modify_IkBeRelA
        p_mod = [
            59 1 0.00003576
            61 1 0.00003576];
   elseif modify_IKKIkBe
        p_mod = [
            45 1 9.228
            46 1 9.228
            65 1 9.228];
   else
        p_mod = [];
   end

    doses = [350]; %350 ng/mL; R848 = 314 g/mol
    dose_scale = 1/314; % convert ng/mL to uM

    output_r84 = [];
    cd(model_path)
    for i = 1:length(doses)
        if isempty(output_r84)
            [t,x,simdata] = nfkbSimulate({'R848',doses(i)*dose_scale},names, p_mod, {},options);
        else
            options.STEADY_STATE = simdata.STEADY_STATE;
            [~,x] = nfkbSimulate({'R848',doses(i)*dose_scale},names, p_mod, {},options);
        end
        output_r84 = cat(3,output_r84,x);
    end
    cd(fpath)

    figure
    hold on

    %RelA Calculations
    rela_cyto = sum(output_r84(:,1:5),2);
    rela_nuc = sum(output_r84(:,6:8),2);
    rela(:,6) = rela_nuc;
    rela_total = (3.5*rela_cyto+rela_nuc)/4.5;
    rela_ratio(:,6) = rela_nuc/(rela_total(1)+untran_r);
    rela_ratio(:,6) = rela_ratio(:,6) - rela_ratio(1,6);

    %cRel Calculations
    crel_cyto = sum(output_r84(:,9:13),2);
    crel_nuc = sum(output_r84(:,14:16),2);
    crel(:,6) = crel_nuc;
    crel_total = (3.5*crel_cyto+crel_nuc)/4.5;
    crel_ratio(:,6) = crel_nuc/(crel_total(1)+untran_c);
    crel_ratio(:,6) = crel_ratio(:,6) - crel_ratio(1,6);
    
    title('R848 Simulation')
    xlabel('Time(hr)')
    ylabel('Nuclear/Total Cell Intensity')
    axis([0 options.SIM_TIME/60 0 2])
    
    plot(time, rela_ratio(:,6),'Color',sscanf('155098','%2x%2x%2x',[1 3])/255, 'LineWidth',4)
    plot(time, crel_ratio(:,6),'Color',sscanf('C4005B','%2x%2x%2x',[1 3])/255,'LineWidth',4)
    
    hold off

end
%% PAM3CSK Simulation
cd(fpath)
if run_pam
    
    if alpha_mut %reduced IkBa induced synthesis by Rela
        p_mod = [
            6 1 1.5e-5];
    elseif epsilon_KO %set translation to zero (hence no IkBe in system)
        p_mod = [
            30 1 0];
    elseif modify_IkBeRelA
        p_mod = [
            59 1 0.00003576
            61 1 0.00003576];
    elseif modify_IKKIkBe
        p_mod = [
            45 1 9.228
            46 1 9.228
            65 1 9.228];
    else
        p_mod = [];
    end

    doses = [40];
    dose_scale = 1/1500; % Convert to uM. Pam3CSK molecular weight: 1.5KDa
    length_t = 12;
    options.SIM_TIME = length_t*60;
    time = linspace(0,options.SIM_TIME/60,options.SIM_TIME+1);
    output_pam = [];
    cd(model_path)
    for i = 1:length(doses)
        if isempty(output_pam)
            [t,x,simdata] = nfkbSimulate({'Pam3CSK',doses(i)*dose_scale},names, p_mod, {},options);
        else
            options.STEADY_STATE = simdata.STEADY_STATE;
            [~,x] = nfkbSimulate({'Pam3CSK',doses(i)*dose_scale},names, p_mod, {},options);
        end
        output_pam = cat(3,output_pam,x);
    end
    cd(fpath)

    figure
    hold on

    %RelA Calculations
    rela_cyto = sum(output_pam(:,1:5),2);
    rela_nuc = sum(output_pam(:,6:8),2);
    rela(:,7) = rela_nuc;
    rela_cyto = [repmat(rela_cyto(1,1),90,1);rela_cyto];
    rela_nuc = [repmat(rela_nuc(1,1),90,1);rela_nuc];
    rela_cyto = rela_cyto(1:size(rela_ratio,1),:);
    rela_nuc = rela_nuc(1:size(rela_ratio,1),:);
    rela_total = (3.5*rela_cyto+rela_nuc)/4.5;
    rela_ratio(:,7) = rela_nuc/(rela_total(1)+untran_r);
    rela_ratio(:,7) = rela_ratio(:,7) - rela_ratio(1,7);

    %cRel Calculations
    crel_cyto = sum(output_pam(:,9:13),2);
    crel_nuc = sum(output_pam(:,14:16),2);
    crel(:,7) = crel_nuc;
    crel_cyto = [repmat(crel_cyto(1,1),90,1);crel_cyto];
    crel_nuc = [repmat(crel_nuc(1,1),90,1);crel_nuc];
    crel_cyto = crel_cyto(1:size(crel_ratio,1),:);
    crel_nuc = crel_nuc(1:size(crel_ratio,1),:);
    crel_total = (3.5*crel_cyto+crel_nuc)/4.5;
    crel_ratio(:,7) = crel_nuc/(crel_total(1)+untran_c);
    crel_ratio(:,7) = crel_ratio(:,7) - crel_ratio(1,7);
    
    title('Pam3CSK Simulation')
    xlabel('Time(hr)')
    ylabel('Nuclear/Total Cell Intensity')
    axis([0 options.SIM_TIME/60 0 2])

    plot(time, rela_ratio(:,7),'Color',sscanf('155098','%2x%2x%2x',[1 3])/255, 'LineWidth',4)
    plot(time, crel_ratio(:,7),'Color',sscanf('C4005B','%2x%2x%2x',[1 3])/255,'LineWidth',4)
    
    hold off
end
