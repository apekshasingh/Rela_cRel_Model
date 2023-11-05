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
run_cpg= true;
run_poly=true;
run_fla=true;
run_r84=true;
run_pam= true;

plot_order = [7, 1, 3, 4, 2, 5, 6];
%% Specify Plotting Options
plot_expt = true;
plot_rmsd = true;
plot_feats = true;
% RMSD options
relative_rmsd = false;
plot_points = false;
plot_errbar = true;

if plot_rmsd
    RMSD_vals = zeros(5, 2, 7);
end
if plot_feats
    feat_vals = zeros(6, 8, 7);
end
%% Plotting experimental data and model simulations
length_t = 12;

untran_r = 0.100;
untran_c = 0.025;

names = {'IkBaRelA','IkBeRelA','IKKIkBaRelA','IKKIkBeRelA','RelA','IkBaRelAn','IkBeRelAn','RelAn','IkBacRel','IkBecRel','IKKIkBacRel','IKKIkBecRel','cRel','IkBacReln','IkBecReln','cReln','IKK'};
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
p_mod = [];
%% TNF SIMULATION
cd(fpath)
if run_TNF
    tnf = readmatrix('tnf.csv');
    tnf(:,1) = [];
    tnf = tnf(tnf(:,1)<length_t,:);
    
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

    %RelA Calculations
    rela = [];
    rela_ratio = [];
    rela_cyto = sum(output_tnf(:,1:5),2);
    rela_nuc = sum(output_tnf(:,6:8),2);
    rela(:,1) = rela_nuc;
    rela_total = (3.5*rela_cyto+rela_nuc)/4.5;
    rela_ratio(:,1) = rela_nuc/(rela_total(1)+untran_r);

    %cRel Calculations
    crel = [];
    crel_ratio = [];
    crel_cyto = sum(output_tnf(:,9:13),2);
    crel_nuc = sum(output_tnf(:,14:16),2);
    crel(:,1) = crel_nuc;
    crel_total = (3.5*crel_cyto+crel_nuc)/4.5;
    crel_ratio(:,1) = crel_nuc/(crel_total(1)+untran_c);

    %Real Data
    if plot_expt
        for cell = 1:5
            plot(tnf(:,(3*(cell-1))+1),tnf(:,(3*(cell-1))+2),'Color', r_color(cell,:),'LineWidth',3)
            plot(tnf(:,(3*(cell-1))+1),tnf(:,(3*(cell-1))+3),'Color',c_color(cell,:),'LineWidth',3)
        end
    end

    %legend('RelA ratio (sim)','cRel ratio (sim)', 'RelA ratio', 'cRel ratio')
    title('TNF Simulation')
    xlabel('Time(hr)')
    ylabel('Nuclear/Total Cell Intensity')
    axis([0 options.SIM_TIME/60 0 1.5])
    
    plot(time, rela_ratio(:,1),'Color',sscanf('155098','%2x%2x%2x',[1 3])/255, 'LineWidth',4)
    plot(time, crel_ratio(:,1),'Color',sscanf('C4005B','%2x%2x%2x',[1 3])/255,'LineWidth',4)

    hold off
    
    if plot_rmsd
        for cell = 1:5
            RMSD_vals(cell, :, plot_order(1)) = calculate_RMSD(tnf(:,(3*(cell-1))+1),tnf(:,(3*(cell-1))+(2:3)),time, rela_ratio(:, 1), crel_ratio(:, 1), relative_rmsd);
        end
    end
    
    if plot_feats
        for cell = 1:5
            feat_vals(cell, :, plot_order(1)) = calculate_feats(tnf(:,(3*(cell-1))+1),tnf(:,(3*(cell-1))+(2:3)),time);
        end
        feat_vals(6, :, plot_order(1)) = calculate_feats(nan,[rela_ratio(:, 1), crel_ratio(:, 1)],time);
    end
end
%% LPS Simulation
cd(fpath)
if run_lps
    lps = readmatrix('lps.csv');
    lps(:,1) = [];
    lps = lps(lps(:,1)<length_t,:);

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

    %cRel Calculations
    crel_cyto = sum(output_lps(:,9:13),2);
    crel_nuc = sum(output_lps(:,14:16),2);
    crel(:,2) = crel_nuc;
    crel_total = (3.5*crel_cyto+crel_nuc)/4.5;
    crel_ratio(:,2) = crel_nuc/(crel_total(1)+untran_c);

    figure
    hold on

    %Real Data
    if plot_expt
        for cell = 1:5
            plot(lps(:,(3*(cell-1))+1),lps(:,(3*(cell-1))+2),'Color', r_color(cell,:),'LineWidth',3)
            hold on
            plot(lps(:,(3*(cell-1))+1),lps(:,(3*(cell-1))+3),'Color',c_color(cell,:),'LineWidth',3)
        end
    end
    
    title('LPS Simulation')
    xlabel('Time(hr)')
    ylabel('Nuclear/Total Cell Intensity')
    axis([0 options.SIM_TIME/60 0 2.5])

    plot(time, rela_ratio(:,2),'Color',sscanf('155098','%2x%2x%2x',[1 3])/255, 'LineWidth',4)
    plot(time, crel_ratio(:,2),'Color',sscanf('C4005B','%2x%2x%2x',[1 3])/255,'LineWidth',4)
    
    hold off
    if plot_rmsd
        for cell = 1:5
            RMSD_vals(cell, :, plot_order(2)) = calculate_RMSD(lps(:,(3*(cell-1))+1),lps(:,(3*(cell-1))+(2:3)),time, rela_ratio(:, 2), crel_ratio(:,2), relative_rmsd);
        end
    end
    
    if plot_feats
        for cell = 1:5
            feat_vals(cell, :, plot_order(2)) = calculate_feats(lps(:,(3*(cell-1))+1),lps(:,(3*(cell-1))+(2:3)),time);
        end
        feat_vals(6, :, plot_order(2)) = calculate_feats(nan,[rela_ratio(:, 2), crel_ratio(:, 2)],time);
    end
end
%% CpG Simulation
cd(fpath)
if run_cpg
    cpg = readmatrix('cpg.csv');
    cpg(:,1) = [];
    cpg = cpg(cpg(:,1)<length_t,:);

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

    %cRel Calculations
    crel_cyto = sum(output_cpg(:,9:13),2);
    crel_nuc = sum(output_cpg(:,14:16),2);
    crel(:,3) = crel_nuc;
    crel_total = (3.5*crel_cyto+crel_nuc)/4.5;
    crel_ratio(:,3) = crel_nuc/(crel_total(1)+untran_c);

    %Real Data
    if plot_expt
        for cell = 1:5
            plot(cpg(:,(3*(cell-1))+1),cpg(:,(3*(cell-1))+2),'Color', r_color(cell,:),'LineWidth',3)
            hold on
            plot(cpg(:,(3*(cell-1))+1),cpg(:,(3*(cell-1))+3),'Color',c_color(cell,:),'LineWidth',3)
        end
    end
    
    title('CpG Simulation')
    xlabel('Time(hr)')
    ylabel('Nuclear/Total Cell Intensity')
    axis([0 options.SIM_TIME/60 0 2.5])

    plot(time, rela_ratio(:,3),'Color',sscanf('155098','%2x%2x%2x',[1 3])/255, 'LineWidth',4)
    plot(time, crel_ratio(:,3),'Color',sscanf('C4005B','%2x%2x%2x',[1 3])/255, 'LineWidth',4)
    
    hold off
    if plot_rmsd
        for cell = 1:5
            RMSD_vals(cell, :, plot_order(3)) = calculate_RMSD(cpg(:,(3*(cell-1))+1),cpg(:,(3*(cell-1))+(2:3)),time, rela_ratio(:, 3), crel_ratio(:, 3), relative_rmsd);
        end
    end
    
    if plot_feats
        for cell = 1:5
            feat_vals(cell, :, plot_order(3)) = calculate_feats(cpg(:,(3*(cell-1))+1),cpg(:,(3*(cell-1))+(2:3)),time);
        end
        feat_vals(6, :, plot_order(3)) = calculate_feats(nan,[rela_ratio(:, 3), crel_ratio(:, 3)],time);
    end
end
%% Poly(I:C) Simulation
cd(fpath)
if run_poly
    poly = readmatrix('poly.csv');
    poly(:,1) = [];
    poly = poly(poly(:,1)<length_t,:);

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
    %rela_cyto = rela_cyto(1:size(rela_ratio,1),:);
    %rela_nuc = rela_nuc(1:size(rela_ratio,1),:);
    rela_total = (3.5*rela_cyto+rela_nuc)/4.5;
    rela_ratio(:,4) = rela_nuc/(rela_total(1)+untran_r);

    %cRel Calculations
    crel_cyto = sum(output_poly(:,9:13),2);
    crel_nuc = sum(output_poly(:,14:16),2);
    crel(:,4) = crel_nuc;
    %crel_cyto = crel_cyto(1:size(crel_ratio,1),:);
    %crel_nuc = crel_nuc(1:size(crel_ratio,1),:);
    crel_total = (3.5*crel_cyto+crel_nuc)/4.5;
    crel_ratio(:,4) = crel_nuc/(crel_total(1)+untran_c);

    %Real Data
    if plot_expt
        for cell = 1:5
            plot(poly(:,(3*(cell-1))+1),poly(:,(3*(cell-1))+2),'Color', r_color(cell,:),'LineWidth',3)
            hold on
            plot(poly(:,(3*(cell-1))+1),poly(:,(3*(cell-1))+3),'Color',c_color(cell,:),'LineWidth',3)
        end
    end
    
    title('Poly(I:C) Simulation')
    xlabel('Time(hr)')
    ylabel('Nuclear/Total Cell Intensity')
    axis([0 options.SIM_TIME/60 0 2])
    
    plot(time, rela_ratio(:,4),'Color',sscanf('155098','%2x%2x%2x',[1 3])/255, 'LineWidth',4)
    plot(time, crel_ratio(:,4),'Color',sscanf('C4005B','%2x%2x%2x',[1 3])/255,'LineWidth',4)
    
    hold off
    if plot_rmsd
        for cell = 1:5
            RMSD_vals(cell, :, plot_order(4)) = calculate_RMSD(poly(:,(3*(cell-1))+1),poly(:,(3*(cell-1))+(2:3)),time, rela_ratio(:, 4), crel_ratio(:, 4), relative_rmsd);
        end
    end
    
    if plot_feats
        for cell = 1:5
            feat_vals(cell, :, plot_order(4)) = calculate_feats(poly(:,(3*(cell-1))+1),poly(:,(3*(cell-1))+(2:3)),time);
        end
        feat_vals(6, :, plot_order(4)) = calculate_feats(nan,[rela_ratio(:, 4), crel_ratio(:, 4)],time);
    end
end
%% FLAGELLIN Simulation
cd(fpath)
if run_fla
    fla = readmatrix('fla.csv');
    fla(:,1) = [];
    fla = fla(fla(:,1)<length_t,:);

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

    %cRel Calculations
    crel_cyto = sum(output_fla(:,9:13),2);
    crel_nuc = sum(output_fla(:,14:16),2);
    crel(:,5) = crel_nuc;
    crel_total = (3.5*crel_cyto+crel_nuc)/4.5;
    crel_ratio(:,5) = crel_nuc/(crel_total(1)+untran_c);

    %Real Data
    if plot_expt
        for cell = 1:5
            plot(fla(:,(3*(cell-1))+1),fla(:,(3*(cell-1))+2),'Color', r_color(cell,:),'LineWidth',3)
            hold on
            plot(fla(:,(3*(cell-1))+1),fla(:,(3*(cell-1))+3),'Color',c_color(cell,:),'LineWidth',3)
        end
    end
    
    title('Flagellin Simulation')
    xlabel('Time(hr)')
    ylabel('Nuclear/Total Cell Intensity')
    axis([0 options.SIM_TIME/60 0 2.5])

    plot(time, rela_ratio(:,5),'Color',sscanf('155098','%2x%2x%2x',[1 3])/255, 'LineWidth',4)
    plot(time, crel_ratio(:,5),'Color',sscanf('C4005B','%2x%2x%2x',[1 3])/255,'LineWidth',4)
    
    hold off
    if plot_rmsd
        for cell = 1:5
            RMSD_vals(cell, :, plot_order(5)) = calculate_RMSD(fla(:,(3*(cell-1))+1),fla(:,(3*(cell-1))+(2:3)),time, rela_ratio(:, 5), crel_ratio(:, 5), relative_rmsd);
        end
    end
    
    if plot_feats
        for cell = 1:5
            feat_vals(cell, :, plot_order(5)) = calculate_feats(fla(:,(3*(cell-1))+1),fla(:,(3*(cell-1))+(2:3)),time);
        end
        feat_vals(6, :, plot_order(5)) = calculate_feats(nan,[rela_ratio(:, 5), crel_ratio(:, 5)],time);
    end
end
%% R848 Simulation
cd(fpath)
if run_r84
    r84 = readmatrix('r84.csv');
    r84(:,1) = [];
    r84 = r84(r84(:,1)<length_t,:);

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

    %cRel Calculations
    crel_cyto = sum(output_r84(:,9:13),2);
    crel_nuc = sum(output_r84(:,14:16),2);
    crel(:,6) = crel_nuc;
    crel_total = (3.5*crel_cyto+crel_nuc)/4.5;
    crel_ratio(:,6) = crel_nuc/(crel_total(1)+untran_c);

    %Real Data
    if plot_expt
        for cell = 1:5
            plot(r84(:,(3*(cell-1))+1),r84(:,(3*(cell-1))+2),'Color', r_color(cell,:),'LineWidth',3)
            hold on
            plot(r84(:,(3*(cell-1))+1),r84(:,(3*(cell-1))+3),'Color',c_color(cell,:),'LineWidth',3)
        end
    end
    
    title('R848 Simulation')
    xlabel('Time(hr)')
    ylabel('Nuclear/Total Cell Intensity')
    axis([0 options.SIM_TIME/60 0 2.5])
    
    plot(time, rela_ratio(:,6),'Color',sscanf('155098','%2x%2x%2x',[1 3])/255, 'LineWidth',4)
    plot(time, crel_ratio(:,6),'Color',sscanf('C4005B','%2x%2x%2x',[1 3])/255,'LineWidth',4)
    
    hold off
    if plot_rmsd
        for cell = 1:5
            RMSD_vals(cell, :, plot_order(6)) = calculate_RMSD(r84(:,(3*(cell-1))+1),r84(:,(3*(cell-1))+(2:3)),time, rela_ratio(:, 6), crel_ratio(:, 6), relative_rmsd);
        end
    end
    
     if plot_feats
        for cell = 1:5
            feat_vals(cell, :, plot_order(6)) = calculate_feats(r84(:,(3*(cell-1))+1),r84(:,(3*(cell-1))+(2:3)),time);
        end
        feat_vals(6, :, plot_order(6)) = calculate_feats(nan,[rela_ratio(:, 6), crel_ratio(:, 6)],time);
    end
end
%% PAM3CSK Simulation
cd(fpath)
if run_pam
    pam = readmatrix('pam.csv');
    pam(:,1) = [];
    pam = pam(pam(:,1)<length_t,:);

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
    rela_cyto = [repmat(rela_cyto(1,1),75,1);rela_cyto];
    rela_nuc = [repmat(rela_nuc(1,1),75,1);rela_nuc];
    rela_cyto = rela_cyto(1:size(rela_ratio,1),:);
    rela_nuc = rela_nuc(1:size(rela_ratio,1),:);
    rela_total = (3.5*rela_cyto+rela_nuc)/4.5;
    rela_ratio(:,7) = rela_nuc/(rela_total(1)+untran_r);

    %cRel Calculations
    crel_cyto = sum(output_pam(:,9:13),2);
    crel_nuc = sum(output_pam(:,14:16),2);
    crel(:,7) = crel_nuc;
    crel_cyto = [repmat(crel_cyto(1,1),75,1);crel_cyto];
    crel_nuc = [repmat(crel_nuc(1,1),75,1);crel_nuc];
    crel_cyto = crel_cyto(1:size(crel_ratio,1),:);
    crel_nuc = crel_nuc(1:size(crel_ratio,1),:);
    crel_total = (3.5*crel_cyto+crel_nuc)/4.5;
    crel_ratio(:,7) = crel_nuc/(crel_total(1)+untran_c);

    %Real Data
    if plot_expt
        for cell = 1:5
            plot(pam(:,(3*(cell-1))+1),pam(:,(3*(cell-1))+2),'Color', r_color(cell,:),'LineWidth',3)
            hold on
            plot(pam(:,(3*(cell-1))+1),pam(:,(3*(cell-1))+3),'Color',c_color(cell,:),'LineWidth',3)
        end
    end
    
    title('Pam3CSK Simulation')
    xlabel('Time(hr)')
    ylabel('Nuclear/Total Cell Intensity')
    axis([0 options.SIM_TIME/60 0 2.5])

    plot(time, rela_ratio(:,7),'Color',sscanf('155098','%2x%2x%2x',[1 3])/255, 'LineWidth',4)
    plot(time, crel_ratio(:,7),'Color',sscanf('C4005B','%2x%2x%2x',[1 3])/255,'LineWidth',4)
    
    hold off
    if plot_rmsd
        for cell = 1:5
            RMSD_vals(cell, :, plot_order(7)) = calculate_RMSD(pam(:,(3*(cell-1))+1),pam(:,(3*(cell-1))+(2:3)),time, rela_ratio(:, 7), crel_ratio(:, 7), relative_rmsd);
        end
    end
    
    if plot_feats
        for cell = 1:5
            feat_vals(cell, :, plot_order(7)) = calculate_feats(pam(:,(3*(cell-1))+1),pam(:,(3*(cell-1))+(2:3)),time);
        end
        feat_vals(6, :, plot_order(7)) = calculate_feats(nan,[rela_ratio(:, 7), crel_ratio(:, 7)],time);
    end
end
%% Plot RMSD of Model Simulations vs Experimental Measurements
if plot_rmsd
    figure()
    hold on
    mean_rmsd = squeeze(mean(RMSD_vals, 1));
    bar((1:7)-0.2, mean_rmsd(1, :), 0.4, 'FaceColor', sscanf('155098','%2x%2x%2x',[1 3])/255);
    bar((1:7)+0.2, mean_rmsd(2, :), 0.4, 'FaceColor', sscanf('C4005B','%2x%2x%2x',[1 3])/255);
    if plot_errbar
        t_val  = 2.776/sqrt(5);
        std_rmsd = squeeze(std(RMSD_vals, 1));
        er = errorbar((1:7)-0.2, mean_rmsd(1, :), t_val*std_rmsd(1, :), t_val*std_rmsd(1, :));
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';  
        er = errorbar((1:7)+0.2, mean_rmsd(2, :), t_val*std_rmsd(2, :), t_val*std_rmsd(2, :));
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';  
    end
    if plot_points
        for stim = 1:7
            scatter(repmat(stim-0.2, 1, 5), RMSD_vals(:,1,stim), 35, 'k', 'filled')
            scatter(repmat(stim+0.2, 1, 5), RMSD_vals(:,2,stim), 35, 'k', 'filled')
        end
    end
    if relative_rmsd
        ylabel("Relative RMSD")
        ylim([0, 3])
        xticklabels({[], "LPS", "Flagellin", "CpG", "Poly(I:C)", "R848", "Pam3CSK4", "TNF", []}) 
        xtickangle(45)
    else
        ylabel("RMSD")
        rng(1)
        index = randsample(7, 7);
        % generate shuffled RMSD
        shuff_RMSD = zeros(5, 2, 7);
        for cell = 1:5
            shuff_RMSD(cell, :, plot_order(1)) = calculate_RMSD(tnf(:,(3*(cell-1))+1),tnf(:,(3*(cell-1))+(2:3)),time, rela_ratio(:, index(1)), crel_ratio(:, index(1)), relative_rmsd);
            shuff_RMSD(cell, :, plot_order(2)) = calculate_RMSD(lps(:,(3*(cell-1))+1),lps(:,(3*(cell-1))+(2:3)),time, rela_ratio(:, index(2)), crel_ratio(:, index(2)), relative_rmsd);
            shuff_RMSD(cell, :, plot_order(3)) = calculate_RMSD(cpg(:,(3*(cell-1))+1),cpg(:,(3*(cell-1))+(2:3)),time, rela_ratio(:, index(3)), crel_ratio(:, index(3)), relative_rmsd);
            shuff_RMSD(cell, :, plot_order(4)) = calculate_RMSD(poly(:,(3*(cell-1))+1),poly(:,(3*(cell-1))+(2:3)),time, rela_ratio(:, index(4)), crel_ratio(:, index(4)), relative_rmsd);
            shuff_RMSD(cell, :, plot_order(5)) = calculate_RMSD(fla(:,(3*(cell-1))+1),fla(:,(3*(cell-1))+(2:3)),time, rela_ratio(:, index(5)), crel_ratio(:, index(5)), relative_rmsd);
            shuff_RMSD(cell, :, plot_order(6)) = calculate_RMSD(r84(:,(3*(cell-1))+1),r84(:,(3*(cell-1))+(2:3)),time, rela_ratio(:, index(6)), crel_ratio(:, index(6)), relative_rmsd);
            shuff_RMSD(cell, :, plot_order(7)) = calculate_RMSD(pam(:,(3*(cell-1))+1),pam(:,(3*(cell-1))+(2:3)),time, rela_ratio(:, index(7)), crel_ratio(:, index(7)), relative_rmsd);
        end
        mean_rmsd = [mean(squeeze(shuff_RMSD(:, 1, :)), 'all'), mean(squeeze(shuff_RMSD(:, 2, :)), 'all')];
        bar(8-0.2, mean_rmsd(1), 0.4, 'FaceColor', sscanf('155098','%2x%2x%2x',[1 3])/255);
        bar(8+0.2, mean_rmsd(2), 0.4, 'FaceColor', sscanf('C4005B','%2x%2x%2x',[1 3])/255);
        
        if plot_errbar
            t_val  = 1.96/sqrt(5*7);
            std_rmsd = [std(squeeze(shuff_RMSD(:, 1, :)), [], 'all'), std(squeeze(shuff_RMSD(:, 2, :)), [], 'all')];
            er = errorbar(8-0.2, mean_rmsd(1), t_val*std_rmsd(1), t_val*std_rmsd(1));
            er.Color = [0 0 0];                            
            er.LineStyle = 'none'; 
            er = errorbar(8+0.2, mean_rmsd(2), t_val*std_rmsd(2), t_val*std_rmsd(2));
            er.Color = [0 0 0];                            
            er.LineStyle = 'none';  
        end
        xticklabels({[], "LPS", "Flagellin", "CpG", "Poly(I:C)", "R848", "Pam3CSK4", "TNF", "shuffled"}) 
        xtickangle(45)
    end
    hold off
end
%% Plot difference in Model Simulations and Experimental Measurements Trajectory Features
if plot_feats
    figure()
    hold on
    %max <= 4 hrs
    for stim = 1:7
        x = repmat(stim-0.25, 1, 5)+(0.5*rand(1, 5));
        y = feat_vals(1:5, 1, stim);
        scatter(x, y, 35, r_color(:, 1:3), "filled", 'MarkerFaceAlpha',0.6)
        scatter(stim, feat_vals(6, 1, stim), 70, sscanf('155098','%2x%2x%2x',[1 3])/255, "filled")
        
        x = repmat(stim-0.25, 1, 5)+(0.5*rand(1, 5));
        y = feat_vals(1:5, 2, stim);
        scatter(x, y, 35, c_color(:, 1:3), "filled", 'MarkerFaceAlpha',0.6)
        scatter(stim, feat_vals(6, 2, stim), 70, sscanf('C4005B','%2x%2x%2x',[1 3])/255, "filled")
    end
    hold off
    ylabel("Peak Amp")
    xticklabels({[], "LPS", "Flagellin", "CpG", "Poly(I:C)", "R848", "Pam3CSK4", "TNF", []}) 
    xtickangle(45)
    
%     figure()
%     hold on
    %time to peak
%     for stim = 1:7
%         x = repmat(stim-0.25, 1, 5)+(0.5*rand(1, 5));
%         y = feat_vals(1:5, 3, stim);
%         scatter(x, y, 35, r_color(:, 1:3), "filled", 'MarkerFaceAlpha',0.6)
%         scatter(stim, feat_vals(6, 3, stim), 70, sscanf('155098','%2x%2x%2x',[1 3])/255, "filled")
%         
%         x = repmat(stim-0.25, 1, 5)+(0.5*rand(1, 5));
%         y = feat_vals(1:5, 4, stim);
%         scatter(x, y, 35, c_color(:, 1:3), "filled", 'MarkerFaceAlpha',0.6)
%         scatter(stim, feat_vals(6, 4, stim), 70, sscanf('C4005B','%2x%2x%2x',[1 3])/255, "filled")
%     end
%     hold off
%     ylabel("Peak Time")
%     xticklabels({[],  "LPS", "Flagellin", "CpG", "Poly(I:C)", "R848", "Pam3CSK4", "TNF", []}) 
%     xtickangle(45)
    
    figure()
    hold on
    %max <= 4 hrs
    for stim = 1:7
        x = repmat(stim-0.25, 1, 5)+(0.5*rand(1, 5));
        y = feat_vals(1:5, 5, stim);
        scatter(x, y, 35, r_color(:, 1:3), "filled", 'MarkerFaceAlpha',0.6)
        scatter(stim, feat_vals(6, 5, stim), 70, sscanf('155098','%2x%2x%2x',[1 3])/255, "filled")
        
        x = repmat(stim-0.25, 1, 5)+(0.5*rand(1, 5));
        y = feat_vals(1:5, 6, stim);
        scatter(x, y, 35, c_color(:, 1:3), "filled", 'MarkerFaceAlpha',0.6)
        scatter(stim, feat_vals(6, 6, stim), 70, sscanf('C4005B','%2x%2x%2x',[1 3])/255, "filled")
    end
    hold off
    ylabel("Activity")
    xticklabels({[], "LPS", "Flagellin", "CpG", "Poly(I:C)", "R848", "Pam3CSK4", "TNF", []}) 
    xtickangle(45)
    
%     figure()
%     hold on
    %time to half max
%     for stim = 1:7
%         x = repmat(stim-0.25, 1, 5)+(0.5*rand(1, 5));
%         y = feat_vals(1:5, 7, stim);
%         scatter(x, y, 35, r_color(:, 1:3), "filled", 'MarkerFaceAlpha',0.6)
%         scatter(stim, feat_vals(6, 7, stim), 70, sscanf('155098','%2x%2x%2x',[1 3])/255, "filled")
%         
%         x = repmat(stim-0.25, 1, 5)+(0.5*rand(1, 5));
%         y = feat_vals(1:5, 8, stim);
%         scatter(x, y, 35, c_color(:, 1:3), "filled", 'MarkerFaceAlpha',0.6)
%         scatter(stim, feat_vals(6, 8, stim), 70, sscanf('C4005B','%2x%2x%2x',[1 3])/255, "filled")
%     end
%     hold off
%     ylabel("Time Half Max")
%     xticklabels({[],  "LPS", "Flagellin", "CpG", "Poly(I:C)", "R848", "Pam3CSK4", "TNF", []}) 
%     xtickangle(45)
end
