%% Get current directory
fpath = mfilename('fullpath');
search = strfind(fpath, '/');
if isempty(search)
    search = strfind(fpath, '\');
end
fpath = fpath(1:search(end));
cd(fpath)
%% Specify Stimulation Conditions
run_TNF = true;
run_lps = true;
run_cpg= true;
run_pam= true;
run_poly= true;
run_fla= true;
run_r84= true;
%% Specify Paramater Modifications
modify_IkBeRelA = false;
modify_IKKIkBe = true;
%% Plotting simulations with parameter modifications
length_t = 12;

untran_r = 0.100;
untran_c = 0.025;

names = {'IkBaRelA','IkBeRelA','IKKIkBaRelA','IKKIkBeRelA','RelA','IkBaRelAn','IkBeRelAn','RelAn','IkBacRel','IkBecRel','IKKIkBacRel','IKKIkBecRel','cRel','IkBacReln','IkBecReln','cReln','IKK'};
options = struct;
options.DEBUG = 1;
options.SIM_TIME = length_t*60;
time = linspace(0,options.SIM_TIME/60,options.SIM_TIME+1);

%set colors
rela_color = sscanf('155098','%2x%2x%2x',[1 3])/255;
crel_color = sscanf('C4005B','%2x%2x%2x',[1 3])/255;

%set modified parameter values
if modify_IkBeRelA %original dissociation 0.08 (KD = 0.0014 to KD = 6.3937e-07)
    mod_vals = linspace(log(0.08), log(0.00003576), 4);
    dp = mod_vals(2)-mod_vals(1);
    mod_vals = exp([mod_vals(1)-dp, mod_vals]);
    alpha = linspace(1, 0.05, length(mod_vals));
elseif modify_IKKIkBe %original dissociation 38 (KD = 0.8235 to KD = 0.20)
    mod_vals = linspace(log(38), log(9.228), 4);
    dp = mod_vals(2)-mod_vals(1);
    mod_vals = exp([mod_vals(1)-dp, mod_vals, mod_vals(length(mod_vals))+dp]);
    alpha = linspace(1, 0.2, length(mod_vals));
else
    mod_vals = nan;
    alpha = 1
end

% Navigate to model script folder
model_path = strcat(fpath(1:search(end-1)), "Model_Script/");
cd(model_path)
%% TNF SIMULATION
if run_TNF
    figure()
    for p = 1:length(mod_vals)
        if modify_IkBeRelA
            p_mod = [
                59 1 mod_vals(p)
                61 1 mod_vals(p)];
        elseif modify_IKKIkBe
            p_mod = [
                45 1 mod_vals(p)
                46 1 mod_vals(p)
                65 1 mod_vals(p)];
        else
            p_mod = [];
        end

        dose_scale = 1/5200; % Dose scaling from Shannon et al (2007)
        doses = [10];

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

        subplot(2, 1, 1)
        plot(time, rela_ratio(:,1),'Color',[rela_color, alpha(p)], 'LineWidth',4)
        hold on
        subplot(2, 1, 2)
        plot(time, crel_ratio(:,1),'Color',[crel_color, alpha(p)],'LineWidth',4)
        hold on
        
        if p ==2
            subplot(2, 1, 1)
            plot(time, rela_ratio(:,1),'Color','k', 'LineWidth',2,'LineStyle', "--")
            subplot(2, 1, 2)
            plot(time, crel_ratio(:,1),'Color','k', 'LineWidth',2,'LineStyle', "--")
        end
        if p == length(mod_vals)-1
            subplot(2, 1, 1)
            plot(time, rela_ratio(:,1),'Color','k', 'LineWidth',2,'LineStyle', ":")
            subplot(2, 1, 2)
            plot(time, crel_ratio(:,1),'Color','k', 'LineWidth',2,'LineStyle', ":")
        end
    end
    ax1 = subplot(2, 1, 1);
    title('TNF Simulation')
    ax2 = subplot(2, 1, 2);
    linkaxes([ax1, ax2], 'xy')
    xlabel('Time(hr)')
    ylabel('Nuclear/Total Cell Intensity')
    %axis([0 options.SIM_TIME/60 0 1.5])
    hold off
end
%% LPS SIMULATION
if run_lps
    figure()
    for p = 1:length(mod_vals)
        if modify_IkBeRelA
            p_mod = [
                59 1 mod_vals(p)
                61 1 mod_vals(p)];
        elseif modify_IKKIkBe
            p_mod = [
                45 1 mod_vals(p)
                46 1 mod_vals(p)
                65 1 mod_vals(p)];
        else
            p_mod = [];
        end

        doses = [10];
        dose_scale = 1/24000; % LPS molecular weight estimated between 10KDa and 100KDa

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

        subplot(2, 1, 1)
        plot(time, rela_ratio(:,2),'Color',[rela_color, alpha(p)], 'LineWidth',4)
        hold on
        subplot(2, 1, 2)
        plot(time, crel_ratio(:,2),'Color',[crel_color, alpha(p)],'LineWidth',4)
        hold on
        
        if p ==2
            subplot(2, 1, 1)
            plot(time, rela_ratio(:,2),'Color','k', 'LineWidth',2,'LineStyle', "--")
            subplot(2, 1, 2)
            plot(time, crel_ratio(:,2),'Color','k', 'LineWidth',2,'LineStyle', "--")
        end
        if p == length(mod_vals)-1
            subplot(2, 1, 1)
            plot(time, rela_ratio(:,2),'Color','k', 'LineWidth',2,'LineStyle', ":")
            subplot(2, 1, 2)
            plot(time, crel_ratio(:,2),'Color','k', 'LineWidth',2,'LineStyle', ":")
        end
    end
    ax1 = subplot(2, 1, 1);
    title('LPS Simulation')
    ax2 = subplot(2, 1, 2);
    linkaxes([ax1, ax2], 'xy')
    xlabel('Time(hr)')
    ylabel('Nuclear/Total Cell Intensity')
    %axis([0 options.SIM_TIME/60 0 2.5])
    hold off
end
%% CpG Simulation
if run_cpg
    figure()
    for p = 1:length(mod_vals)
        if modify_IkBeRelA
            p_mod = [
                59 1 mod_vals(p)
                61 1 mod_vals(p)];
        elseif modify_IKKIkBe
            p_mod = [
                45 1 mod_vals(p)
                46 1 mod_vals(p)
                65 1 mod_vals(p)];
        else
            p_mod = [];
        end
        doses = [25];
        dose_scale = 0.14185; % Convert to uM (from ug/ml)
        output_cpg = [];
        for i = 1:length(doses)
            if isempty(output_cpg)
                [t,x,simdata] = nfkbSimulate({'CpG',doses(i)*dose_scale},names, p_mod, {},options);
            else
                options.STEADY_STATE = simdata.STEADY_STATE;
                [~,x] = nfkbSimulate({'CpG',doses(i)*dose_scale},names, p_mod, {},options);
            end
            output_cpg = cat(3,output_cpg,x);
        end

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
        
        subplot(2, 1, 1)
        plot(time, rela_ratio(:,3),'Color',[rela_color, alpha(p)], 'LineWidth',4)
        hold on
        subplot(2, 1, 2)
        plot(time, crel_ratio(:,3),'Color',[crel_color, alpha(p)],'LineWidth',4)
        hold on
        
        if p ==2
            subplot(2, 1, 1)
            plot(time, rela_ratio(:,3),'Color','k', 'LineWidth',2,'LineStyle', "--")
            subplot(2, 1, 2)
            plot(time, crel_ratio(:,3),'Color','k', 'LineWidth',2,'LineStyle', "--")
        end
        if p == length(mod_vals)-1
            subplot(2, 1, 1)
            plot(time, rela_ratio(:,3),'Color','k', 'LineWidth',2,'LineStyle', ":")
            subplot(2, 1, 2)
            plot(time, crel_ratio(:,3),'Color','k', 'LineWidth',2,'LineStyle', ":")
        end
    end
    ax1 = subplot(2, 1, 1);
    title('CpG Simulation')
    ax2 = subplot(2, 1, 2);
    linkaxes([ax1, ax2], 'xy')
    xlabel('Time(hr)')
    ylabel('Nuclear/Total Cell Intensity')
    %axis([0 options.SIM_TIME/60 0 2.5])
    hold off
end
%% Poly(I:C) Simulation
if run_poly
    figure()
    for p = 1:length(mod_vals)
        if modify_IkBeRelA
            p_mod = [
                59 1 mod_vals(p)
                61 1 mod_vals(p)];
        elseif modify_IKKIkBe
            p_mod = [
                45 1 mod_vals(p)
                46 1 mod_vals(p)
                65 1 mod_vals(p)];
        else
            p_mod = [];
        end

        doses = 1000*[20];
        dose_scale = 1/5e6; % Convert to uM. PolyI:C molecular weight: 1000KDa(+)
        output_poly = [];
        for i = 1:length(doses)
            if isempty(output_poly)
                [t,x,simdata] = nfkbSimulate({'polyIC',doses(i)*dose_scale},names, p_mod, {},options);
            else
                options.STEADY_STATE = simdata.STEADY_STATE;
                [~,x] = nfkbSimulate({'polyIC',doses(i)*dose_scale},names, p_mod, {},options);
            end
            output_poly = cat(3,output_poly,x);
        end

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
    
        subplot(2, 1, 1)
        plot(time, rela_ratio(:,4),'Color',[rela_color, alpha(p)], 'LineWidth',4)
        hold on
        subplot(2, 1, 2)
        plot(time, crel_ratio(:,4),'Color',[crel_color, alpha(p)],'LineWidth',4)
        hold on
        
        if p ==2
            subplot(2, 1, 1)
            plot(time, rela_ratio(:,4),'Color','k', 'LineWidth',2,'LineStyle', "--")
            subplot(2, 1, 2)
            plot(time, crel_ratio(:,4),'Color','k', 'LineWidth',2,'LineStyle', "--")
        end
        if p == length(mod_vals)-1
            subplot(2, 1, 1)
            plot(time, rela_ratio(:,4),'Color','k', 'LineWidth',2,'LineStyle', ":")
            subplot(2, 1, 2)
            plot(time, crel_ratio(:,4),'Color','k', 'LineWidth',2,'LineStyle', ":")
        end
    end
    ax1 = subplot(2, 1, 1);
    title('Poly(I:C) Simulation')
    ax2 = subplot(2, 1, 2);
    linkaxes([ax1, ax2], 'xy')
    xlabel('Time(hr)')
    ylabel('Nuclear/Total Cell Intensity')
    %axis([0 options.SIM_TIME/60 0 2])
    hold off
end
%% FLAGELLIN Simulation
if run_fla
   figure()
    for p = 1:length(mod_vals)
        if modify_IkBeRelA
            p_mod = [
                59 1 mod_vals(p)
                61 1 mod_vals(p)];
        elseif modify_IKKIkBe
            p_mod = [
                45 1 mod_vals(p)
                46 1 mod_vals(p)
                65 1 mod_vals(p)];
        else
            p_mod = [];
        end

        doses = [250]; %250 ng/mL; Flagellin = 50kDa
        dose_scale = 1/50000; % convert ng/mL to uM

        output_fla = [];
        for i = 1:length(doses)
            if isempty(output_fla)
                [t,x,simdata] = nfkbSimulate({'FLA',doses(i)*dose_scale},names, p_mod, {},options);
            else
                options.STEADY_STATE = simdata.STEADY_STATE;
                [~,x] = nfkbSimulate({'FLA',doses(i)*dose_scale},names, p_mod, {},options);
            end
            output_fla = cat(3,output_fla,x);
        end

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

        subplot(2, 1, 1)
        plot(time, rela_ratio(:,5),'Color',[rela_color, alpha(p)], 'LineWidth',4)
        hold on
        subplot(2, 1, 2)
        plot(time, crel_ratio(:,5),'Color',[crel_color, alpha(p)],'LineWidth',4)
        hold on
        
        if p ==2
            subplot(2, 1, 1)
            plot(time, rela_ratio(:,5),'Color','k', 'LineWidth',2,'LineStyle', "--")
            subplot(2, 1, 2)
            plot(time, crel_ratio(:,5),'Color','k', 'LineWidth',2,'LineStyle', "--")
        end
        if p == length(mod_vals)-1
            subplot(2, 1, 1)
            plot(time, rela_ratio(:,5),'Color','k', 'LineWidth',2,'LineStyle', ":")
            subplot(2, 1, 2)
            plot(time, crel_ratio(:,5),'Color','k', 'LineWidth',2,'LineStyle', ":")
        end
    end
    ax1 = subplot(2, 1, 1);
    title('Flagellin Simulation')
    ax2 = subplot(2, 1, 2);
    linkaxes([ax1, ax2], 'xy')
    xlabel('Time(hr)')
    ylabel('Nuclear/Total Cell Intensity')
    %axis([0 options.SIM_TIME/60 0 2.5])
    hold off
end
%% R848 Simulation
if run_r84
   figure()
    for p = 1:length(mod_vals)
        if modify_IkBeRelA
            p_mod = [
                59 1 mod_vals(p)
                61 1 mod_vals(p)];
        elseif modify_IKKIkBe
            p_mod = [
                45 1 mod_vals(p)
                46 1 mod_vals(p)
                65 1 mod_vals(p)];
        else
            p_mod = [];
        end

        doses = [350]; %350 ng/mL; R848 = 314 g/mol
        dose_scale = 1/314; % convert ng/mL to uM

        output_r84 = [];
        for i = 1:length(doses)
            if isempty(output_r84)
                [t,x,simdata] = nfkbSimulate({'R848',doses(i)*dose_scale},names, p_mod, {},options);
            else
                options.STEADY_STATE = simdata.STEADY_STATE;
                [~,x] = nfkbSimulate({'R848',doses(i)*dose_scale},names, p_mod, {},options);
            end
            output_r84 = cat(3,output_r84,x);
        end

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
    
        subplot(2, 1, 1)
        plot(time, rela_ratio(:,6),'Color',[rela_color, alpha(p)], 'LineWidth',4)
        hold on
        subplot(2, 1, 2)
        plot(time, crel_ratio(:,6),'Color',[crel_color, alpha(p)],'LineWidth',4)
        hold on
        
        if p ==2
            subplot(2, 1, 1)
            plot(time, rela_ratio(:,6),'Color','k', 'LineWidth',2,'LineStyle', "--")
            subplot(2, 1, 2)
            plot(time, crel_ratio(:,6),'Color','k', 'LineWidth',2,'LineStyle', "--")
        end
        if p == length(mod_vals)-1
            subplot(2, 1, 1)
            plot(time, rela_ratio(:,6),'Color','k', 'LineWidth',2,'LineStyle', ":")
            subplot(2, 1, 2)
            plot(time, crel_ratio(:,6),'Color','k', 'LineWidth',2,'LineStyle', ":")
        end
    end
    ax1 = subplot(2, 1, 1);
    title('R848 Simulation')
    ax2 = subplot(2, 1, 2);
    linkaxes([ax1, ax2], 'xy')
    xlabel('Time(hr)')
    ylabel('Nuclear/Total Cell Intensity')
    %axis([0 options.SIM_TIME/60 0 2.5])
    hold off
end
%% PAM3CSK Simulation
if run_pam
    figure()
    for p = 1:length(mod_vals)
        if modify_IkBeRelA
            p_mod = [
                59 1 mod_vals(p)
                61 1 mod_vals(p)];
        elseif modify_IKKIkBe
            p_mod = [
                45 1 mod_vals(p)
                46 1 mod_vals(p)
                65 1 mod_vals(p)];
        else
            p_mod = [];
        end
        doses = [40];
        dose_scale = 1/1500; % Convert to uM. Pam3CSK molecular weight: 1.5KDa
        length_t = 12;
        options.SIM_TIME = length_t*60;
        time = linspace(0,options.SIM_TIME/60,options.SIM_TIME+1);
        output_pam = [];
        for i = 1:length(doses)
            if isempty(output_pam)
                [t,x,simdata] = nfkbSimulate({'Pam3CSK',doses(i)*dose_scale},names, p_mod, {},options);
            else
                options.STEADY_STATE = simdata.STEADY_STATE;
                [~,x] = nfkbSimulate({'Pam3CSK',doses(i)*dose_scale},names, p_mod, {},options);
            end
            output_pam = cat(3,output_pam,x);
        end

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
    
        subplot(2, 1, 1)
        plot(time, rela_ratio(:,7),'Color',[rela_color, alpha(p)], 'LineWidth',4)
        hold on
        subplot(2, 1, 2)
        plot(time, crel_ratio(:,7),'Color',[crel_color, alpha(p)],'LineWidth',4)
        hold on
        
        if p ==2
            subplot(2, 1, 1)
            plot(time, rela_ratio(:,7),'Color','k', 'LineWidth',2,'LineStyle', "--")
            subplot(2, 1, 2)
            plot(time, crel_ratio(:,7),'Color','k', 'LineWidth',2,'LineStyle', "--")
        end
        if p == length(mod_vals)-1
            subplot(2, 1, 1)
            plot(time, rela_ratio(:,7),'Color','k', 'LineWidth',2,'LineStyle', ":")
            subplot(2, 1, 2)
            plot(time, crel_ratio(:,7),'Color','k', 'LineWidth',2,'LineStyle', ":")
        end
    end
    ax1 = subplot(2, 1, 1);
    title('Pam3CSK Simulation')
    ax2 = subplot(2, 1, 2);
    linkaxes([ax1, ax2], 'xy')
    xlabel('Time(hr)')
    ylabel('Nuclear/Total Cell Intensity')
    %axis([0 options.SIM_TIME/60 0 2.5])
    hold off
end