function aurora_ppg_quality_assessment2
%
% AURORA_PPG_QUALITY_ASSESSMENT
%
% Code for analysing data for the following publication:
%    Charlton PH et al., 'Determinants of photoplethysmography signal
%    quality at the wrist'.
%

%% setup universal parameters
up = setup_up;

%% Cycle through each protocol
for protocol_no = 1 : length(up.protocols)

    clear data

    curr_protocol = up.protocols{protocol_no};

    %% setup universal parameters for this protocol
    up = setup_up(curr_protocol);
    
    %% Extract data from dataset and analyse signals
    % Extract the data from the dataset, analyse signals, and save to file.
    % Or, if this has already been done, load from the file.
    data = extract_data_from_dataset_and_analyse_signals(up);
    
    %% Extract subject characteristics
    % Print a latex table of subject characteristicsg
    extract_subj_chars(data, up);

    %% Extract data recording characteristics
    % Print recording characteristics
    extract_recording_chars(data, up);

    %% Report levels of signal quality
    % Print a summary of the levels of signal quality, make histograms of
    % levels of signal quality, and print a latex table.
    report_levels_signal_quality(data, up);

    %% linear mixed effects model
    % Perform mixed effects analysis, including creating plots to visualise
    % the associations between factors and signal quality
    linear_mixed_effects_analysis(data, up);

    %% Make plots comparing postures
    make_plots_comparing_postures(data, up);

    %% Create demo figures
    create_all_demo_figures(data, up);

    clear data

end

end

function data = load_extracted_data(up)

fprintf('\n - Loading converted data')
load(up.conv_data_filepath);

end

function make_plots_comparing_postures(data, up)

% only do this for the oscillometric protocol
if ~strcmp(up.protocol, 'oscillometric')
    return
end

% identify postures
postures = unique(data.posture);

% cycle through each signal quality metric
all_p = [];
for sig_char_no = 1 : length(up.signal_quality_metrics)
    curr_sig_char = up.signal_quality_metrics{sig_char_no};

    % plot comparing different postures
    rel_hts = {'up', 'natural'};
    for ht_no = 1 : length(rel_hts)
        rel_arm_height = rel_hts{ht_no};
        all_p = make_plot_comparing_postures(data, curr_sig_char, postures, rel_arm_height, all_p, up);
    end

    % plots comparing different sensor heights for each posture
    for posture_no = 1 : length(postures)
        posture = string(postures(posture_no));
        all_p = make_plot_at_certain_posture(data, curr_sig_char, posture, all_p, up);
    end

end

% check whether differences remain significant after correction for
% multiple comparisons:
sig_diffs = correct_multiple_comparisons(all_p, up);
if sum(~sig_diffs)
    warning('significant differences werent all significant after correction for multiple comparisons')
end

% To control for dc_amp repeat with: rel_els = data.dc_amp>-4e5 & data.dc_amp<-2e5; data = data(rel_els,:);

end

function sig_diffs = correct_multiple_comparisons(p, up)

% Adapted from:
% https://github.com/peterhcharlton/RRest/blob/master/RRest_v3.0/Publication_Specific_Scripts/run_vortal_determinants_analysis.m
% (as originally written by me, PH Charlton, it can be resued without a GNU GPL)

% extract required data, and order according to p-value (ascending)
[p_sort, order] = sort(p);

% setup variable
temp_sig_diffs = false(length(p_sort),1);

% find significant differences
for comparison_no = 1 : length(p_sort)
    % Holm's correction for only the number of remaining comparisons
    no_tests = length(p_sort)-comparison_no+1;
    % Sidak's correction for the number of remaining comparisons
    alpha_sidak = 1 - ((1-up.alpha)^(1/no_tests));
    if p_sort(comparison_no) < alpha_sidak
        temp_sig_diffs(comparison_no) = true;
    else
        break
    end
end

% return significant differences back to original order
orig_p = zeros(size(p));
orig_p(order) = p_sort;

sig_diffs(order) = temp_sig_diffs;

% check
if ~isequal(p, orig_p)
    error('these should be equal')
end

end

function data = extract_data_from_dataset_and_analyse_signals(up)

if ~exist(up.conv_data_filepath, 'file') || up.redo_extraction

    %% load data
    data = load_data(up);

    %% Identify subjects for inclusion
    % only include these subjects in the analysis
    data = identify_subjs_for_inclusion(data, up);

    %% Extract signal characteristics
    data = extract_sig_chars(data, up);

    %% Import BPs
    % need to note that excluded BPs with high error
    data = import_bps(data, up);

    %% save data
    fprintf('\n - Saving converted data')
    save(up.conv_data_filepath, 'data');

else

    data = load_extracted_data(up);
    
end

end

function linear_mixed_effects_analysis(data, up)
% Performs mixed effects analysis for each signal quality metric, and with
% and without dc amplitude as a fixed effect.

fprintf('\n\n ~~~~~~~~ ')
fprintf('\n - Performing mixed effects analysis');
fprintf('\n ~~~~~~~~ \n')

% cycle through each signal quality metric
for ref_var_no = 1 : length(up.signal_quality_metrics)
    curr_ref_var = up.signal_quality_metrics{ref_var_no};

    % Repeat with and without including dc amplitude as a fixed effect
    for include_dc_amp = 0:1
        
        % Perform mixed effects analysis (including creating plots if dc
        % amplitude is included)
        make_mixed_effects_models(data, curr_ref_var, include_dc_amp, up);

    end
end

end

function create_all_demo_figures(data, up)

if strcmp(up.protocol, 'oscillometric')
    if ~exist('data', 'var')
        fprintf('\n - Loading converted data')
        load(up.conv_data_filepath);
    end
    create_ppg_quality_figure(data, up);
    create_exceptions_ppg_quality_figure(data, up);
    create_ppg_quality_figure_age(data, up);
end
% create_demo_figures(data, up); - not used

end

function report_levels_signal_quality(data, up)

fprintf('\n\n ~~~~~~~~ ')
fprintf('\n - Reporting levels of signal quality for %s protocol', up.protocol);
fprintf('\n ~~~~~~~~ \n')

% - no participants
n_pts = length(unique(data.subj_id));
fprintf('\n - %d Participants', n_pts)

% - recordings per participant
[unique_subjs, ~, J]=unique(data.subj_id);
occ = histc(J, 1:numel(unique_subjs));
no_recs_per_pt = unique(occ);
if length(no_recs_per_pt) == 1
    fprintf('\n - %d recordings per participant:', no_recs_per_pt)
else
    error('this shouldn''t happen')
end
[unique_activities, ~, J]=unique(data.activity);
occ = histc(J, 1:numel(unique_activities))./n_pts;
for activity_no = 1 : length(unique_activities)
    fprintf('\n   - %d recordings per participant for %s activity', occ(activity_no), unique_activities{activity_no});
end

metrics = {'snr', 'ac_dc_ratio', 'tm_cc'};
for metric_no = 1 : length(metrics)
    eval(['rel_data = data.' metrics{metric_no} ';']);
    fprintf('\n There was a wide range of %s, from %.2f to %.2f.', metrics{metric_no}, min(rel_data), max(rel_data));
    fprintf('\n   with typical values of %.2f to %.2f.', quantile(rel_data, 0.25), quantile(rel_data, 0.75));

    % create histogram
    create_histogram(rel_data, metrics{metric_no}, up);

end

% - signal quality at different positions
%activities = {'standingarmup', 'standingarmdown', 'sittingarmup', 'sittingarmlap', 'sittingarmdown', 'supine'};
for activity_no = 1 : length(unique_activities)
    curr_activity = unique_activities{activity_no};
    for metric_no = 1 : length(metrics)
        [med(activity_no, metric_no), lq(activity_no, metric_no), uq(activity_no, metric_no), ...
            mu(activity_no, metric_no), sd(activity_no, metric_no), min_val(activity_no, metric_no), max_val(activity_no, metric_no)] ...
            = print_metric(data, curr_activity, metrics{metric_no});
    end
end

output_observed_signal_quality_table(metrics, unique_activities, med, lq, uq, mu, sd, up);


% assess colinearity

% extract matrix of data
var_data = [];
for metric_no = 1 : length(metrics)
    curr_var = metrics{metric_no};
    eval(['var_data(:, metric_no) = data.' curr_var ';']);
end

% exclude rows with nans
var_data = var_data(sum(isnan(var_data),2)==0,:);

% calculate correlations between variables
corrMatrix = corr(var_data);
rsquaredMatrix = corrMatrix.^2;

end

function create_histogram(rel_data, metric_name, up)

ftsize = 24;
figure('Position', [20, 20, 700, 400])

histogram(rel_data)

% labels
set(gca, 'FontSize', ftsize-4, 'YGrid', 'on')
if strcmp(metric_name, 'snr')
    xlabel('SNR (dB)', 'FontSize', ftsize)
    xlim([-10 65])
elseif strcmp(metric_name, 'ac_dc_ratio')
    xlabel('PI (%)', 'FontSize', ftsize)
    xlim([0 7])
elseif strcmp(metric_name, 'tm_cc')
    xlabel('TMCC', 'FontSize', ftsize)
    xlim([-0.2 1])
end
ylim([0 500])
ylabel('No. recordings', 'FontSize', ftsize)

% tidy up
box off

% save figure
filepath = [up.plots_path, 'quality_metric_histogram_', metric_name, '_', up.protocol];
save_figure(filepath)

end

function output_observed_signal_quality_table(metrics, activities, med, lq, uq, mu, sd, up)

fprintf('\n - Observed signal quality table for %s protocol:\n', up.protocol)
for activity_no = 1 : length(activities)
    curr_activity = activities{activity_no};
    curr_activity = strrep(curr_activity, 'arm', '\textunderscore arm\textunderscore ');
    curr_activity = strrep(curr_activity, 'challenge', '\textunderscore challenge');
    curr_activity = strrep(curr_activity, 'start', '\textunderscore start');
    %param_label = find_param_label(lmeModel.Coefficients.Name{coeff_no});
    fprintf('\n\\texttt{%s} ', curr_activity);
    for metric_no = 1 : length(metrics)
        if strcmp(metrics{metric_no}, 'snr')
            fprintf('& %.1f (%.1f-%.1f) & %.1f (%.1f)', med(activity_no, metric_no), lq(activity_no, metric_no), uq(activity_no, metric_no), mu(activity_no, metric_no), sd(activity_no, metric_no));
        else
            fprintf('& %.2f (%.2f-%.2f) & %.2f (%.2f)', med(activity_no, metric_no), lq(activity_no, metric_no), uq(activity_no, metric_no), mu(activity_no, metric_no), sd(activity_no, metric_no));
        end
    end
    fprintf('  \\\\')

end
fprintf('\n')

end

function make_mixed_effects_models(data, ref_var, include_dc_amp, up)

%% Print summary of mixed effects modelling
fprintf('\n\n ~~~~~~~~ ')
fprintf('\n - Making linear mixed effects model for %s protocol using %s as a reference variable ', up.protocol, ref_var);
if include_dc_amp == 0
    fprintf('(not including DC amplitude)')
else
    fprintf('(including DC amplitude)')
end
fprintf('\n ~~~~~~~~ \n')

%% extract relevant variables

% specify variables to extract from the dataset
vars_to_extract = {ref_var, 'subj_id', 'age', 'gender', 'bmi', 'diabetes', 'fitzpatrick', 'sbp', 'dbp', 'pp', 'healthy', 'dc_amp'};

% add additional posture variables for oscillometric protocol
if strcmp(up.protocol, 'oscillometric')
    additional_oscillometric_vars = {'posture', 'arm_height'};
    vars_to_extract(end+1:end+length(additional_oscillometric_vars)) = additional_oscillometric_vars;
end

% extract table of relevant variables
tbl = data(:, vars_to_extract);

%% Specify mixed effects formula

% specify fixed effects
fixedEffectsFormula = [ref_var ' ~ 1 + age + bmi + gender + diabetes + fitzpatrick + pp + sbp'];

% add in additional posture variables for oscillometric protocol
if strcmp(up.protocol, 'oscillometric')
    for var_no = 1 : length(additional_oscillometric_vars)
        fixedEffectsFormula = [fixedEffectsFormula ' + ' additional_oscillometric_vars{var_no}];
    end
end

% add in DC amplitude if required
if include_dc_amp == 1
    fixedEffectsFormula = [fixedEffectsFormula ' + dc_amp'];
end

% add in random effect of subjects
mixedEffectsFormula = [fixedEffectsFormula ' +  (1|subj_id)'];

% Fit the linear mixed-effects model
lmeModel = fitlme(tbl, mixedEffectsFormula);

% Display the model summary
output_lme_results_latex_table(lmeModel, ref_var, up);

% assess multicolinearity
assess_multicolinearity(lmeModel, up)

% create signal quality determinants figures
if include_dc_amp
    create_signal_quality_determinants_figures(lmeModel, data, up, ref_var);
end

end

function create_ppg_quality_figure_age(data, up)

%% PPG AC and DC components

% settings
rel_subjs = [];
%tol = 0.05;
rel_dc_amp_els = data.dc_amp>2.1e4 & data.dc_amp<3.9e4;
all_rel_young_els = find(data.age<35 & data.healthy & rel_dc_amp_els);
all_rel_old_els = find(data.age>55 & data.healthy & rel_dc_amp_els);
for cat = {'young', 'old'}
    curr_cat = cat{1,1};
    eval(['curr_rel_els = all_rel_' curr_cat '_els;']);
    curr_rel_subjs = [];
    for snr_val_no = -1 : 1 : 4
        target_snr = snr_val_no*5;
        [~, temp] = min(abs(data.snr(curr_rel_els)-target_snr));
        curr_rel_subjs(end+1,1) = curr_rel_els(temp);
    end
    eval(['rel_' curr_cat '_subjs = curr_rel_subjs;']);
end

start_el = 1;
durn = 10; %nan;
    
figure('Position', [20, 20, 1000, 900])
for rel_subj_no = 1 : 2*length(rel_young_subjs)
    if rem(rel_subj_no,2) == 0
        subj_cat = 'old';
    else
        subj_cat = 'young';
    end
    eval(['curr_rel_subjs = rel_' subj_cat '_subjs;'])
    subj_no = curr_rel_subjs(ceil(rel_subj_no/2));
    
    % extract PPG signal
    ppg.fs = data.ppg(subj_no).fs;
    if ~isnan(durn)
        rel_els = start_el:start_el+(durn*ppg.fs);
        ppg.v = data.ppg(subj_no).v(rel_els);
    else
        ppg.v = data.ppg(subj_no).v;
    end
    ppg.t = [0:length(ppg.v)-1]./ppg.fs;
    ppg.snr = data.snr(subj_no);

    % filter
    [ppg_bpf.fs, ppg_hpf.fs, ppg_lpf.fs, ppg_higherpf.fs] = deal(ppg.fs);
    ppg_bpf.v = filtfilt(up.snr_bpf.b, up.snr_bpf.a, detrend(ppg.v));
    
    % Make figure
    if rem(rel_subj_no,2) == 1
        plot_options.do_ylabel = 1;
    else
        plot_options.do_ylabel = 0;
    end
    % - zoomed in original PPG
    subplot(length(rel_young_subjs),2,rel_subj_no)
    plot_options.type = 'ac_component';
    plot_options.snr_val = find_snr(ppg_bpf, up); % ppg.snr; %
    plot_options.do_title = 1;
    plot_options.title_prefix = [subj_cat ' subject'];
    plot_options.do_xlabel = 0;
    plot_options.auto_ylims = true;
    plot_ppg(ppg_bpf, plot_options)
    
    % tidy up
    box off
end

% save figure
filepath = [up.plots_path, 'PPG_quality_demo_age'];
save_figure(filepath)

end

function create_ppg_quality_figure(data, up)

%% PPG AC and DC components

% settings
rel_subjs = [];
%tol = 0.05;
rel_els = find(data.dc_amp>0e3 & data.dc_amp<50e3);
rel_els = find(data.dc_amp>2.1e4 & data.dc_amp<3.9e4);
%no_examples = 3;
target_snrs = [-5.7,-0.8,6.0];
for snr_val_no = 1:length(target_snrs)
    %target_snr = quantile(data.snr, (2*snr_val_no-1)/(2*no_examples));
    target_snr = target_snrs(snr_val_no);
    [~, temp] = min(abs(data.snr(rel_els)-target_snr));
    rel_subjs(end+1,1) = rel_els(temp);
end

start_el = 1;
durn = 10; %nan;

title_prefixes = {'a','b','c'};
    
figure('Position', [20, 20, 1200, 800])
for rel_subj_no = 1 : 3 %length(rel_subjs)
    subj_no = rel_subjs(rel_subj_no);
    % extract PPG signal
    ppg.fs = data.ppg(subj_no).fs;
    if ~isnan(durn)
        rel_els = start_el:start_el+(durn*ppg.fs);
        ppg.v = data.ppg(subj_no).v(rel_els);
    else
        ppg.v = data.ppg(subj_no).v;
    end
    ppg.t = [0:length(ppg.v)-1]./ppg.fs;
    ppg.snr = data.snr(subj_no);

    % filter
    [ppg_bpf.fs, ppg_hpf.fs, ppg_lpf.fs, ppg_higherpf.fs] = deal(ppg.fs);
    ppg_bpf.v = filtfilt(up.snr_bpf.b, up.snr_bpf.a, detrend(ppg.v));
    ppg_hpf.v = filtfilt(up.hpf.b, up.hpf.a, ppg.v); % Here high pass filter refers to above 0.5 Hz
    ppg_lpf.v = ppg.v - ppg_hpf.v; % Here LPF refers to below 0.5 Hz.
    ppg_higherpf.v = filtfilt(up.higherpf.b, up.higherpf.a, ppg.v); % here higher PF refers to above 12 Hz

    % Make figure
    if rel_subj_no == 1
        plot_options.do_ylabel = 1;
    else
        plot_options.do_ylabel = 0;
    end
    plot_options.title_prefix = title_prefixes{rel_subj_no};
    % - zoomed in original PPG
    subplot(4,length(rel_subjs),rel_subj_no)
    plot_options.type = 'zoomed_in';
    plot_options.snr_val = ppg.snr;
    plot_options.do_title = 1;
    plot_options.do_xlabel = 0;
    plot_options.auto_ylims = 1;
    plot_ppg(ppg, plot_options)
    % - DC component only
    subplot(4,length(rel_subjs),rel_subj_no+length(rel_subjs))
    plot_options.type = 'dc_component';
    plot_options.do_title = 0;
    plot_options.do_xlabel = 0;
    plot_ppg(ppg_lpf, plot_options)
    % - ac component only
    subplot(4,length(rel_subjs),rel_subj_no+2*length(rel_subjs))
    plot_options.type = 'ac_component';
    plot_options.do_title = 0;
    plot_options.do_xlabel = 0;
    plot_options.auto_ylims = 0;
    plot_ppg(ppg_bpf, plot_options)
    % - high freq noise only
    subplot(4,length(rel_subjs),rel_subj_no+3*length(rel_subjs))
    plot_options.type = 'hf_noise';
    plot_options.do_title = 0;
    plot_options.do_xlabel = 1;
    plot_ppg(ppg_higherpf, plot_options)
    
    % tidy up
    box off
end

% save figure
filepath = [up.plots_path, 'PPG_quality_demo'];
save_figure(filepath)

end


function create_exceptions_ppg_quality_figure(data, up)

%% PPG AC and DC components

% settings
rel_subjs = [];
rel_els = find(data.dc_amp>2.1e4 & data.dc_amp<3.9e4);
target_snrs = [-2,0,8.2]; % 9.5
for snr_val_no = 1:length(target_snrs)
    target_snr = target_snrs(snr_val_no);
    [~, temp] = min(abs(data.snr(rel_els)-target_snr));
    rel_subjs(end+1,1) = rel_els(temp);
end

start_el = 1;
durn = 10; %nan;

title_prefixes = {'a','b','c'};
    
figure('Position', [20, 20, 1200, 500])
for rel_subj_no = 1 : 3 %length(rel_subjs)
    subj_no = rel_subjs(rel_subj_no);
    % extract PPG signal
    ppg.fs = data.ppg(subj_no).fs;
    if ~isnan(durn)
        rel_els = start_el:start_el+(durn*ppg.fs);
        ppg.v = data.ppg(subj_no).v(rel_els);
    else
        ppg.v = data.ppg(subj_no).v;
    end
    ppg.t = [0:length(ppg.v)-1]./ppg.fs;
    ppg.snr = data.snr(subj_no);

    % filter
    [ppg_bpf.fs, ppg_hpf.fs, ppg_lpf.fs, ppg_higherpf.fs] = deal(ppg.fs);
    ppg_bpf.v = filtfilt(up.snr_bpf.b, up.snr_bpf.a, detrend(ppg.v));
    
    % Make figure
    if rel_subj_no == 1
        plot_options.do_ylabel = 1;
    else
        plot_options.do_ylabel = 0;
    end
    plot_options.title_prefix = title_prefixes{rel_subj_no};
    % - zoomed in original PPG
    subplot(2,length(rel_subjs),rel_subj_no)
    plot_options.type = 'zoomed_in';
    plot_options.snr_val = ppg.snr;
    plot_options.do_title = 1;
    plot_options.qual_type = 0;
    plot_options.do_xlabel = 0;
    plot_options.auto_ylims = 1;
    plot_ppg(ppg, plot_options)
    % - ac component only
    subplot(2,length(rel_subjs),rel_subj_no+length(rel_subjs))
    plot_options.type = 'ac_component';
    plot_options.do_title = 0;
    plot_options.do_xlabel = 0;
    plot_options.auto_ylims = 0;
    if round(ppg.snr) == 0
        plot_options.auto_ylims = 1;
    end
    plot_ppg(ppg_bpf, plot_options)
    
    % tidy up
    box off
end

% save figure
filepath = [up.plots_path, 'PPG_quality_exceptions'];
save_figure(filepath)

end

function plot_ppg(orig_ppg, plot_options)

plot_options.ftsize = 12;
plot_options.lwidth = 1;
orig_ppg.t = [0:length(orig_ppg.v)-1]./orig_ppg.fs;

switch plot_options.type
    case 'entire_signal'
        ylims = [-6e4 0];
        ylab_txt = {'Original','PPG'};
    case 'dc_component'
        %ylims = [0 4e4];
        ylims = [-4e4 0];
        ylab_txt = {'DC','component (au)'};
        %ylab_txt = {'Absolute DC','component'};
        %orig_ppg.v = abs(orig_ppg.v);
    case 'zoomed_in'
        if sum(strcmp(fieldnames(plot_options), 'auto_ylims')) && plot_options.auto_ylims
            ylims = [-inf inf];
        else
            ylim_centre = mean([min(orig_ppg.v), max(orig_ppg.v)]);
            total_range = 0.03e4;
            ylims = [ylim_centre-(total_range/2), ylim_centre+(total_range/2)];
        end
        ylab_txt = {'Original','PPG (au)'};
    case 'ac_component'
        if sum(strcmp(fieldnames(plot_options), 'auto_ylims')) && plot_options.auto_ylims
            ylims = [-inf inf];
        else
            ylims = [-400 300];
        end
        ylab_txt = {'AC','component (au)'};
    case 'hf_noise'
        ylims = [-30 30];
        ylab_txt = {'High-frequency','noise (au)'};
end

plot(orig_ppg.t, orig_ppg.v, 'b', 'LineWidth', plot_options.lwidth)
if ~strcmp(plot_options.type, 'auto')
    ylim(ylims)
end

set(gca, 'FontSize', plot_options.ftsize, 'XGrid', 'on') %, 'XTick', 0:durn)
%set(gca, 'XAxisLocation', 'origin');
if strcmp(plot_options.type, 'entire_signal') || strcmp(plot_options.type, 'dc_component') || strcmp(plot_options.type, 'zoomed_in')
    ytickformat('%.2f');
end
if plot_options.do_ylabel
    ylabel(ylab_txt,'FontSize', plot_options.ftsize)
end
if plot_options.do_xlabel
    xlabel('Time (s)', 'FontSize', plot_options.ftsize)
end
box off

if plot_options.do_title
    if plot_options.snr_val < -1
        qual_type = 'Low';
    elseif plot_options.snr_val < 5
        qual_type = 'Medium';
    else
        qual_type = 'High';
    end
    if sum(strcmp(fieldnames(plot_options), 'title_prefix'))
        title_txt = ['(', plot_options.title_prefix, ')  '];
    else
        title_txt = '';
    end
    if sum(strcmp(fieldnames(plot_options), 'qual_type')) && ~plot_options.qual_type
        title_txt = [title_txt, 'SNR = ', num2str(plot_options.snr_val, '%.1f'), ' dB'];
    else
        title_txt = [title_txt, qual_type, ' Quality, SNR = ', num2str(plot_options.snr_val, '%.1f'), ' dB'];
    end
    title({title_txt, ''},'FontSize', plot_options.ftsize+4)
end

end

function create_demo_figures(data, up)

%% PPG AC and DC components

% settings
rel_subjs = [];
%tol = 0.05;
rel_els = find(data.dc_amp>0e3 & data.dc_amp<50e3);
for snr_val_no = -1 : 5
    target_snr = snr_val_no*5;
    %snr_lower_lim = snr_val_no*5-tol;
    %snr_upper_lim = snr_val_no*5+tol;
    %rel_subjs(end+1,1) = find(data.snr>snr_lower_lim & data.snr<snr_upper_lim,1);
    [~, temp] = min(abs(data.snr(rel_els)-target_snr));
    rel_subjs(end+1,1) = rel_els(temp);
end
subj_no = 1;
ftsize = 12;
lwidth = 1;
start_el = 1;
durn = nan;


figure('Position', [20, 20, 1200, 800])
for rel_subj_no = 1 : length(rel_subjs)
    subj_no = rel_subjs(rel_subj_no);
    % extract PPG signal
    ppg.fs = data.ppg(subj_no).fs;
    if ~isnan(durn)
        rel_els = start_el:start_el+(durn*ppg.fs);
        ppg.v = data.ppg(subj_no).v(rel_els);
    else
        ppg.v = data.ppg(subj_no).v;
    end
    ppg.t = [0:length(ppg.v)-1]./ppg.fs;
    ppg.snr = data.snr(subj_no);

    % filter
    ppg.v = filtfilt(up.snr_bpf.b, up.snr_bpf.a, ppg.v);
    
    % calc snr
    % - old SNR (dumb)
    sig = filtfilt(up.pk_detect_bpf.b, up.pk_detect_bpf.a, data.ppg(subj_no).v);
    snr_old = snr(sig);
    % SNR high-freq (i.e. calc cardiac SNR after eliminating low freq resp)
    hf_noise = detrend(data.ppg(subj_no).v)-ppg.v;
    snr_new_hf = snr(ppg.v, hf_noise); % This isn't yet correct - low freqs haven't been removed
    % SNR low-freq (i.e. resp SNR)
    sig = filtfilt(up.resp_bpf.b, up.resp_bpf.a, detrend(data.ppg(subj_no).v));
    noise = detrend(data.ppg(subj_no).v)-sig;
    snr_new_lf = snr(ppg.v, noise); % this isn't yet correct

    % Make figure
    subplot(ceil(length(rel_subjs)/2),2,rel_subj_no)
    plot(ppg.t, ppg.v, 'LineWidth', lwidth), hold on,
    %plot(ppg.t, noise, 'r')

    % labels
    xlabel('Time (s)', 'FontSize', ftsize)
    ylabel('PPG (au)','FontSize', ftsize)
    title(['SNR old = ', num2str(snr_old, '%.1f') ', SNR new = ', num2str(snr_new_lf, '%.1f')],'FontSize', ftsize)
    set(gca, 'FontSize', ftsize) %, 'XTick', 0:durn)

    % tidy up
    box off
end


end

function extract_recording_chars(data, up)

fprintf('\n ~~~~~~~~ ')
fprintf('\n - Recording characteristics for %s protocol', up.protocol);
fprintf('\n ~~~~~~~~ \n')

activities = unique(data.activity);

for act_no = 1 : length(activities)
    curr_act = activities{act_no};
    rel_rows = strcmp(data.activity, curr_act);
    durns = nan(height(data),1);
    for s = 1 : height(data)
        if strcmp(data.activity(s), curr_act)
            durns(s,1) = length(data.ppg(s).v)/data.ppg(s).fs;
        end
    end
    durns = durns(~isnan(durns));
    fprintf('\n   - %s lasted %.1f (%.1f-%.1f) secs, with %d repetitions per subject, and a total of %d occurrences', curr_act, median(durns), quantile(durns,0.25), quantile(durns,0.75), sum(rel_rows)/length(unique(data.subj_id)), sum(rel_rows))
end



end

function data = import_bps(data, up)

fprintf('\n - Importing BPs')

measurements_data = importfile_measurements(up.measurements_filepath);

tol = 10; % error tolerance between two observers (mmHg)

if ~strcmp(up.protocol, 'oscillometric')
    fprintf('\n   - Only importing BPs where both SBP and DBPs were both within a tolerance of less than %d mmHg between the two observers', tol)
end

%% SBP and DBP
[data.sbp, data.dbp] = deal(nan(height(data),1));
for meas_no = 1 : height(data)
    curr_subj = data.subj_id{meas_no};
    curr_file = data.rec_id{meas_no};
    filepath = ['measurements_', up.protocol, '/', curr_subj, '/', curr_file];
    rel_measurements_data_row = strcmp(measurements_data.waveform_file_path, filepath);

    if strcmp(up.protocol, 'oscillometric') || (measurements_data.consensus_systolic_error(rel_measurements_data_row)<=tol && measurements_data.consensus_diastolic_error(rel_measurements_data_row)<=tol)
        data.sbp(meas_no) = measurements_data.sbp(rel_measurements_data_row);
        data.dbp(meas_no) = measurements_data.dbp(rel_measurements_data_row);
    end
    % leave as nan if the consensus isn't within the tolerance.
end

%% add PP
data.pp = data.sbp-data.dbp;

end

function create_signal_quality_determinants_figures(lmeModel, data, up, ref_var)

pval_to_inc = 1.01; %up.alpha; %1.01; % 
factors = lmeModel.Coefficients.Name(lmeModel.Coefficients.pValue<pval_to_inc, 1);
factors(strcmp(factors, '(Intercept)')) = [];
factors(strcmp(factors, 'arm_height_up')) = [];
factors(strcmp(factors, 'arm_height_lap')) = []; %{'arm_height'};
factors(strcmp(factors, 'posture_supine')) = []; %{'posture'};
factors(strcmp(factors, 'posture_sitting')) = []; %{'posture'};
% create a plot of effect of each factor
for factor_no = 1 : length(factors)
    curr_factor = factors{factor_no};
    % - make table with default values
    tbl_new = table();
    no_rows = 100;
    tbl_new.age = 45*ones(no_rows,1);
    tbl_new.gender = 1*ones(no_rows,1);
    tbl_new.bmi = 22.5*ones(no_rows,1);
    tbl_new.diabetes = 0*ones(no_rows,1);
    tbl_new.fitzpatrick = 3.5*ones(no_rows,1);
    tbl_new.dc_amp = 3e5*ones(no_rows,1);
    tbl_new.subj_id = repmat({'temp'}, [no_rows,1]);
    tbl_new.sbp = 120*ones(no_rows,1);
    tbl_new.dbp = 80*ones(no_rows,1);
    tbl_new.pp = 40*ones(no_rows,1); % from 120/80
    tbl_new.posture = repmat(categorical({'sitting'}), size(tbl_new,1), 1);
    tbl_new.arm_height = repmat(categorical({'lap'}), size(tbl_new,1), 1);
    tbl_new.activity = repmat({'temp'}, [no_rows,1]);
    if strcmp(up.protocol, 'oscillometric')
        rel_postures = {'supine', 'sitting', 'standing', 'sitting', 'sitting'};
        rel_hts = {'up', 'up', 'up', 'lap', 'down'};
        rel_colors = {'b', 'k', 'r', 'c', 'm'};
        rel_leg_labels = {'supine, arm up', 'sitting, arm up', 'standing, arm up', 'sitting, arm in lap', 'sitting, arm down'};
    else
        rel_postures = {'supine'};
        rel_hts = {'up'};
        rel_colors = {'b'};
        rel_leg_labels = {'supine, arm up'};
    end
    
    % make figure
    figure();
    ftsize = 20;
    
    for posture_no = 1 : length(rel_postures)
        curr_posture = rel_postures{posture_no};
        curr_ht = rel_hts{posture_no};
        tbl_new.posture = repmat(categorical({curr_posture}), size(tbl_new,1), 1);
        tbl_new.arm_height = repmat(categorical({curr_ht}), size(tbl_new,1), 1);

        % add in variable of interest
        eval(['tbl_new.' curr_factor ' = linspace(min(data.' curr_factor '), max(data.' curr_factor '), no_rows)'';']);
        % predict values
        [ypred,yCI,DF] = predict(lmeModel,tbl_new); %,'Simultaneous',true); % ie non-simultaneous CIs, see https://uk.mathworks.com/help/curvefit/confidence-and-prediction-bounds.html
        lwidth = 2;
        eval(['x_var = tbl_new.' curr_factor ';']);
        plot(x_var,ypred, 'Color', rel_colors{posture_no}, 'LineWidth', lwidth);
        hold on;
        if strcmp(up.protocol, 'auscultatory')
            plot(x_var,yCI, '-.', 'Color', rel_colors{posture_no}, 'LineWidth', 1);
        end

    end

    legend(rel_leg_labels,'NumColumns',2,'Location','north','FontSize', ftsize-4)

    set(gca, 'FontSize', ftsize)
    %xlim([-inf, max(x_var)])
    if strcmp(ref_var, 'snr')
        ylabel('SNR (dB)', 'FontSize', ftsize)
    elseif strcmp(ref_var, 'ac_dc_ratio')
        ylabel('PI (%)', 'FontSize', ftsize)
    elseif strcmp(ref_var, 'dc_amp')
        ylabel('DC amplitude (au)', 'FontSize', ftsize)
    end
    param_label = find_param_label(curr_factor);
    xlabel(strrep(param_label, '$', ''), 'FontSize', ftsize)
    ax = gca;
    ax.YGrid = 'on';
    box off
    if strcmp(ref_var, 'snr')
        ylim([0,30])
        switch curr_factor
            case 'age'
                xlim([20, 85])
            case 'fitzpatrick'
                xlim([1,6])
            case 'dc_amp'
                xlim([0, 12e5])
        end
    elseif strcmp(ref_var, 'ac_dc_ratio')
        ylim([-0.1 0.9])
    end
    % report change in SNR with this factor
    eval(['min_factor_val = min(data.' curr_factor '); max_factor_val = max(data.' curr_factor ');']);
    eval(['min_snr_val = ypred(1);'])
    eval(['max_snr_val = ypred(end);'])
    if strcmp(curr_factor, 'age')
        step_size = 10;
    elseif strcmp(curr_factor, 'fitzpatrick')
        step_size = 1;
    elseif strcmp(curr_factor, 'dc_amp')
        step_size = 2e5;
    elseif strcmp(curr_factor, 'gender')
        step_size = 1;
    elseif strcmp(curr_factor, 'pp')
        step_size = 1;
    end
    change_per_step = (max_snr_val-min_snr_val)/((max_factor_val-min_factor_val)/step_size);
    fprintf('\n   - As %s increased from %.0f to %.0f, %s went from %.2f to %.2f, which is a change of %.1f per step of %.0f', curr_factor, min_factor_val, max_factor_val, ref_var, min_snr_val, max_snr_val, change_per_step, step_size)
    % save figure
    filepath = [up.plots_path, upper(ref_var), '_', curr_factor, '_', up.protocol];
    save_figure(filepath)


end

end

function assess_multicolinearity(lmeModel, up)


fprintf('\n - Assessing multicolinearity for %s protocol:\n', up.protocol)

% identify relevant (statistically significant) variables
vars = {};
for coeff_no = 1 : height(lmeModel.Coefficients)
    nom = lmeModel.Coefficients.Name{coeff_no};
    if strcmp(nom, '(Intercept)') || strcmp(nom, 'posture_supine') || strcmp(nom, 'posture_sitting') || strcmp(nom, 'arm_height_up') || strcmp(nom, 'arm_height_lap')
        continue
    end
    if lmeModel.Coefficients.pValue(coeff_no) > up.alpha
        continue
    end
    vars{end+1} = nom;

end

% extract matrix of data
var_data = [];
for var_no = 1 : length(vars)
    curr_var = vars{var_no};
    eval(['var_data(:, var_no) = lmeModel.Variables.' curr_var ';']);
end

% exclude rows with nans
var_data = var_data(sum(isnan(var_data),2)==0,:);

% calculate correlations between variables
corrMatrix = corr(var_data);
rsquaredMatrix = corrMatrix.^2;

% output correlations
T = array2table(rsquaredMatrix);
T.Properties.VariableNames = vars;
T.Properties.RowNames = vars;
T

end

function output_lme_results_latex_table(lmeModel, ref_var, up)

fprintf('\n - LME Results:\n')
for coeff_no = 1 : height(lmeModel.Coefficients)
    nom = lmeModel.Coefficients{1,1};
    if strcmp(lmeModel.Coefficients.Name(coeff_no), '(Intercept)')
        continue
    end
    param_label = find_param_label(lmeModel.Coefficients.Name{coeff_no});
    if strcmp(param_label, 'DC amplitude (arbitrary units)')
        fprintf('\n%s & %.6f (%.6f-%.6f) & %.6f & %.2f & %.3f', param_label, lmeModel.Coefficients.Estimate(coeff_no), lmeModel.Coefficients.Lower(coeff_no), lmeModel.Coefficients.Upper(coeff_no), lmeModel.Coefficients.SE(coeff_no), lmeModel.Coefficients.tStat(coeff_no), lmeModel.Coefficients.pValue(coeff_no))
    elseif strcmp(ref_var, 'tm_cc')
        fprintf('\n%s & %.3f (%.3f-%.3f) & %.3f & %.2f & %.3f', param_label, lmeModel.Coefficients.Estimate(coeff_no), lmeModel.Coefficients.Lower(coeff_no), lmeModel.Coefficients.Upper(coeff_no), lmeModel.Coefficients.SE(coeff_no), lmeModel.Coefficients.tStat(coeff_no), lmeModel.Coefficients.pValue(coeff_no))
    else
        fprintf('\n%s & %.2f (%.2f-%.2f) & %.2f & %.2f & %.3f', param_label, lmeModel.Coefficients.Estimate(coeff_no), lmeModel.Coefficients.Lower(coeff_no), lmeModel.Coefficients.Upper(coeff_no), lmeModel.Coefficients.SE(coeff_no), lmeModel.Coefficients.tStat(coeff_no), lmeModel.Coefficients.pValue(coeff_no))
    end
    if lmeModel.Coefficients.pValue(coeff_no) < up.alpha
        fprintf('*')
    end
    fprintf('  \\\\')

end

fprintf('\n\n   - Proportion of variability explained by fitted model: %.2f (ordinary), %.2f (adjusted) R-squared\n', lmeModel.Rsquared.Ordinary, lmeModel.Rsquared.Adjusted)


end

function param_label = find_param_label(coeff_name)

switch coeff_name
    case 'age'
        param_label = 'Age (years)';
    case 'gender'
        param_label = 'Gender (F=1, M=2)';
    case 'bmi'
        param_label = 'BMI (kgm$^{-3}$)';
    case 'diabetes'
        param_label = 'Diabetes (present = 1, absent = 0)';
    case 'healthy_1'
        param_label = 'Health status (no disease = 1, disease present = 0)';
    case 'fitzpatrick'
        param_label = 'Skin colour (Fitzpatrick skin type, 1-6)';
    case 'dc_amp'
        param_label = 'DC amplitude (arbitrary units)';
    case 'dc_amp_cat'
        param_label = 'DC amplitude category';
    case 'dc_amp_from_cat_med'
        param_label = 'Deviation of DC amplitude from category median';
    case 'sbp'
        param_label = 'Systolic blood pressure (mmHg)';
    case 'dbp'
        param_label = 'Diastolic blood pressure (mmHg)';
    case 'pp'
        param_label = 'Pulse pressure (mmHg)';
    case 'posture_sitting'
        param_label = 'Sitting posture (vs. standing)';
    case 'posture_supine'
        param_label = 'Supine posture (vs. standing)';
    case 'arm_height_up'
        param_label = 'Arm up (vs. arm down)';
    case 'arm_height_lap'
        param_label = 'Hand in lap (vs. arm down)';
end

end

function [med, lq, uq, mu, sd, min_val, max_val] = print_metric(data, activity, metric_type)

eval(['rel_metric_data = data.' metric_type ';']);
rel_rows = strcmp(data.activity, activity);
fprintf('\n   - %s for %s: ', metric_type, activity)
med = median(rel_metric_data(rel_rows)); lq = quantile(rel_metric_data(rel_rows), 0.25); uq = quantile(rel_metric_data(rel_rows), 0.75);
fprintf('med (quartiles): %.2f (%.2f - %.2f); ', med, lq, uq)
mu = mean(rel_metric_data(rel_rows)); sd = std(rel_metric_data(rel_rows));
fprintf('mean (SD): %.2f (%.2f); ', mu, sd)
min_val = min(rel_metric_data(rel_rows)); max_val = max(rel_metric_data(rel_rows));
fprintf('range: %.2f - %.2f', min_val, max_val)

end

function all_p = make_plot_comparing_postures(data, sig_char, postures, rel_arm_height, all_p, up)

%% Identifying relevant postures

% identify relevant postures for arm at natural height
if strcmp(rel_arm_height, 'natural')
    postures = unique(data.posture(data.posture == 'standing' & data.arm_height=='down' | ...
        data.posture == 'sitting' & data.arm_height=='lap' | ...
        data.posture == 'supine' & data.arm_height=='up' ));
end
% re-format relevant postures
if isequal(postures, {'standing';'sitting';'supine'})
    postures = categorical({'supine';'sitting';'standing'});
end

%% Extract data

% setup variables
rel_vals = nan(2*length(unique(data.subj_id)),length(postures));
rel_vals_subj_id = cell(2*length(unique(data.subj_id)),length(postures));

% extract reference signal quality data
eval(['rel_var = data.' sig_char ';'])
rel_subj_id = data.subj_id;

% cycle through each posture
for posture_no = 1 : length(postures)
    curr_posture = postures(posture_no);

    % Identify current arm height
    if strcmp(rel_arm_height, 'natural')
        if curr_posture == 'standing'
            curr_rel_arm_height = 'down';
        elseif curr_posture == 'sitting'
            curr_rel_arm_height = 'lap';
        elseif curr_posture == 'supine'
            curr_rel_arm_height = 'up';
        end
    else
        curr_rel_arm_height = rel_arm_height;
    end

    % extract data for this posture and this arm height
    rel_rows = data.posture == curr_posture & data.arm_height == curr_rel_arm_height;
    rel_vals(1:sum(rel_rows),posture_no) = rel_var(rel_rows);
    rel_vals_subj_id(1:sum(rel_rows),posture_no) = rel_subj_id(rel_rows);
end

%% Perform comparison between postures

fprintf('\n\n ~~~~~~ ')
fprintf('\n Comparison between postures with arm at ''%s'' height:', rel_arm_height)
fprintf('\n ~~~~~~ \n')

% specify filepath to save plot comparing postures
filepath = [up.plots_path, sig_char, '_postures'];
if strcmp(rel_arm_height, 'natural')
    filepath = [filepath, '_natural'];
end

% make actual plot
plot_title = '';
all_p = make_actual_plot_comparing_postures(rel_vals, rel_vals_subj_id, plot_title, postures, filepath, sig_char, all_p, up);

end

function all_p = make_plot_at_certain_posture(data, sig_char, posture, all_p, up)

posture_cat = categorical(cellstr(posture));
arm_heights = unique(data.arm_height(data.posture==posture));

if length(arm_heights) == 1
    return
end

if length(arm_heights) == 3
    arm_heights = categorical({'down';'lap';'up'});
end
rel_vals = nan(length(unique(data.subj_id)),length(arm_heights));
rel_vals_subj_id = cell(2*length(unique(data.subj_id)),length(arm_heights));
eval(['rel_var = data.' sig_char ';'])
rel_subj_id = data.subj_id;
for height_no = 1 : length(arm_heights)
    curr_height = arm_heights(height_no);
    rel_rows = data.arm_height == curr_height & data.posture == posture_cat;
    temp = rel_var(rel_rows);
    rel_vals(1:length(temp),height_no) = temp;
    rel_vals_subj_id(1:sum(rel_rows),height_no) = rel_subj_id(rel_rows);
end

fprintf('\n ~~~~~~ ')
fprintf('\n Comparison of %s between arm heights whilst %s:', sig_char, posture);
fprintf('\n ~~~~~~ \n')

plot_title = posture;
filepath = [up.plots_path, sig_char, '_' char(posture)];
all_p = make_actual_plot_comparing_postures(rel_vals, rel_vals_subj_id, plot_title, arm_heights, filepath, sig_char, all_p, up);
end

function save_figure(filepath)

print(filepath, '-depsc');
print(filepath, '-dpng');
close all

end

function all_p = make_actual_plot_comparing_postures(rel_vals, rel_vals_subj_id, plot_title, arm_heights, filepath, plot_type, all_p, up)

fprintf('\n - Making %s plot comparing postures: %s', plot_type, plot_title)

%% Print summary of signal qualities at each posture / height
for height_no = 1 : length(arm_heights)
    fprintf('\n  - %s at %s: %.1f (%.1f - %.1f) dB', plot_type, arm_heights(height_no), nanmedian(rel_vals(:,height_no)), quantile(rel_vals(:,height_no), 0.25), quantile(rel_vals(:,height_no), 0.75))
end

%% Perform statistical tests comparing postures or arm heights

% check whether this is comparing postures or sensor heights
[extra_str, xlab_txt] = check_whether_comparing_postures(plot_title);

% Perform statistical tests
xticklabs = cell(length(arm_heights),1);
sig_bar_p = nan(length(arm_heights),length(arm_heights));

% cycle through each posture or arm height
for height_no = 1 : length(arm_heights)
    
    % identify the other postures / heights
    all_arm_heights = 1:length(arm_heights);
    other_heights = all_arm_heights(all_arm_heights>height_no);
    %other_heights = setxor(1:length(arm_heights), height_no);

    % go through each other posture / height
    for other_height_no = 1 : length(other_heights)
        curr_other_height_no = other_heights(other_height_no);

        % compare the signal quality values at the specific posture /
        % height against the current other posture / height
        vals_group1 = rel_vals(:,height_no);
        vals_group2 = rel_vals(:,curr_other_height_no);
        subjs_group1 = rel_vals_subj_id(:,height_no);
        subjs_group2 = rel_vals_subj_id(:,curr_other_height_no);
        p = compare_quality_between_postures_or_heights(vals_group1, vals_group2, subjs_group1, subjs_group2);
        
        % store result
        sig_bar_p(height_no, curr_other_height_no) = p;
        all_p(end+1,1) = p;

        % print result
        fprintf('\n  - p=%.5f for %s vs. %s', p, arm_heights(height_no), arm_heights(curr_other_height_no))
    end

end

%% Create plot

% make plot
h = boxchart(rel_vals);

% colour in
h.BoxMedianLineColor = [1,0,0];

% Labels
ftsize = 16;
if strcmp(plot_type, 'ac_amp')
    ylab_txt = 'Pulsatile (AC) amplitude (au)';
    ylim([0 800])
elseif strcmp(plot_type, 'dc_amp')
    ylab_txt = 'Baseline (DC) amplitude (au)';
    ylim([-4e5 1e5])
elseif strcmp(plot_type, 'ac_dc_ratio')
    ylab_txt = 'PI (%)';
    ylim([0 3])
elseif strcmp(plot_type, 'tm_cc')
    ylab_txt = 'TMCC';
    ylim([0.65 1.08])
elseif strcmp(plot_type, 'snr')
    ylab_txt = 'SNR (dB)';
    ylim([-10 40])
end
ylabel(ylab_txt, 'FontSize', ftsize)
xlabel(xlab_txt, 'FontSize', ftsize)
set(gca, 'FontSize', ftsize)

% create x-tick labels
for height_no = 1 : length(arm_heights)
    xticklabs{height_no} = [extra_str, char(arm_heights(height_no))];
end
xticklabels(xticklabs);

% additional settings
box off
ax = gca;
ax.YGrid = 'on';
ax.TickLabelInterpreter = 'tex';

%% add significance bars
hold on

% set up
ylims = ylim;
max_y_val = ylims(2);
y_offset = range(ylims)/10;
min_y_val = max_y_val - y_offset;

no_comparisons = sum(sum(~isnan(sig_bar_p)));
if no_comparisons > 1
    inc = (max_y_val-min_y_val)/(no_comparisons-1);
else
    inc = (max_y_val-min_y_val);
end
counter = 0;

for height_no = 1 : length(arm_heights)
    all_arm_heights = 1:length(arm_heights);
    other_heights = all_arm_heights(all_arm_heights>height_no);
    for other_height_no = 1 : length(other_heights)
        curr_other_height_no = other_heights(other_height_no);
        curr_p = sig_bar_p(height_no, curr_other_height_no);
        counter = counter+1;

        % plot sig bar
        y_val = min_y_val + (counter-1)*inc;
        plot([height_no, curr_other_height_no], [y_val, y_val], 'k');

        % annotate
        lab_txt = 'n.s.';
        if curr_p < 0.05
            lab_txt = '*';
        end
        if isnan(curr_p)
            error('this shouldnt happen')
        end
        text(mean([height_no, curr_other_height_no]), y_val+(inc/5), lab_txt, 'FontSize', ftsize, 'HorizontalAlignment', 'center');
    end
end

%% save figure
save_figure(filepath)

end

function p = compare_quality_between_postures_or_heights(vals_group1, vals_group2, subjs_group1, subjs_group2)
%
% Compares signal quality values between two groups
%
% Inputs:
%  - vals_group1 - a column vector of signal quality values for group 1
%  - vals_group2 - a column vector of signal quality values for group 2
%  - subjs_group1 - a column cell of strings, where each string is the subject ID corresponding to the value in vals_group1 (with the same dimensions as vals_group1)
%  - subjs_group2 - a column cell of strings, where each string is the subject ID corresponding to the value in vals_group2 (with the same dimensions as vals_group2)
%
% Outputs:
%  - p - p-value for comparison between groups

% Setting
stats_test = 'mixed_effects'; % one of 'wilcoxon_signed_rank' or 'mixed_effects'

% Perform comparison between groups
if strcmp(stats_test, 'wilcoxon_signed_rank')
    
    [p, ~, ~] = signrank(vals_group1, vals_group2);
    
elseif strcmp(stats_test, 'mixed_effects')

    % remove excess nans
    rel_els1 = 1: find(~isnan(vals_group1),1,'last');
    rel_els2 = 1: find(~isnan(vals_group2),1,'last');
    vals_group1 = vals_group1(rel_els1);
    vals_group2 = vals_group2(rel_els2);
    subjs_group1 = subjs_group1(rel_els1);
    subjs_group2 = subjs_group2(rel_els2);

    % Create a table for mixed-effects model
    group_labels = [repmat({'Group1'}, size(vals_group1)); repmat({'Group2'}, size(vals_group2))];
    subject_ids = [subjs_group1; subjs_group2];
    signal_quality = [vals_group1; vals_group2];

    tbl = table(signal_quality, group_labels, subject_ids, 'VariableNames', {'SNR', 'Group', 'Subject'});

    % Convert categorical variables
    tbl.Group = categorical(tbl.Group);
    tbl.Subject = categorical(tbl.Subject);

    % Fit linear mixed-effects model
    lme = fitlme(tbl, 'SNR ~ Group + (1|Subject)');

    % Extract p-value for group comparison
    coeffs = lme.Coefficients;
    p = coeffs.pValue(2); % Group comparison is the second row in coefficients

end

end


function [extra_str, xlab_txt] = check_whether_comparing_postures(plot_title)
if isempty(plot_title) || contains(plot_title, 'posture')
    extra_str = '';
    xlab_txt = {'Posture'};
else
    extra_str = 'Arm ';
    xlab_txt = {'Sensor height'};
end
end

function snr_val = find_snr(sig, up)
filtered.v = filtfilt(up.snr_bpf.b, up.snr_bpf.a, sig.v); % 16-Jan: used to use peak detection filtering, now changed code (but not re-run analysis) to snr filtering
filtered.fs = sig.fs;
snr_val = snr(filtered.v);
end

function [ac_amp, dc_amp] = find_amp_ac_dc(up, sig)

bpf_filtered.v = filtfilt(up.pk_detect_bpf.b, up.pk_detect_bpf.a, sig.v);
bpf_filtered.fs = sig.fs;
[peaks, onsets, mid_amps] = detect_ppg_beats(bpf_filtered, up.beat_detector);
ac_amp = median(sig.v(peaks)-sig.v(onsets));
dc_amp = -1*median((sig.v(mid_amps)));

end

function pt_data = importfile(filename)
%IMPORTFILE Import data from a text file
%  PARTICIPANTSDELETE = IMPORTFILE(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a table.
%
%  PARTICIPANTSDELETE = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  participantsdelete = importfile("/Users/petercharlton/Documents/Data/Aurora/raw_data/participants_delete.txt", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 20-Dec-2023 17:59:33

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 23);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["pid", "n_meas_inlab", "n_meas_ambulatory", "aurora_size", "fitzpatrick_scale", "bp_cuff_arm", "in_feature_table", "age", "height", "weight", "gender", "self_report_htn", "high_blood_pressure", "coronary_artery_disease", "diabetes", "arrythmia", "previous_heart_attack", "previous_stroke", "heart_failure", "aortic_stenosis", "valvular_heart_disease", "other_cv_diseases", "cvd_meds"];
opts.VariableTypes = ["string", "double", "double", "categorical", "double", "categorical", "double", "double", "double", "double", "categorical", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "pid", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["pid", "aurora_size", "bp_cuff_arm", "gender", "self_report_htn"], "EmptyFieldRule", "auto");

% Import the data
pt_data = readtable(filename, opts);

end

function measurementsauscultatory = importfile_measurements(filename)
%IMPORTFILE Import data from a text file
%  MEASUREMENTSAUSCULTATORY = IMPORTFILE(FILENAME) reads data from text
%  file FILENAME for the default selection.  Returns the data as a table.
%
%  MEASUREMENTSAUSCULTATORY = IMPORTFILE(FILE, DATALINES) reads data for
%  the specified row interval(s) of text file FILENAME. Specify
%  DATALINES as a positive scalar integer or a N-by-2 array of positive
%  scalar integers for dis-contiguous row intervals.
%
%  Example:
%  measurementsauscultatory = importfile("/Users/petercharlton/Documents/Data/Aurora/raw_data/measurements_auscultatory.txt", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 04-Jan-2024 11:47:29

dataLines = [2, Inf];

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 17);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["pid", "phase", "measurement", "date_time", "sbp", "dbp", "duration", "pressure_quality", "optical_quality", "waveform_file_path", "waveforms_generated", "primary_systolic", "primary_diastolic", "secondary_systolic", "secondary_diastolic", "consensus_systolic_error", "consensus_diastolic_error"];
opts.VariableTypes = ["double", "categorical", "categorical", "datetime", "double", "double", "double", "double", "double", "string", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "waveform_file_path", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["phase", "measurement", "waveform_file_path"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "date_time", "InputFormat", "yyyy-MM-dd HH:mm:ss");
opts = setvaropts(opts, "pid", "TrimNonNumeric", true);
opts = setvaropts(opts, "pid", "ThousandsSeparator", ",");

% Import the data
measurementsauscultatory = readtable(filename, opts);

end

function extract_subj_chars(data, up)

fprintf('\n\n ~~~~~~~~ ')
fprintf('\n - Subject characteristics for %s protocol:', up.protocol)
fprintf('\n ~~~~~~~~ \n')

[rel_subjs, rel_inds] = unique(data.subj_id);
data = data(rel_inds,:);
no_pts = height(data);
fprintf('\nNo. subjects & %d \\\\', no_pts);
fprintf('\nFemale, n (\\%%) & %d (%.1f) \\\\', sum(data.gender==1), 100*sum(data.gender==1)/no_pts);
fprintf('\nAge (years), mean (SD) & %.1f (%.1f) \\\\', mean(data.age), std(data.age));
bmi = (0.453592*data.weight)./((data.ht*2.54/100).^2);
fprintf('\nBMI (kgm$^{-2}$), mean (SD) & %.1f (%.1f) \\\\', mean(bmi(~isnan(bmi))), std(bmi(~isnan(bmi))));
fprintf('\nSystolic blood pressure (mmHg), mean (SD) & %.0f (%.0f) \\\\', mean(data.sbp(~isnan(data.sbp))), std(data.sbp(~isnan(data.sbp))));
fprintf('\nDiastolic blood pressure (mmHg), mean (SD) & %.0f (%.0f) \\\\', mean(data.dbp(~isnan(data.dbp))), std(data.dbp(~isnan(data.dbp))));
fprintf('\nSkin colour (Fitzpatrick scale), n (\\%%): & \\\\')
colors = unique(data.fitzpatrick(~isnan(data.fitzpatrick))); 
colors = [colors; nan];
for col_no = 1 : length(colors)
    curr_col = colors(col_no);
    if ~isnan(curr_col)
        n = sum(data.fitzpatrick == curr_col);
        fprintf('\n - %d & %d (%.1f) \\\\', curr_col, n, 100*n/no_pts);
    else
        n = sum(isnan(data.fitzpatrick));
        fprintf('\n - unknown & %d (%.1f) \\\\', n, 100*n/no_pts);
    end
end
fprintf('\nSelf-reported history of, n (\\%%): & \\\\')
diseases = up.diseases; 
for dis_no = 1 : length(diseases)
    curr_dis = diseases{dis_no};
    eval(['dis_data = data.' curr_dis ';'])
    fprintf('\n - %s & %d (%.1f) \\\\', strrep(strrep(curr_dis, '_', ' '), 'cv diseases', 'cardiovascular diseases'), sum(dis_data), 100*sum(dis_data)/no_pts);
end
fprintf('\n - currently taking medicine for cardiovascular disease & %d (%.1f) \\\\', sum(data.cvd_meds), 100*sum(data.cvd_meds)/no_pts);
fprintf('\n - No self-reported cardiovascular disease or medication & %d (%.1f) \\\\', sum(data.healthy), 100*sum(data.healthy)/no_pts);
fprintf('\nDC amplitude (arbitrary units), n (\\%%): & \\\\')
dc_amps = [0, 50000; 295000, 330000; 570000, 605000; 850000, 900000; 1135000, 1140000]; 
for amp_no = 1 : size(dc_amps,1)
    n = sum(data.dc_amp > dc_amps(amp_no,1) & data.dc_amp < dc_amps(amp_no,2));
    fprintf('\n - %.0fk-%.0fk & %d (%.1f) \\\\', dc_amps(amp_no,1)/1000, dc_amps(amp_no,2)/1000, n, 100*n/no_pts);
end

fprintf('\n')

end

function up = setup_up(curr_protocol)
% Setup universal parameters

fprintf('\n - Setting up general universal params')

close all

%% general settings
up.protocols = {'auscultatory', 'oscillometric'};
up.signal_quality_metrics = {'snr','tm_cc','ac_dc_ratio'};
up.redo_extraction = false; % whether to re-do the extraction of data from the dataset.
up.beat_detector = 'MSPTD';
up.diseases = {'managed_hypertension', 'unmanaged_hypertension' , 'high_blood_pressure','diabetes','arrythmia','previous_stroke','previous_heart_attack','coronary_artery_disease','heart_failure','aortic_stenosis','valvular_heart_disease','other_cv_diseases'};
up.alpha = 0.05;  % significance level

% paths
up.root_data_path = '/Users/petercharlton/Documents/Data/Aurora/';
up.plots_path = '/Users/petercharlton/Library/CloudStorage/GoogleDrive-peterhcharlton@gmail.com/My Drive/Work/Publications/In Preparation/2024 Determinants wrist PPG signal quality/PLOS DH reformatting 202505/source_files/figures_reproduction/';

%% Design filters
fprintf('\n - Designing filters')
Fs = 500;
% Design a high-pass Butterworth filter
up.hpf.order = 4; % Choose the filter order (adjust as needed)
cutoff_frequency = 0.5; % Cutoff frequency in Hz
[up.hpf.b, up.hpf.a] = butter(up.hpf.order, cutoff_frequency/(Fs/2), 'high');
% Design a higher-pass Butterworth filter
up.lpf.order = 4; % Choose the filter order (adjust as needed)
cutoff_frequency = 12; % Cutoff frequency in Hz
[up.higherpf.b, up.higherpf.a] = butter(up.lpf.order, cutoff_frequency/(Fs/2), 'high');
% Design a band-pass Butterworth filter
up.pk_detect_bpf.order = 4;
Nyquist = Fs / 2; % Nyquist frequency
low_freq = 0.5; % Lower cutoff frequency in Hz
high_freq = 8; % Upper cutoff frequency in Hz
Wn = [low_freq, high_freq] / Nyquist; % Normalize the frequencies by the Nyquist frequency
[up.pk_detect_bpf.b, up.pk_detect_bpf.a] = butter(up.pk_detect_bpf.order, Wn, 'bandpass'); % Design the Butterworth bandpass filter
% Design a band-pass Chebyshev filter
up.snr_bpf.order = 4;
Nyquist = Fs / 2; % Nyquist frequency
low_freq = 0.5; % Lower cutoff frequency in Hz
high_freq = 12; % Upper cutoff frequency in Hz
Wn = [low_freq, high_freq] / Nyquist; % Normalize the frequencies by the Nyquist frequency
[up.snr_bpf.b, up.snr_bpf.a] = cheby2(up.snr_bpf.order, 20, Wn); % Design the Chebyshev II bandpass filter
% Design a band-pass Butterworth filter (resp)
up.resp_bpf.order = 4;
Nyquist = Fs / 2; % Nyquist frequency
low_freq = (4/60); % Lower cutoff frequency in Hz
high_freq = (45/60); % Upper cutoff frequency in Hz
Wn = [low_freq, high_freq] / Nyquist; % Normalize the frequencies by the Nyquist frequency
[up.resp_bpf.b, up.resp_bpf.a] = butter(up.resp_bpf.order, Wn, 'bandpass'); % Design the Butterworth bandpass filter

% stop here if no protocol specified
if nargin<1
    return
end

%% protocol-specific settings
up.protocol = curr_protocol;
fprintf('\n - Setting up universal params for %s protocol', up.protocol)

% paths
up.conv_data_folder = [up.root_data_path, 'conv_data_reproduction/'];
up.conv_data_filepath = [up.conv_data_folder, curr_protocol, filesep, 'data.mat'];
up.measurements_filepath = [up.root_data_path, 'raw_data/measurements_', curr_protocol, '.txt'];

end

function all_data = load_data(up)

fprintf('\n - Loading extracted data (previously extracted using collate_aurora_data.m)')
all_data = table;

if strcmp(up.protocol, 'oscillometric')
    
    root_folder = [up.conv_data_folder, '/oscillometric/'];

    % - standing
    filename = [root_folder, 'aurora_standingarmdown_data.mat'];
    posture = 'standing';
    arm_height = 'down';
    all_data = load_specific_data(all_data,filename,posture,arm_height);
    
    filename = [root_folder, 'aurora_standingarmup_data.mat'];
    posture = 'standing';
    arm_height = 'up';
    all_data = load_specific_data(all_data,filename,posture,arm_height);

    % - sitting
    filename = [root_folder, 'aurora_sittingarmdown_data.mat'];
    posture = 'sitting';
    arm_height = 'down';
    all_data = load_specific_data(all_data,filename,posture,arm_height);

    filename = [root_folder, 'aurora_sittingarmlap_data.mat'];
    posture = 'sitting';
    arm_height = 'lap';
    all_data = load_specific_data(all_data,filename,posture,arm_height);

    filename = [root_folder, 'aurora_sittingarmup_data.mat'];
    posture = 'sitting';
    arm_height = 'up';
    all_data = load_specific_data(all_data,filename,posture,arm_height);

    % - posture
    filename = [root_folder, 'aurora_supine_data.mat'];
    posture = 'supine';
    arm_height = 'up';
    all_data = load_specific_data(all_data,filename,posture,arm_height);


elseif strcmp(up.protocol, 'auscultatory')
    
    root_folder = [up.conv_data_folder, '/auscultatory/'];

    % - supine
    filename = [root_folder, 'aurora_staticchallengestart_data.mat'];
    posture = 'supine';
    arm_height = 'up';
    all_data = load_specific_data(all_data,filename,posture,arm_height);

    % - supine
    filename = [root_folder, 'aurora_calibrationstart_data.mat'];
    posture = 'supine';
    arm_height = 'up';
    all_data = load_specific_data(all_data,filename,posture,arm_height);

end

end

function all_data = load_specific_data(all_data,filename,posture,arm_height)

fprintf('\n   - Loading data for %s and arm %s', posture, arm_height);
load(filename);
data = struct2table(data);
posture_arm_height = categorical(cellstr(repmat([posture, arm_height], [height(data),1])));
data = addvars(data, posture_arm_height);
new_vars = {'posture', 'arm_height'};
for var_no = 1 : length(new_vars)
    curr_var = new_vars{var_no};
    eval([curr_var ' = categorical(cellstr(repmat(' curr_var ', [height(data),1])));']);
    eval(['data = addvars(data,' curr_var ');']);
end
all_data = [all_data; data];
end

function data = identify_subjs_for_inclusion(data, up)

fprintf('\n - Identifying subjects for inclusion:\n    - out of %d original subjs,', length(unique(data.subj_id)))

if strcmp(up.protocol, 'oscillometric')
    % need to have one recording for each data subset, except supine which should have two recordings

    % determine whether or not to include each subject
    unique_subj_ids = unique(data.subj_id);
    unique_posture_arm_heights = unique(data.posture_arm_height);
    for id_no = 1 : length(unique_subj_ids)
        curr_id = unique_subj_ids{id_no};
        curr_rows = strcmp(data.subj_id, curr_id);
        curr_posture_arm_heights = data.posture_arm_height(curr_rows);
        if length(unique(curr_posture_arm_heights)) ~= length(unique_posture_arm_heights) || ...
                sum(curr_posture_arm_heights == 'supineup') ~=2
            fprintf('\n   - Removing subj %s', curr_id)
            data(curr_rows,:) = [];
            continue
        end
        %     if id_no > 30 %%% CHANGE
        %         data(curr_rows,:) = [];
        %     end

    end

elseif strcmp(up.protocol, 'auscultatory')

    % need to have three static supine challenges (i.e. activity of 'staticchallengestart') for a subject to be included
    unique_subj_ids = unique(data.subj_id);
    no_expected_recordings = 6;
    for id_no = 1 : length(unique_subj_ids)
        curr_id = unique_subj_ids{id_no};
        curr_rows = strcmp(data.subj_id, curr_id);
        if sum(curr_rows) ~= no_expected_recordings
            fprintf('\n   - Removing subj %s as it had %d out of %d recordings', curr_id, sum(curr_rows), no_expected_recordings)
            data(curr_rows,:) = [];
            continue
        end
    end

end

fprintf('\n    - %d subjects retained as they had the expected number of recordings,', length(unique(data.subj_id)))

end

function data = extract_sig_chars(data, up)

fprintf('\n - Extracting signal quality metrics:')

% for each metric
tic
for row_no = 1 : height(data)
    
    if rem(row_no,100) == 0
        txt = ['Estimated time remaining: ' num2str((height(data)-row_no)*toc/row_no) ' seconds'];
        if row_no<200
            f = waitbar(row_no/height(data), txt);
        else
            waitbar(row_no/height(data), f, txt);
        end
    end
    
    curr_ppg = data.ppg(row_no);

    % extract PPG signal quality metrics using ppg-quality toolbox
    % new approach (using ppg-quality toolbox)
    options.quality_metrics = {'snr', 'amp_metrics', 'tm_cc'};
    options.win_durn = (length(curr_ppg.v)-1)/curr_ppg.fs;
    qual = assess_ppg_quality(curr_ppg.v, curr_ppg.fs, options);
    qual.snr = mean(qual.snr(~isnan(qual.snr)));    % the old script didn't provide beat-by-beat snr, so just chose the mean
    qual.tm_cc = mean(qual.tm_cc(~isnan(qual.tm_cc)));    % the mean was used in the old script
    qual.ac_amp = median(qual.ac_amp(~isnan(qual.ac_amp)));   % the median was used in the old script
    qual.dc_amp = median(qual.dc_amp(~isnan(qual.dc_amp)));
    qual.ac_dc_ratio = 100*abs(qual.ac_amp)/abs(qual.dc_amp);
    % old approach
    %qual_old = assess_ppg_quality_old(curr_ppg.v, curr_ppg.fs);
    
    % store PPG signal quality metrics
    fields = fieldnames(qual);
    for field_no = 1 : length(fields)
        curr_field = fields{field_no};
        eval(['data.' curr_field '(row_no) = qual.' curr_field ';'])
    end

end

% Use positive DC amplitude
data.dc_amp = abs(data.dc_amp);

end
