function collate_aurora_data
%
% COLLATE_AURORA_DATA
%
% Code for curating the Aurora-BP Dataset in preparation for the analysis
% reported in:
%    Charlton PH et al., 'Determinants of photoplethysmography signal
%    quality at the wrist'.
%

fprintf('\n ~~~ Extracting Aurora data ~~~')

% Setup
protocols = {'oscillometric', 'auscultatory'};
root_folder = '/Users/petercharlton/Documents/Data/Aurora/';
max_subjs = inf; % set to inf to do whole dataset
sigs = {'ekg','ppg','acc_ppg_site'};
up.pt_data_filename = '/Users/petercharlton/Documents/Data/Aurora/raw_data/participants_delete.txt';
up.diseases = {'self_report_htn', 'high_blood_pressure','diabetes','arrythmia','previous_stroke','previous_heart_attack','coronary_artery_disease','heart_failure','aortic_stenosis','valvular_heart_disease','other_cv_diseases', 'cvd_meds'};

for protocol_no = 1 : length(protocols)
    curr_protocol = protocols{protocol_no};

    fprintf('\n Extracting data for %s protocol', curr_protocol);

    % setup paths
    raw_data_folder = [root_folder, 'raw_data', filesep, 'measurements_', curr_protocol, filesep];
    save_folder = [root_folder, 'conv_data_reproduction', filesep, curr_protocol, filesep];

    % identify subjects
    subjs = identify_subjs(raw_data_folder);

    % load data and store in single structure
    data = load_data(raw_data_folder, subjs, max_subjs);

    % remove any data with flat signals
    data = remove_flat_signals(data, sigs);

    % insert demographics
    data = insert_demographics(data, up);

    % save data for each activity
    save_data(data, save_folder);

end

fprintf('\n ~~~ Finished ~~~')

end

function data = insert_demographics(data, up)

fprintf('\n Extracting subject characteristics:\n')

% import demographic data from file
pt_data = importfile(up.pt_data_filename);

% insert demographic data into data structure
subj_data = extractfield(data, 'subj_id');
rel_subjs = unique(subj_data);
demogs.orig = {'fitzpatrick_scale', 'age', 'height', 'weight', 'gender'};
demogs.orig = [demogs.orig, up.diseases];
demogs.new = {'fitzpatrick', 'age', 'ht', 'weight', 'gender'};
demogs.new = [demogs.new, up.diseases];
for subj_no = 1 : length(rel_subjs)
    curr_subj = rel_subjs{subj_no};
    subj_rows = find(strcmp(subj_data, curr_subj));
    rel_pt_data_row = strcmp(pt_data.pid, curr_subj);
    curr_subj_demogs = pt_data(rel_pt_data_row,:);
    for demog_no = 1 : length(demogs.orig)
        for subj_row_no = 1 : length(subj_rows)
            eval(['data(subj_rows(subj_row_no)).' demogs.new{demog_no} ' = curr_subj_demogs.' demogs.orig{demog_no} ';']);
        end
    end

end

% convert variable types
for row_no = 1 : length(data)
    % gender
    if data(row_no).gender == 'F'
        data(row_no).gender = 1;
    elseif data(row_no).gender == 'M'
        data(row_no).gender = 2;
    end
    % self_report_htn
    if data(row_no).self_report_htn == 'managed'
        data(row_no).managed_hypertension = 1;
    else
        data(row_no).managed_hypertension = 0;
    end
    if data(row_no).self_report_htn == 'unmanaged'
        data(row_no).unmanaged_hypertension = 1;
    else
        data(row_no).unmanaged_hypertension = 0;
    end
end

% changing diseases
up.diseases = up.diseases(~strcmp(up.diseases, 'self_report_htn'));
up.diseases = ['managed_hypertension', 'unmanaged_hypertension' , up.diseases];
% adding derived variables
for row_no = 1 : length(data)
    % insert bmi
    data(row_no).bmi = (0.453592*data(row_no).weight)./((data(row_no).ht*2.54/100).^2);
    % insert whether healthy or not
    healthy = true;
    for dis_no = 1 : length(up.diseases)
        eval(['curr_dis_data = data(row_no).' up.diseases{dis_no} ';']);
        healthy = healthy & ~curr_dis_data;
    end
    data(row_no).healthy = healthy;
end

end

function subjs = identify_subjs(raw_data_folder)
subjs = dir(raw_data_folder);
subjs = extractfield(subjs, 'name');
subjs = subjs(~strcmp(subjs, '.') & ~strcmp(subjs, '..') & ~strcmp(subjs, '.DS_Store'));
end

function filedata = load_file(curr_path)

opts = delimitedTextImportOptions("NumVariables", 7);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["t", "ekg", "optical", "pressure", "accel_x", "accel_y", "accel_z"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
filedata = readtable(curr_path, opts);

end

function filedata = reformat_data(filedata, counter_no)

fs = round(1/ ((filedata.t(end)-filedata.t(1))/(length(filedata.t)-1)));
temp.ekg.v = filedata.ekg;
temp.ekg.fs = fs;
temp.ppg.v = filedata.optical;
temp.ppg.fs = fs;
if counter_no == 0
    warning('assuming that PPG and accel are measured at the same site')
end
%    - convert to milligravitational units
filedata.accel_x = 1000*filedata.accel_x;
filedata.accel_y = 1000*filedata.accel_y;
filedata.accel_z = 1000*filedata.accel_z;
temp.acc_ppg_site.v = sqrt(filedata.accel_x.^2 + filedata.accel_y.^2 + filedata.accel_z.^2);
temp.acc_ppg_site.fs = fs;

filedata = temp;

end

function data = load_data(raw_data_folder, subjs, max_subjs)
    
data = struct;
counter_no = 0;
subj_counter = 0;
fprintf('\n - Loading data from %d out of %d subjects:', min([length(subjs), max_subjs]), length(subjs))
for subj_no = 1 : length(subjs)
    if subj_counter >= max_subjs
        continue
    else
        subj_counter = subj_counter+1;
    end
    fprintf('\n   - Subject %d', subj_no)
    subj_folder = [raw_data_folder, subjs{subj_no}, filesep];
    subj_files = dir([subj_folder, '*.tsv']);
    if isempty(subj_files)
        continue
    end
    subj_files = extractfield(subj_files, 'name');
    % cycle through each file for this subject
    for file_no = 1 : length(subj_files)
        fprintf(', %d', file_no)
        curr_filename = subj_files{file_no};
        % skip if this was measured in the return session (so that everyone contributes a maximum of one recording for each posture / sensor height)
        if contains(curr_filename, 'return')
            fprintf(' (skipped as return session)')
            continue
        end
        % identify activity of this file
        activity = curr_filename(1:end-4); % gets rid of extension
        temp = strfind(activity, '.');
        activity = activity(temp(2)+1:end);
        activity = strrep(activity, 'Cool_down_1', 'Cool_down_one');
        activity = strrep(activity, 'Cool_down_2', 'Cool_down_two');
        activity = regexprep(activity, '[0-9_]', ''); % get rid of numbers and underscores
        activity = lower(strrep(activity, 'measurement', 'ambulatory'));
        % load data from this file
        curr_path = [subj_folder, curr_filename];
        filedata = load_file(curr_path);
        % reformat data
        filedata = reformat_data(filedata, counter_no);
        % store data
        counter_no = counter_no+1;
        filedata.subj_id = subjs{subj_no};
        filedata.rec_id = subj_files{file_no};
        filedata.activity = activity;
        if subj_no == 1 & file_no == 1
            data = filedata;
        else
            data(counter_no) = filedata;
        end
        clear filedata

    end

end

end

function data = remove_flat_signals(data, sigs)
rows_to_exc = false(length(data),1);
subj_ids = extractfield(data, 'subj_id');
for row_no = 1 : length(data)
    for sig_no = 1 : length(sigs)
        curr_sig = sigs{sig_no};
        eval(['curr_sig_data = data(row_no).' curr_sig ';']);
        if length(unique(curr_sig_data.v)) == 1
            curr_subj = data(row_no).subj_id;
            curr_row_nos_to_exc = strcmp(subj_ids, curr_subj);
            rows_to_exc(curr_row_nos_to_exc) = true;
        end
    end
end
fprintf('\n - Removing recordings for whom at least one of the following was a flat line: ')
fprintf(1,'%s, ', sigs{:})
fprintf('\n   - Originally %d recordings from %d subjects', length(data), length(unique(extractfield(data, 'subj_id'))))
data = data(~rows_to_exc);
fprintf('\n   - Removed data from %d recordings', sum(rows_to_exc))
fprintf('\n   - Leaving %d recordings from %d subjects', length(data), length(unique(extractfield(data, 'subj_id'))))

end

function save_data(data, save_folder)
all_data = data; clear data
activities = extractfield(all_data, 'activity');
unique_activities = unique(activities);
for activity_no = 1: length(unique_activities)
    curr_activity = unique_activities{activity_no};
    fprintf('\n - Saving data for %s:', curr_activity)
    % identify data for this activity
    rel_rows = strcmp(activities, curr_activity);
    data = all_data(rel_rows);
    fprintf(' %d recordings from %d subjects', sum(rel_rows), length(unique(extractfield(data, 'subj_id'))))
    % check save folder exists
    if ~exist(save_folder,'dir')
        mkdir(save_folder)
    end
    % save data for this activity
    savename = ['aurora_' curr_activity '_data'];
    save([save_folder,savename], 'data')
end

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