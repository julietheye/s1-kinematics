function [currentMuscleTable,currentNeurEval,columnNames] = kinModelAllMusc(trial_data_cell,td_trim,params)
    
    %default params
    arrayname = 'S1';
    model_aliases = {'allMuscEMG'};
    models_to_plot = {'S1_FR','allMuscEMG'};
    postBumpWindow = false;
    bumpDuration = 0.125/0.005; %not in this TD, but can get it from CDS
    %found in actPasAnalysis JH
    endofMove = 18; %window (bins) between end of bump and end of movement 
    %removed 'lat_dorsi_sup','lat_dorsi_inf' bec only one lat emg muscle
    muscleArray = {'bicep_lh','bicep_sh','brachialis',...
        'brachioradialis','deltoid_ant','deltoid_med','deltoid_pos',...
        'ext_carp_rad_brevis','ext_carpi_ulnaris',...
        'ext_digitorum','flex_carpi_radialis','flex_carpi_ulnaris',...
        'flex_digit_superficialis','infraspinatus',...
        'lat_dorsi_cen','pectoralis_sup',...
        'pectoralis_inf',...
        'teres_major','tricep_lat','tricep_lon','tricep_sho'};
    %removed 'trap' bec no equivalent muscle kinematics
    notMusc = [];
    assignParams(who,params); %overwrite params
    
    %drop muscle to get to 20
    muscleArray(notMusc) = [];
    
    %get muscle indeces
    [muscleIndeces] = find(ismember(td_trim(1).muscle_names,muscleArray));

    neuron_eval_cell =cell(length(trial_data_cell),1);
    fileclock = tic;
    fprintf('Starting analysis for %d files. This will take a while...',length(trial_data_cell))
    for filenum = 1:length(trial_data_cell)
        td = trial_data_cell{filenum};

        if ~postBumpWindow
            % trim to just movements
            td = trimTD(td,{'idx_movement_on',0},{'idx_movement_on',11});
            [~,td_act] = getTDidx(td,'ctrHoldBump',false);
            [~,td_pas] = getTDidx(td,'ctrHoldBump',true);
        end
        
        if postBumpWindow
            % trim to just post-bump movements for passive movements 
            % (and equivalent time window for active reaches)
            [~,td_act] = getTDidx(td,'ctrHoldBump',false);
            [~,td_pas] = getTDidx(td,'ctrHoldBump',true);
            td_pas = trimTD(td_pas,{'idx_bumpTime',bumpDuration+1},{'idx_bumpTime',bumpDuration+1+endofMove-1});
            td_act = trimTD(td_act,{'idx_endTime',-1*numel(td_pas(1).pos(:,1))+1},{'idx_endTime',0});
            td = cat(2,td_pas,td_act);
        end

        % check to make sure all neurons fire at least once in each condition (pretty rare that one doesn't)
        firing_units = mean(getSig(td_act,'S1_spikes'))~=0 & mean(getSig(td_pas,'S1_spikes'))~=0;
        if any(~firing_units)
            unit_ids = td(1).([arrayname '_unit_guide']);
            new_unit_guide = unit_ids(firing_units,:);

            for trialnum = 1:length(td)
                td(trialnum).(sprintf('%s_unit_guide',arrayname)) = new_unit_guide;

                spikes = td(trialnum).(sprintf('%s_spikes',arrayname));
                spikes(:,~firing_units) = [];
                td(trialnum).(sprintf('%s_spikes',arrayname)) = spikes;
            end
            fprintf('Removed %d neurons for not firing in one condition\n',sum(~firing_units))
        end

        % add firing rates in addition to spike counts
        td = addFiringRates(td,struct('array',arrayname));

        % find average over the movement
        td = binTD(td,'average');
        
        %find L0 and normalize for all muscles (1 at a time)
        for muscleNum=1:length(muscleArray)
            
            %find l-naught for muscle length/velocity normalization
            trialsL0 = zeros(1,numel(td_trim));
            muscleIndex = find(ismember(td_trim(1).muscle_names,muscleArray(muscleNum)));
            for t=1:numel(td_trim)
                %for bump trials
                if td_trim(t).ctrHoldBump==1
                    trialL0 = mean(td_trim(t).muscle_len(td_trim(t).idx_startTime:td_trim(t).idx_bumpTime,muscleIndex));
                elseif td_trim(t).ctrHoldBump==0
                    trialL0 = mean(td_trim(t).muscle_len(td_trim(t).idx_startTime:td_trim(t).idx_goCueTime,muscleIndex));
                end
                trialsL0(1,t) = trialL0;
            end
            L0 = mean(trialsL0);
            
            for tr = 1:numel(td)
                %normalize muscle length and velocity
                td(tr).muscle_len(:,muscleIndex) = td(tr).muscle_len(:,muscleIndex)./L0;
                td(tr).muscle_vel(:,muscleIndex) = td(tr).muscle_vel(:,muscleIndex)./L0;

                %nonlinearize if muscNonlin model
                if strcmpi(model_aliases,'muscNonlin')
                    if curExp<1
                        td(tr).muscle_vel(:,muscleIndex) = sign(td(tr).muscle_vel(:,muscleIndex)).*(abs(td(tr).muscle_vel(:,muscleIndex))).^curExp;
                    elseif curExp>1
                        td(tr).muscle_vel(:,muscleIndex) = (td(tr).muscle_vel(:,muscleIndex)).^curExp;
                    end
                end
            end
        end

        %% find separabilities
        % suppress getTDfields warning...
        getTDfields(td,'time');
        onetime_warn = warning('query','last'); 
        warning('off',onetime_warn.identifier)

        sepResults = actpasSep(td,struct(...
            'neural_signals',[arrayname '_FR'],...
            'model_aliases',{model_aliases},...
            'muscleIndeces',muscleIndeces));

        % turn warning back on
        warning('on',onetime_warn.identifier)

        % extract neuron_eval_table and trial_table
        % replace infs with nans
        numeric_cols = strcmpi(sepResults.neuron_eval_table.Properties.VariableDescriptions,'linear');
        numeric_vals = sepResults.neuron_eval_table(:,numeric_cols).Variables;
        infidx = isinf(numeric_vals);
        numeric_vals(infidx) = NaN;
        sepResults.neuron_eval_table(:,numeric_cols).Variables = numeric_vals;

        % compile neuron eval table together
        neuron_eval_cell{filenum} = sepResults.neuron_eval_table;

        % extract only the columns we want to keep
        neuron_eval_cell{filenum}.Properties.VariableNames = strrep(neuron_eval_cell{filenum}.Properties.VariableNames,'glm_','');
        neuron_eval_cell{filenum}.Properties.VariableNames = strrep(neuron_eval_cell{filenum}.Properties.VariableNames,'model_','');
        cols_to_keep = [...
            {'monkey','date','task','signalID','crossvalID'},...
            strcat(models_to_plot(2:end),'_eval'),...
            strcat(models_to_plot(2:end),'_act_eval'),...
            strcat(models_to_plot(2:end),'_pas_eval'),...
            strcat(models_to_plot(2:end),'_train_act_eval'),...
            strcat(models_to_plot(2:end),'_train_pas_eval'),...
            strcat(models_to_plot(2:end),'_half_full_train_act_eval'),...
            strcat(models_to_plot(2:end),'_half_full_train_pas_eval'),...
            strcat(models_to_plot(1),'_indiv_sep')];

        neuron_eval_cell{filenum} = neuron_eval_cell{filenum}(:,cols_to_keep);

        % output a counter
        fprintf('Processed file %d of %d at time %f\n',filenum,length(trial_data_cell),toc(fileclock))
    end

    % compile and average
    neuron_eval = vertcat(neuron_eval_cell{:});
    avg_neuron_eval = neuronAverage(neuron_eval,struct(...
        'keycols',{{'monkey','date','task','signalID'}},...
        'do_ci',false,...
        'do_nanmean',true));
    
    % save each muscle neuron_eval into 3d array
    currentNeurEval = table2array(neuron_eval(:,6:end));

    % save each muscle avg_neuron_eval into 3d array
    currentMuscleTable = table2array(avg_neuron_eval(:,5:end));
    
    temp = avg_neuron_eval;
    temp.monkey=[]; temp.date=[]; temp.task=[]; temp.signalID=[];
    columnNames = temp.Properties.VariableNames;
    