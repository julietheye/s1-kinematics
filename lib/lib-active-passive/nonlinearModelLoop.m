function [multiExp4dArray,multiExp4dNeurEval,columnNames] = nonlinearModelLoop(trial_data_cell,td_trim,params)

    %default params
    arrayname = 'S1';
    model_aliases = {'muscNonlin'};
    models_to_plot = {'S1_FR','muscNonlin'};
    muscleArray = {'bicep_lh','bicep_sh','brachialis',...
            'brachioradialis','deltoid_ant','deltoid_med','deltoid_pos',...
            'ext_carp_rad_brevis','ext_carpi_ulnaris',...
            'ext_digitorum','flex_carpi_radialis','flex_carpi_ulnaris',...
            'flex_digit_superficialis','infraspinatus',...
            'lat_dorsi_sup','lat_dorsi_cen','lat_dorsi_inf','pectoralis_sup',...
            'pectoralis_inf',...
            'teres_major','tricep_lat','tricep_lon','tricep_sho'};
    exponents = [1];
    assignParams(who,params); %overwrite params
    
    td = trial_data_cell{1};
    multiExp4dArray = zeros(size(td(1).S1_spikes,2),8,numel(muscleArray),numel(exponents));
    multiExp4dNeurEval = zeros(size(td(1).S1_spikes,2)*100,8,numel(muscleArray),numel(exponents));
    for d=1:numel(exponents)
        curExp = exponents(d);

        allMuscAvgNeurEval = [];
        allMuscNeurEval = [];

        for muscleNum=1:length(muscleArray)

            fprintf('exp: %i',exponents(d))
            fprintf('muscle: %i',muscleNum)

            muscle = muscleArray(muscleNum);

            neuron_eval_cell =cell(length(trial_data_cell),1);
            fileclock = tic;
            fprintf('Starting analysis for %d files. This will take a while...',length(trial_data_cell))
            for filenum = 1:length(trial_data_cell)
                td = trial_data_cell{filenum};

                % trim to just movements
                td = trimTD(td,{'idx_movement_on',0},{'idx_movement_on',11});

                % check to make sure all neurons fire at least once in each condition (pretty rare that one doesn't)
                [~,td_act] = getTDidx(td,'ctrHoldBump',false);
                [~,td_pas] = getTDidx(td,'ctrHoldBump',true);
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

                %find l-naught for muscle length/velocity normalization
                listL0 = zeros(1,numel(td_trim));
                muscleIndex = find(ismember(td_trim(1).muscle_names,muscleArray(muscleNum)));
                for t=1:numel(td_trim)
                    %for bump trials
                    if td_trim(t).ctrHoldBump==1
                        trialL0 = mean(td_trim(t).muscle_len(td_trim(t).idx_startTime:td_trim(t).idx_bumpTime,muscleIndex));
                    elseif td_trim(t).ctrHoldBump==0
                        trialL0 = mean(td_trim(t).muscle_len(td_trim(t).idx_startTime:td_trim(t).idx_goCueTime,muscleIndex));
                    end
                    listL0(1,t) = trialL0;
                end
                L0 = mean(listL0);

                %% find separabilities
                % suppress getTDfields warning...
                getTDfields(td,'time');
                onetime_warn = warning('query','last'); 
                warning('off',onetime_warn.identifier)

                %normalize
                muscleIndex = find(ismember(td_trim(1).muscle_names,muscleArray(muscleNum)));
                for tr=1:numel(td)
                    td(tr).muscle_len(:,muscleNum) = td(tr).muscle_len(:,muscleNum)./L0;
                    td(tr).muscle_vel(:,muscleNum) = td(tr).muscle_vel(:,muscleNum)./L0;

                    %nonlinearize if muscNonlin model
                    if strcmpi(model_aliases,'muscNonlin')
                        if curExp<1
                            td(tr).muscle_vel(:,muscleNum) = sign(td(tr).muscle_vel(:,muscleNum)).*(abs(td(tr).muscle_vel(:,muscleNum))).^curExp;
                        elseif curExp>1
                            td(tr).muscle_vel(:,muscleNum) = (td(tr).muscle_vel(:,muscleNum)).^curExp;
                        end
                    end
                end

                sepResults = actpasSep(td,struct(...
                    'neural_signals',[arrayname '_FR'],...
                    'model_aliases',{model_aliases},...
                    'muscle',muscle));

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
            
            % save each muscle neureval
            currentMuscleNeurEval = table2array(neuron_eval(:,6:end));
            allMuscNeurEval = cat(3,allMuscNeurEval,currentMuscleNeurEval);

            % save each muscle into 3d array
            currentMuscleTable = table2array(avg_neuron_eval(:,5:end));
            allMuscAvgNeurEval = cat(3,allMuscAvgNeurEval,currentMuscleTable);
        end

        multiExp4dArray(:,:,:,d) = allMuscAvgNeurEval(:,:,:,1);
        multiExp4dNeurEval(:,:,:,d) = allMuscNeurEval(:,:,:,1);
    end
    
    temp = avg_neuron_eval;
    temp.monkey=[]; temp.date=[]; temp.task=[]; temp.signalID=[];
    columnNames = temp.Properties.VariableNames;