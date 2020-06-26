%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% actpas_analysis - 
%%      script to run active-passive task analyses in manuscript
%%      and plot out relevant details
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up meta info
    filenames = {'Han_20171201_COactpas_5ms.mat'};

    % plotting variables
%     monkey_names = {'C','Han'};
    monkey_names = {'Han'};
    arrayname = 'S1';
    neural_signals = [arrayname '_FR'];
    session_colors = [...
        102,194,165;...
        252,141,98;...
        141,160,203]/255;
    
%% Loop through trial data files to clean them up
    trial_data_cell = cell(1,length(filenames));
    for filenum = 1:length(filenames)
        %% load and preprocess data
%         td = load(fullfile(datadir,[filenames{filenum}]));
        td = load(filenames{filenum});
    
        % rename trial_data for ease
        td = td.trial_data;
        
        %add ctrHoldBump (boolean)
        for i=1:numel(td.bumpDir)  
            if ~isnan(td.bumpDir(i))
                td.ctrHoldBump(i) = true;
            else
                td.ctrHoldBump(i) = false;
            end
        end
    
        % first process marker data
        % find times when markers are NaN and replace with zeros temporarily
        for trialnum = 1:length(td)
            markernans = isnan(td(trialnum).markers);
            td(trialnum).markers(markernans) = 0;
            td(trialnum) = smoothSignals(td(trialnum),struct('signals','markers'));
            td(trialnum).markers(markernans) = NaN;
            clear markernans
        end
    
        % get marker velocity and acceleration
        td = getDifferential(td,struct('signals','markers','alias','marker_vel'));
        td = getDifferential(td,struct('signals','marker_vel','alias','marker_accel'));
        
        % get speed and ds
        td = getNorm(td,struct('signals','vel','field_extra','_norm'));
        td = getDifferential(td,struct('signals','vel_norm','alias','dvel_norm'));
        
        % remove unsorted neurons
        unit_ids = td(1).([arrayname '_unit_guide']);
        unsorted_units = (unit_ids(:,2)==0);
        new_unit_guide = unit_ids(~unsorted_units,:);
        for trialnum = 1:length(td)
            td(trialnum).(sprintf('%s_unit_guide',arrayname)) = new_unit_guide;
            
            spikes = td(trialnum).(sprintf('%s_spikes',arrayname));
            spikes(:,unsorted_units) = [];
            td(trialnum).(sprintf('%s_spikes',arrayname)) = spikes;
        end
    
        % prep trial data by getting only rewards and trimming to only movements
        % split into trials
        td = splitTD(...
            td,...
            struct(...
                'split_idx_name','idx_startTime',...
                'linked_fields',{{...
                    'trialID',...
                    'result',...
                    'bumpDir',...
                    'tgtDir',...
                    'ctrHoldBump'...
%                     'ctrHold',...
                    }},...
                'start_name','idx_startTime',...
                'end_name','idx_endTime'));
        [~,td] = getTDidx(td,'result','R');
        td = reorderTDfields(td);
        
        % clean nans out...?
        nanners = isnan(cat(1,td.tgtDir));
        td = td(~nanners);
        fprintf('Removed %d trials because of missing target direction\n',sum(nanners))
        biggers = cat(1,td.ctrHoldBump) & abs(cat(1,td.bumpDir))>360;
        td = td(~biggers);
        fprintf('Removed %d trials because bump direction makes no sense\n',sum(biggers))
    
        % remove trials where markers aren't present
        bad_trial = false(length(td),1);
        for trialnum = 1:length(td)
            if any(any(isnan(td(trialnum).markers)))
                bad_trial(trialnum) = true;
            end
        end
        td(bad_trial) = [];
        fprintf('Removed %d trials because of missing markers\n',sum(bad_trial))
        
        % remove trials where muscles aren't present
        bad_trial = false(length(td),1);
        for trialnum = 1:length(td)
            if any(any(isnan(td(trialnum).muscle_len) | isnan(td(trialnum).muscle_vel)))
                bad_trial(trialnum) = true;
            end
        end
        td(bad_trial) = [];
        fprintf('Removed %d trials because of missing muscles\n',sum(bad_trial))
        
        % for C_20170912, trial structure is such that active and passive are part of the same trial--split it up
        if strcmpi(td(1).monkey,'C') && contains(td(1).date_time,'2017/9/12')
            td_copy = td;
            [td_copy.ctrHoldBump] = deal(false);
            td = cat(2,td,td_copy);
            clear td_copy
        end
        
        % split into active and passive
        [~,td_act] = getTDidx(td,'ctrHoldBump',false);
        [~,td_pas] = getTDidx(td,'ctrHoldBump',true);
    
        % find the relevant movmement onsets
        td_act = getMoveOnsetAndPeak(td_act,struct(...
            'start_idx','idx_goCueTime',...
            'start_idx_offset',20,...
            'peak_idx_offset',20,...
            'end_idx','idx_endTime',...
            'method','peak',...
            'peak_divisor',10,...
            'min_ds',1));
        td_pas = getMoveOnsetAndPeak(td_pas,struct(...
            'start_idx','idx_bumpTime',...
            'start_idx_offset',-5,... % give it some wiggle room
            'peak_idx_offset',-5,... % give it some wiggle room
            'end_idx','idx_goCueTime',...
            'method','peak',...
            'peak_divisor',10,...
            'min_ds',1));
        % throw out all trials where bumpTime and movement_on are more than 3 bins apart
        bad_trial = isnan(cat(1,td_pas.idx_movement_on)) | abs(cat(1,td_pas.idx_movement_on)-cat(1,td_pas.idx_bumpTime))>3;
        td_pas = td_pas(~bad_trial);
        fprintf('Removed %d trials because of bad movement onset\n',sum(bad_trial))
    
        % even out sizes and put back together
        minsize = min(length(td_act),length(td_pas));
        td_act = td_act(1:minsize);
        td_pas = td_pas(1:minsize);
        td_trim = cat(2,td_act,td_pas);
    
        % remove low firing neurons
        td_trim = removeBadNeurons(td_trim,struct(...
            'min_fr',1,...
            'fr_window',{{'idx_movement_on',0;'idx_movement_on',11}},...
            'calc_fr',true));
    
        trial_data_cell{filenum} = td_trim;
    end

%% Plot trial info (hand speed and example rasters)
    for filenum = 1:length(trial_data_cell)
        %% load and preprocess data
        td = trial_data_cell{filenum};
    
        % trim to just movements
        td = trimTD(td,{'idx_movement_on',-50},{'idx_movement_on',60});
    
        %% Plot out hand speed
        figure('defaultaxesfontsize',18)
        for trial = 1:length(td)
            timevec = ((1:length(td(trial).vel_norm))-td(trial).idx_movement_on)*td(trial).bin_size;
            if td(trial).ctrHoldBump
                plot(timevec,td(trial).vel_norm,'r')
            else
                plot(timevec,td(trial).vel_norm,'k')
            end
            hold on
        end
        plot(zeros(2,1),ylim,'--k','linewidth',2)
        hold on
        plot(repmat(0.12,2,1),ylim,'--k','linewidth',2)
        xlabel('Time from movement onset (s)')
        ylabel('Hand speed (cm/s)')
        set(gca,'box','off','tickdir','out','xtick',[-0.5 0 0.12 0.5])
        set(gcf,'renderer','Painters')
        suptitle(sprintf('Monkey %s %s',td(1).monkey, td(1).date_time))
        
        %% Plot out hand acceleration
        figure('defaultaxesfontsize',18)
        for trial = 1:length(td)
            timevec = ((1:length(td(trial).vel_norm))-td(trial).idx_movement_on)*td(trial).bin_size;
            if td(trial).ctrHoldBump
%                 subplot(2,1,1);
                plot(timevec,td(trial).dvel_norm,'r')
            end
            hold on
        end
        plot(zeros(2,1),ylim,'--k','linewidth',2)
        hold on
        plot(repmat(0.12,2,1),ylim,'--k','linewidth',2)
        xlabel('Time from movement onset (s)')
        ylabel('Hand acceleration (cm/s2)')
        set(gca,'box','off','tickdir','out','xtick',[-0.5 0 0.12 0.5])
        set(gcf,'renderer','Painters')
        suptitle(sprintf('Monkey %s %s',td(1).monkey, td(1).date_time))
        
        figure('defaultaxesfontsize',18)
        for trial = 1:length(td)
            timevec = ((1:length(td(trial).vel_norm))-td(trial).idx_movement_on)*td(trial).bin_size;
            if ~td(trial).ctrHoldBump
%                 subplot(2,1,2);
                plot(timevec,td(trial).dvel_norm,'k')
            end
            hold on
        end
        plot(zeros(2,1),ylim,'--k','linewidth',2)
        hold on
        plot(repmat(0.12,2,1),ylim,'--k','linewidth',2)
        xlabel('Time from movement onset (s)')
        ylabel('Hand acceleration (cm/s2)')
        set(gca,'box','off','tickdir','out','xtick',[-0.5 0 0.12 0.5])
        set(gcf,'renderer','Painters')
        suptitle(sprintf('Monkey %s %s',td(1).monkey, td(1).date_time))
    
        %% Plot out example rasters for each direction
        dirs = unique(cat(1,td.tgtDir));
        figure('defaultaxesfontsize',18)
        for dirnum = 1:length(dirs)
            % pick a random active and random passive trial with this direction
            act_idx = getTDidx(td,'tgtDir',dirs(dirnum),'ctrHoldBump',false,'rand',1);
            pas_idx = getTDidx(td,'bumpDir',dirs(dirnum),'ctrHoldBump',true,'rand',1);
            td_temp = td([act_idx pas_idx]);
    
            for trialnum = 1:length(td_temp)
                spikes = getSig(td_temp(trialnum),'S1_spikes')';
                timevec = ((1:size(spikes,2))-td_temp(trialnum).idx_movement_on)*td_temp(trialnum).bin_size;
                % active on left, passive on right
                subplot(length(dirs),length(td_temp),(dirnum-1)*length(td_temp)+trialnum)
                % neurons
                for neuronnum = 1:size(spikes,1)
                    spike_times = timevec(spikes(neuronnum,:)>0);
                    scatter(spike_times,repmat(neuronnum,size(spike_times)),5,'k','filled')
                    hold on
                end
                plot(zeros(1,2),[0 size(spikes,1)+1],'--k')
                plot(ones(1,2)*0.12,[0 size(spikes,1)+1],'--k')
                xlabel('Time from movement onset (s)')
                set(gca,'box','off','tickdir','out','xtick',[-0.5 0 0.12 0.5],'ytick',[])
            end
            subplot(length(dirs),length(td_temp),(dirnum-1)*length(td_temp)+1)
            ylabel(sprintf('Direction %f',dirs(dirnum)))
        end
        suptitle(sprintf('Monkey %s %s',td(1).monkey, td(1).date_time))
    end

%% Loop through results to pull out relevant info

model_aliases = {'muscLin'};
models_to_plot = {neural_signals,'muscLin'};
model_titles = {'Actual Firing','Muscle Kin'};

%for linear 
if max(strcmp(model_aliases,{'muscLin'}))==1
    [multiExp4dArray,multiExp4dNeurEval,linColumnNames] = nonlinearModelLoop(trial_data_cell,td_trim,struct(...
                        'arrayname',arrayname,...
                        'model_aliases',{model_aliases},...
                        'models_to_plot',{models_to_plot}));
end

%for nonlinear kinematics (normalized) loop
if max(strcmp(model_aliases,{'muscNonlin'}))==1
    exponents = [0.5];
    [multiExp4dArray,multiExp4dNeurEval,nonlinColumnNames] = nonlinearModelLoop(trial_data_cell,td_trim,struct(...
                        'arrayname',arrayname,...
                        'model_aliases',{model_aliases},...
                        'models_to_plot',{models_to_plot},...
                        'exponents',exponents));
end

%for emg loop
if max(strcmp(model_aliases,{'muscEMG'}))==1 || max(strcmp(model_aliases,{'EMGonly'}))==1
    %muscleArray must match emgMuscleArray
    bumpDuration = 0.125/td(1).bin_size; %(not in this TD, but can get it from CDS)
    endofMove = 18; %window (bins) between end of bump and end of movement
    [emgMuscAvgNeurEval,emgMuscNeurEval,columnNames] = emgModelLoop(trial_data_cell,td_trim,struct(...
                        'arrayname',arrayname,...
                        'model_aliases',{model_aliases},...
                        'models_to_plot',{models_to_plot},...
                        'bumpDuration',bumpDuration,...
                        'endofMove',endofMove));
end

%% evaluate kin vs emg

% if resulting table from previous section is for one model only,
% will need to combine with another model's table to compare the two

    %load the saved tables you want (same dimensions, 53x8x21)
    %load('linearNormalizedOnlyTable.mat')
    %load('emgMuscAvgNeurEval-kinVsEMG.mat')
    
    %load the corresponding neur_eval tables (5300x8x21)
    
    %make any necessary changes
    %removed 'lat_dorsi_sup','lat_dorsi_inf' bec only one lat emg muscle
    linearMuscleTableMod = muscLinAvgNeurEval;
    linearMuscleTableMod(:,:,[15,17]) = [];
    linearMuscleNeurEvalMod = muscLinNeurEval;
    linearMuscleNeurEvalMod(:,:,[15,17]) = [];
    
    %two models to be combined
    model1table = linearMuscleTableMod;
    model1neureval = linearMuscleNeurEvalMod;
    model1vars = linColumnNames;
    model2table = emgMuscAvgNeurEval;
    model2neureval = emgMuscNeurEval;
    model2vars = columnNames;
    
    combinedTable = zeros(size(model2table,1),16,size(model2table,3));
    combinedNeurEval = zeros(size(model2neureval,1),16,size(model2neureval,3));
    combinedVars = strings(1,16);
    j = 1;
    for i=1:8
        combinedTable(:,j,:) = model1table(:,i,:);
        combinedNeurEval(:,j,:) = model1neureval(:,i,:);
        combinedVars(1,j) = model1vars(1,i);
        j=j+1;
        combinedTable(:,j,:) = model2table(:,i,:);
        combinedNeurEval(:,j,:) = model2neureval(:,i,:);
        combinedVars(1,j) = model2vars(1,i);
        j=j+1;
    end
    
% evaluate kinematics vs EMG
    muscAvgNeurEval = combinedTable; %set correct table (containing eval from 2 models combined)
    muscNeurEval = combinedNeurEval; %set correct neur eval (for stats)
    columnNames = combinedVars; %set correct vars (16 cols)
    numMuscles = 5; %find numMuscles best muscles pR2 for each unit
    doPlots = true; %scatter plots, no stats

    [emgNeurCondensedTable,emgCondensedNeurEval] = kinVsEMG(muscAvgNeurEval,muscNeurEval,columnNames,struct(...
                                'numMuscles',numMuscles,...
                                'doPlots',doPlots));
                            
%save condensedTable, condensedNeurEval, and combinedVars

%% evaluate lin vs nonlinear

% if resulting table from previous section is for one model only,
% will need to combine with another model's table to compare the two

    %load the saved tables you want (same dimensions, 53x8x21xexponents)
    %load('linearNormalizedOnlyTable.mat')
    %load('multiExp4dArray_052820.mat')
    multiExp4dArray2Models = zeros(size(multiExp4dArray,1),16,...
                                    size(multiExp4dArray,3),size(multiExp4dArray,4));
    for d=1:size(multiExp4dArray,4)
        combinedTable = zeros(size(multiExp4dArray,1),16,size(multiExp4dArray,3));
        j = 1;
        for i=1:8
            combinedTable(:,j,:) = linearMuscleTable(:,i,:);
            j=j+1;
            combinedTable(:,j,:) = multiExp4dArray(:,i,:,d);
            j=j+1;
        end
        multiExp4dArray2Models(:,:,:,d) = combinedTable;
    end

% make necessary changes to table
    dontWant = [1,2,7,11,12,16,17,20,22,27,30,31,32,33,34,36]; %REMOVE non-emg muscles from allMuscAvgNeurEval
    multiExp4dArray2ModelsMod = multiExp4dArray2Models;
    multiExp4dArray2ModelsMod(:,:,dontWant,:) = [];
    
% show top muscles pR2 per category bar plots (pick one exponent at a time)
    numMuscles = 5; %find numMuscles best muscles pR2 for each unit
    expInd = 1;
    muscAvgNeurEval = multiExp4dArray2ModelsMod(:,:,:,expInd);

    %bar plots showing top 5 muscles per unit per category (for specified exponent)
    topMuscBarPlots(muscAvgNeurEval,struct(...
                        'dontWant',dontWant,...
                        'numMuscles',numMuscles));

% evaluate kinematics vs EMG  
    multiExpArrayInput = multiExp4dArray2ModelsMod;
    numMuscles = 5; %find numMuscles best muscles pR2 for each unit
    doPlots = true; %scatter plots, no stats
    
    linVsNonlin(multiExpArrayInput,columnNames,struct(...
                    'numMuscles',numMuscles,...
                    'doPlots',doPlots));

%% make plots for condensed neurEval tables w stats!

model_pairs = {'muscLin','muscEMG'};
models_to_plot = {neural_signals,'muscLin','muscEMG'};
model_titles = {'Muscle Kinematics','Muscle Kin and EMG'};

% compare pR2 of handelbow vs ext
figure('defaultaxesfontsize',18)
for pairnum = 1:size(model_pairs,1)
    for monkeynum = 1:length(monkey_names)
        % set subplot...
        subplot(size(model_pairs,1),length(monkey_names),...
            (pairnum-1)*length(monkey_names)+monkeynum)
        plot([-0.4 0.6],[-0.4 0.6],'k--','linewidth',0.5)
        hold on
        plot([0 0],[-0.4 0.6],'k-','linewidth',0.5)
        plot([-0.4 0.6],[0 0],'k-','linewidth',0.5)

%             % get sessions
%             [~,monkey_evals] = getNTidx(neuron_eval,'monkey',monkey_names{monkeynum});
%             session_dates = unique(monkey_evals.date);
% 
%             for sessionnum = 1:length(session_dates)
%                 [~,session_evals] = getNTidx(monkey_evals,'date',session_dates{sessionnum});
            sessionnum = 1;
            avgNeurEval = emgNeurCondensedTable;
            neurEval = emgCondensedNeurEval;
            colNames = combinedVars;
            pr2_winners = compareEncoderMetricsMod(avgNeurEval,neurEval,colNames,struct(...
                'bonferroni_correction',6,...
                'models',{models_to_plot},...
                'model_pairs',{model_pairs},...
                'postfix','_eval',...
                'num_repeats',20,...
                'num_folds',5));

%                 [~,avg_pR2] = getNTidx(avg_neuron_eval,'monkey',monkey_names{monkeynum},'date',session_dates{sessionnum});
            avg_pR2 = avgNeurEval;

            % scatter filled circles if there's a winner, empty circles if not
            no_winner =  cellfun(@isempty,pr2_winners(pairnum,:));
            idx1 = strcmp(colNames,strcat(model_pairs{1},'_eval'));
            idx2 = strcmp(colNames,strcat(model_pairs{2},'_eval'));
            scatterlims(...
                [-0.4 0.6],...
                [-0.4 0.6],...
                avg_pR2(no_winner,idx1),...
                avg_pR2(no_winner,idx2),...
                [],session_colors(sessionnum,:))
            scatterlims(...
                [-0.4 0.6],...
                [-0.4 0.6],...
                avg_pR2(~no_winner,idx1),...
                avg_pR2(~no_winner,idx2),...
                [],session_colors(sessionnum,:),'filled')
%             end

        % make axes pretty
        set(gca,'box','off','tickdir','out')
        axis image
        if monkeynum ~= 1 || pairnum ~= 1
            set(gca,'box','off','tickdir','out',...
                'xtick',[],'ytick',[])
        end
        xlabel(sprintf('%s pR2',model_titles{pairnum,1}))
        ylabel(sprintf('%s pR2',model_titles{pairnum,2}))
    end
end

%% Plot within condition vs across condition pR2 for each neuron in all sessions
conds = {'act','pas'};
model_pairs = {'muscLin','muscEMG'};
models_to_plot = {neural_signals,'muscLin','muscEMG'};
model_titles = {'Muscle Kinematics','Muscle Kin and EMG'};
    
for modelnum = 2:length(models_to_plot)
    figure('defaultaxesfontsize',18)
    for monkeynum = 1:length(monkey_names)
        for condnum = 1:2
            % set subplot
            subplot(2,length(monkey_names),(condnum-1)*length(monkey_names)+monkeynum)
            plot([-0.7 0.7],[-0.7 0.7],'k--','linewidth',0.5)
            hold on
            plot([0 0],[-0.7 0.7],'k-','linewidth',0.5)
            plot([-0.7 0.7],[0 0],'k-','linewidth',0.5)

            % get sessions
%             [~,monkey_evals] = getNTidx(neuron_eval,'monkey',monkey_names{monkeynum});
%             session_dates = unique(monkey_evals.date);
% 
%             % plot out each session
%             for sessionnum = 1:length(session_dates)
%                 [~,avg_pR2] = getNTidx(avg_neuron_eval,'monkey',monkey_names{monkeynum},'date',session_dates{sessionnum});
% 
%                 [~,session_evals] = getNTidx(monkey_evals,'date',session_dates{sessionnum});
                sessionnum = 1;
                avgNeurEval = emgNeurCondensedTable;
                neurEval = emgCondensedNeurEval;
                colNames = combinedVars;

                pr2_winners = compareEncoderMetricsMod(avgNeurEval,neurEval,colNames,struct(...
                    'bonferroni_correction',6,...
                    'models',{models_to_plot},...
                    'model_pairs',{model_pairs},...
                    'postfix','_eval',...
                    'num_repeats',20,...
                    'num_folds',5));

                % fill by whether separable or not?
%                 sig_seps = avg_pR2.S1_FR_indiv_sep_CI_lo > 0.5;

                % fill by whether winner of comparison with ext or not?
                no_winner =  cellfun(@isempty,pr2_winners(1,:));
                avg_pR2 = avgNeurEval;
                idx1 = strcmp(colNames,sprintf('%s_%s_eval',models_to_plot{modelnum},conds{condnum}));
                idx2 = strcmp(colNames,sprintf('%s_train_%s_eval',models_to_plot{modelnum},conds{condnum}));

                scatterlims(...
                    [-0.7 0.7],...
                    [-0.7 0.7],...
                    avg_pR2(no_winner,idx1),...
                    avg_pR2(no_winner,idx2),...
                    [],session_colors(sessionnum,:),'filled')
                scatterlims(...
                    [-0.7 0.7],...
                    [-0.7 0.7],...
                    avg_pR2(~no_winner,idx1),...
                    avg_pR2(~no_winner,idx2),...
                    [],session_colors(sessionnum,:),'filled')
%             end
            % make axes pretty
            set(gca,'box','off','tickdir','out',...
                'xlim',[-0.7 0.7],'ylim',[-0.7 0.7])
            axis equal
            % if monkeynum ~= 1 || condnum ~= 1
            %     set(gca,'box','off','tickdir','out',...
            %         'xtick',[],'ytick',[])
            % end
            xlabel(sprintf('%s pR2, trained full, tested %s',model_titles{modelnum-1},conds{condnum}),'FontSize',15)
            ylabel(sprintf('%s pR2, trained %s, tested %s',model_titles{modelnum-1},conds{condnum},conds{condnum}),'FontSize',10)
        end
    end
    % suptitle('Full pR^2 vs within condition pR^2')
end

    
    
%% make plots
    % plot separability of each neuron and save CIs into avg_neuron_eval
    for monkeynum = 1:length(monkey_names)
        [~,monkey_evals] = getNTidx(neuron_eval,'monkey',monkey_names{monkeynum});
        session_dates = unique(monkey_evals.date);
        for sessionnum = 1:length(session_dates)
            [~,session_evals] = getNTidx(monkey_evals,'date',session_dates{sessionnum});
            [~,session_avg] = getNTidx(avg_neuron_eval,'monkey',monkey_names{monkeynum},'date',session_dates{sessionnum});

            % create place to save CIs
            signalID = session_avg.signalID;
            S1_FR_indiv_sep_CI_lo = zeros(size(signalID,1),1);
            S1_FR_indiv_sep_CI_hi = zeros(size(signalID,1),1);

            figure('defaultaxesfontsize',18)
            plot([0 size(signalID,1)+1],[0.5 0.5],'--k','linewidth',2)
            hold on
            for neuronnum = 1:size(signalID,1)
                [~,single_neuron_eval] = getNTidx(session_evals,'signalID',signalID(neuronnum,:));

                % figure out error bars and save
                [CI_lo,CI_hi] = crossval_errorbars(single_neuron_eval.S1_FR_indiv_sep,struct(...
                    'num_repeats',double(max(single_neuron_eval.crossvalID(:,1))),...
                    'num_folds',double(max(single_neuron_eval.crossvalID(:,2)))));
                S1_FR_indiv_sep_CI_lo(neuronnum,:) = CI_lo;
                S1_FR_indiv_sep_CI_hi(neuronnum,:) = CI_hi;
                
                % scatter(repmat(neuronnum,1,height(single_neuron_eval)),single_neuron_eval.S1_FR_indiv_sep,25,'k','filled','markerfacealpha',0.2)
                plot(repmat(neuronnum,1,2),[CI_lo CI_hi],'-k','linewidth',2)
                scatter(neuronnum,session_avg.S1_FR_indiv_sep(neuronnum,:),100,'k','filled')
            end
            set(gca,'box','off','tickdir','out','ylim',[0 1],'xlim',[0 size(signalID,1)+1])

            % save error bars into avg_neuron_eval table
            CI_cell{monkeynum,sessionnum} = table(...
                session_avg.monkey,...
                session_avg.date,...
                signalID,...
                S1_FR_indiv_sep_CI_lo,...
                S1_FR_indiv_sep_CI_hi,...
                'VariableNames',{'monkey','date','signalID','S1_FR_indiv_sep_CI_lo','S1_FR_indiv_sep_CI_hi'});
        end
    end
    CI_table = vertcat(CI_cell{:});
    avg_neuron_eval = join(avg_neuron_eval,CI_table);

    % compare pR2 of handelbow vs ext
    figure('defaultaxesfontsize',18)
    model_pairs = {'extAccel','markers_pcaAccel'};
    for pairnum = 1:size(model_pairs,1)
        for monkeynum = 1:length(monkey_names)
            % set subplot...
            subplot(size(model_pairs,1),length(monkey_names),...
                (pairnum-1)*length(monkey_names)+monkeynum)
            plot([-0.4 0.6],[-0.4 0.6],'k--','linewidth',0.5)
            hold on
            plot([0 0],[-0.4 0.6],'k-','linewidth',0.5)
            plot([-0.4 0.6],[0 0],'k-','linewidth',0.5)

            % get sessions
            [~,monkey_evals] = getNTidx(neuron_eval,'monkey',monkey_names{monkeynum});
            session_dates = unique(monkey_evals.date);

            for sessionnum = 1:length(session_dates)
                [~,session_evals] = getNTidx(monkey_evals,'date',session_dates{sessionnum});
                pr2_winners = compareEncoderMetrics(session_evals,struct(...
                    'bonferroni_correction',6,...
                    'models',{models_to_plot},...
                    'model_pairs',{model_pairs},...
                    'postfix','_eval'));

                [~,avg_pR2] = getNTidx(avg_neuron_eval,'monkey',monkey_names{monkeynum},'date',session_dates{sessionnum});
                % scatter filled circles if there's a winner, empty circles if not
                no_winner =  cellfun(@isempty,pr2_winners(pairnum,:));
                scatterlims(...
                    [-0.4 0.6],...
                    [-0.4 0.6],...
                    avg_pR2.(strcat(model_pairs{pairnum,1},'_eval'))(no_winner),...
                    avg_pR2.(strcat(model_pairs{pairnum,2},'_eval'))(no_winner),...
                    [],session_colors(sessionnum,:))
                scatterlims(...
                    [-0.4 0.6],...
                    [-0.4 0.6],...
                    avg_pR2.(strcat(model_pairs{pairnum,1},'_eval'))(~no_winner),...
                    avg_pR2.(strcat(model_pairs{pairnum,2},'_eval'))(~no_winner),...
                    [],session_colors(sessionnum,:),'filled')
            end

            % make axes pretty
            set(gca,'box','off','tickdir','out')
            axis image
            if monkeynum ~= 1 || pairnum ~= 1
                set(gca,'box','off','tickdir','out',...
                    'xtick',[],'ytick',[])
            end
            xlabel(sprintf('%s pR2',getModelTitles(model_pairs{pairnum,1})))
            ylabel(sprintf('%s pR2',getModelTitles(model_pairs{pairnum,2})))
        end
    end

    % Plot within condition vs across condition pR2 for each neuron in all sessions
    conds = {'act','pas'};
    model_pairs = {'extAccel','markers_pcaAccel'};
    for modelnum = 2:length(models_to_plot)
        figure('defaultaxesfontsize',18)
        for monkeynum = 1:length(monkey_names)
            for condnum = 1:2
                % set subplot
                subplot(2,length(monkey_names),(condnum-1)*length(monkey_names)+monkeynum)
                plot([-0.7 0.7],[-0.7 0.7],'k--','linewidth',0.5)
                hold on
                plot([0 0],[-0.7 0.7],'k-','linewidth',0.5)
                plot([-0.7 0.7],[0 0],'k-','linewidth',0.5)

                % get sessions
                [~,monkey_evals] = getNTidx(neuron_eval,'monkey',monkey_names{monkeynum});
                session_dates = unique(monkey_evals.date);

                % plot out each session
                for sessionnum = 1:length(session_dates)
                    [~,avg_pR2] = getNTidx(avg_neuron_eval,'monkey',monkey_names{monkeynum},'date',session_dates{sessionnum});
                    
                    [~,session_evals] = getNTidx(monkey_evals,'date',session_dates{sessionnum});
                    pr2_winners = compareEncoderMetrics(session_evals,struct(...
                        'bonferroni_correction',6,...
                        'models',{models_to_plot},...
                        'model_pairs',{model_pairs},...
                        'postfix','_eval'));

                    % fill by whether separable or not?
                    sig_seps = avg_pR2.S1_FR_indiv_sep_CI_lo > 0.5;
                    
                    % fill by whether winner of comparison with ext or not?
                    no_winner =  cellfun(@isempty,pr2_winners(1,:));

                    scatterlims(...
                        [-0.7 0.7],...
                        [-0.7 0.7],...
                        avg_pR2.(sprintf('%s_%s_eval',models_to_plot{modelnum},conds{condnum}))(no_winner),...
                        avg_pR2.(sprintf('%s_train_%s_eval',models_to_plot{modelnum},conds{condnum}))(no_winner),...
                        [],session_colors(sessionnum,:),'filled')
                    scatterlims(...
                        [-0.7 0.7],...
                        [-0.7 0.7],...
                        avg_pR2.(sprintf('%s_%s_eval',models_to_plot{modelnum},conds{condnum}))(~no_winner),...
                        avg_pR2.(sprintf('%s_train_%s_eval',models_to_plot{modelnum},conds{condnum}))(~no_winner),...
                        [],session_colors(sessionnum,:),'filled')
                end
                % make axes pretty
                set(gca,'box','off','tickdir','out',...
                    'xlim',[-0.7 0.7],'ylim',[-0.7 0.7])
                axis equal
                % if monkeynum ~= 1 || condnum ~= 1
                %     set(gca,'box','off','tickdir','out',...
                %         'xtick',[],'ytick',[])
                % end
                xlabel(sprintf('%s pR2, trained full, tested %s',getModelTitles(models_to_plot{modelnum}),conds{condnum}))
                ylabel(sprintf('%s pR2, trained %s, tested %s',getModelTitles(models_to_plot{modelnum}),conds{condnum},conds{condnum}))
            end
        end
        % suptitle('Full pR^2 vs within condition pR^2')
    end

    % plot separability against full and within condition pR2
    conds = {'','act_','pas_'};
    for modelnum = 2:length(models_to_plot)
        figure('defaultaxesfontsize',18)
        for monkeynum = 1:length(monkey_names)
            for condnum = 1:length(conds)
                % set subplot
                subplot(length(monkey_names),length(conds),(monkeynum-1)*length(conds)+condnum)
                plot([0 0],[0 1],'k-','linewidth',0.5)
                hold on
                plot([-0.7 0.7],[0 0],'k-','linewidth',0.5)
                plot([-0.7 0.7],[0.5 0.5],'k--','linewidth',0.5)

                % get sessions
                [~,monkey_evals] = getNTidx(neuron_eval,'monkey',monkey_names{monkeynum});
                session_dates = unique(monkey_evals.date);

                % plot out each session
                for sessionnum = 1:length(session_dates)
                    [~,session_evals] = getNTidx(monkey_evals,'date',session_dates{sessionnum});
                    [~,avg_pR2] = getNTidx(avg_neuron_eval,'monkey',monkey_names{monkeynum},'date',session_dates{sessionnum});

                    sig_seps = avg_pR2.S1_FR_indiv_sep_CI_lo > 0.5;

                    scatterlims(...
                        [-0.7 0.7],...
                        [0.4 1],...
                        avg_pR2.(sprintf('%s_%seval',models_to_plot{modelnum},conds{condnum}))(~sig_seps),...
                        avg_pR2.S1_FR_indiv_sep(~sig_seps),...
                        [],session_colors(sessionnum,:),'filled')
                    scatterlims(...
                        [-0.7 0.7],...
                        [0.4 1],...
                        avg_pR2.(sprintf('%s_%seval',models_to_plot{modelnum},conds{condnum}))(sig_seps),...
                        avg_pR2.S1_FR_indiv_sep(sig_seps),...
                        [],session_colors(sessionnum,:),'filled')

                    % fit quick linear model to plot fit line
                    % lm = fitlm(...
                    %     avg_pR2.(sprintf('%s_%seval',models_to_plot{modelnum},conds{condnum})),...
                    %     avg_pR2.S1_FR_indiv_sep);
                    % plot([-1;1],lm.predict([-1;1]),...
                    %     '--','color',session_colors(sessionnum,:),'linewidth',1)
                end
                % make axes pretty
                set(gca,'box','off','tickdir','out','xtick',[-0.7 0.7])
                if condnum ~= 1
                    set(gca,'ytick',[])
                else
                    ylabel('Neural Separability')
                end
                if monkeynum == length(monkey_names)
                    xlabel(sprintf('%s %s pR2',getModelTitles(models_to_plot{modelnum}),conds{condnum}))
                end
                axis image
                set(gca,'ylim',[0.4 1])
            end
        end
        suptitle('Neural separability vs pR^2')
    end

    % get correlation values for each crossval
    keycols = {'monkey','date','task','crossvalID'};
    keyTable = unique(neuron_eval(:,keycols));
    corr_cell = cell(height(keyTable),1);
    for key_idx = 1:height(keyTable)
        key = keyTable(key_idx,:);
        cond_idx = ismember(neuron_eval(:,keycols),key);
        neuron_eval_select = neuron_eval(cond_idx,:);

        % get correlations
        model_corr = cell(1,length(models_to_plot)-1);
        for modelnum = 2:length(models_to_plot)
            corr_pr2_neuronsep = corr(neuron_eval_select.S1_FR_indiv_sep,neuron_eval_select.(sprintf('%s_eval',models_to_plot{modelnum})),'rows','complete');
            corr_actpr2_neuronsep = corr(neuron_eval_select.S1_FR_indiv_sep,neuron_eval_select.(sprintf('%s_act_eval',models_to_plot{modelnum})),'rows','complete');
            corr_paspr2_neuronsep = corr(neuron_eval_select.S1_FR_indiv_sep,neuron_eval_select.(sprintf('%s_pas_eval',models_to_plot{modelnum})),'rows','complete');

            model_corr{modelnum-1} = table(...
                corr_pr2_neuronsep,...
                corr_actpr2_neuronsep,...
                corr_paspr2_neuronsep,...
                'VariableNames',strcat(models_to_plot{modelnum},{...
                    '_corr_pr2_neuronsep',...
                    '_corr_actpr2_neuronsep',...
                    '_corr_paspr2_neuronsep'}));
        end

        % put together in table
        corr_cell{key_idx} = horzcat(model_corr{:});
    end
    neuron_corr_table = horzcat(keyTable,vertcat(corr_cell{:}));

    % make figure for correlations of pR2 handelbow model with separability
    alpha = 0.05;
    xvals = [2 5 8]/10;
    for modelnum = 2:length(models_to_plot)
        figure('defaultaxesfontsize',18)
        for monkeynum = 1:length(monkey_names)
            subplot(1,length(monkey_names),monkeynum)
            plot([min(xvals)-0.2 max(xvals)+0.2],[0 0],'-k','linewidth',2)
            hold on
            
            % figure out what sessions we have for this monkey
            [~,monkey_corrs] = getNTidx(neuron_corr_table,'monkey',monkey_names{monkeynum});
            session_dates = unique(monkey_corrs.date);

            for sessionnum = 1:length(session_dates)
                [~,session_corrs] = getNTidx(monkey_corrs,'date',session_dates{sessionnum});

                % estimate error bars
                [~,cols] = ismember(...
                    strcat(models_to_plot{modelnum},{'_corr_pr2_neuronsep','_corr_actpr2_neuronsep','_corr_paspr2_neuronsep'}),...
                    session_corrs.Properties.VariableNames);
                num_repeats = double(max(session_corrs.crossvalID(:,1)));
                num_folds = double(max(session_corrs.crossvalID(:,2)));
                crossval_correction = 1/(num_folds*num_repeats) + 1/(num_folds-1);
                yvals = mean(session_corrs{:,cols});
                var_corrs = var(session_corrs{:,cols});
                upp = tinv(1-alpha/2,num_folds*num_repeats-1);
                low = tinv(alpha/2,num_folds*num_repeats-1);
                CI_lo = yvals + low * sqrt(crossval_correction*var_corrs);
                CI_hi = yvals + upp * sqrt(crossval_correction*var_corrs);
                
                % plot dots and lines
                plot(repmat(xvals,2,1),[CI_lo;CI_hi],'-','linewidth',2,'color',session_colors(sessionnum,:))
                scatter(xvals(:),yvals(:),50,session_colors(sessionnum,:),'filled')
            end
            ylabel('Correlation with active/passive separability')
            title(sprintf('Monkey %s',monkey_names{monkeynum}))
            set(gca,'box','off','tickdir','out',...
                'xlim',[min(xvals)-0.2 max(xvals)+0.2],...
                'xtick',xvals,'xticklabel',{'Full pR^2','Active pR^2','Passive pR^2'},...
                'ylim',[-1 1],'ytick',[-1 -0.5 0 0.5 1])
        end
        suptitle(models_to_plot{modelnum})
    end

