function [emgNeurCondensedTable,emgCondensedNeurEval] = kinVsEMG(muscAvgNeurEval,muscNeurEval,columnNames,params)

    %default params
    muscleArray = {'bicep_lh','bicep_sh','brachialis',...
        'brachioradialis','deltoid_ant','deltoid_med','deltoid_pos',...
        'ext_carp_rad_brevis','ext_carpi_ulnaris',...
        'ext_digitorum','flex_carpi_radialis','flex_carpi_ulnaris',...
        'flex_digit_superficialis','infraspinatus',...
        'lat_dorsi_cen','pectoralis_sup',...
        'pectoralis_inf',...
        'teres_major','tricep_lat','tricep_lon','tricep_sho'};
    numMuscles = 5; %find numMuscles best muscles pR2 for each unit
    doPlots = true;
    assignParams(who,params); %overwrite params

    %pre-allocate arrays
    kinAct = zeros(1,numel(muscAvgNeurEval(:,1,1)));
    emgAct = zeros(1,numel(muscAvgNeurEval(:,1,1)));
    kinPas = zeros(1,numel(muscAvgNeurEval(:,1,1)));
    emgPas = zeros(1,numel(muscAvgNeurEval(:,1,1)));
    unitMuscles = strings(numMuscles,numel(muscAvgNeurEval(:,1,1)));
    kinEmgCondensedTable = zeros(numel(muscAvgNeurEval(:,1,1)),16,numMuscles);
    kinEmgCondensedNeurEval = zeros(numel(muscNeurEval(:,1,1)),16,numMuscles);

    %loop through each unit
    for n=1:numel(muscAvgNeurEval(:,1,1))

        %find best muscles pR2 evaluated in act
        [maxKinAct, maxKinActInd] = maxk(muscAvgNeurEval(n,3,:),numMuscles);
        bestMuscKinAct = reshape(permute(maxKinAct,[3,2,1]),[1,numMuscles]);
        maxKinActInd = reshape(permute(maxKinActInd,[3,2,1]),[1,numMuscles]);
        [maxEMGAct, maxEMGActInd] = maxk(muscAvgNeurEval(n,4,:),numMuscles);
        bestMuscEMGAct = reshape(permute(maxEMGAct,[3,2,1]),[1,numMuscles]);
        maxEMGActInd = reshape(permute(maxEMGActInd,[3,2,1]),[1,numMuscles]);
        %find best muscles pR2 evaluated in pas
        [maxKinPas, maxKinPasInd] = maxk(muscAvgNeurEval(n,5,:),numMuscles);
        bestMuscKinPas = reshape(permute(maxKinPas,[3,2,1]),[1,numMuscles]);
        maxKinPasInd = reshape(permute(maxKinPasInd,[3,2,1]),[1,numMuscles]);
        [maxEMGPas, maxEMGPasInd] = maxk(muscAvgNeurEval(n,6,:),numMuscles);
        bestMuscEMGPas = reshape(permute(maxEMGPas,[3,2,1]),[1,numMuscles]);
        maxEMGPasInd = reshape(permute(maxEMGPasInd,[3,2,1]),[1,numMuscles]);

        %%FOR COMPARING MEAN PR2 OF SAME MUSCLES ACROSS 4 CATEGORIES
        %find best group of muscles
        meanKinAct = mean(bestMuscKinAct);
        meanEMGAct = mean(bestMuscEMGAct);
        meanKinPas = mean(bestMuscKinPas);
        meanEMGPas = mean(bestMuscEMGPas);
        [~,bestGroupIdx] = maxk([meanKinAct meanEMGAct meanKinPas meanEMGPas],1);
        if bestGroupIdx==1
            checkMuscles = maxKinActInd;
        elseif bestGroupIdx==2
            checkMuscles = maxEMGActInd;
        elseif bestGroupIdx==3
            checkMuscles = maxKinPasInd;
        elseif bestGroupIdx==4
            checkMuscles = maxEMGPasInd;
        end
        bestMuscKinAct = reshape(permute(muscAvgNeurEval(n,3,checkMuscles),[3,2,1]),[1,numMuscles]);
        bestMuscEMGAct = reshape(permute(muscAvgNeurEval(n,4,checkMuscles),[3,2,1]),[1,numMuscles]);
        bestMuscKinPas = reshape(permute(muscAvgNeurEval(n,5,checkMuscles),[3,2,1]),[1,numMuscles]);
        bestMuscEMGPas = reshape(permute(muscAvgNeurEval(n,6,checkMuscles),[3,2,1]),[1,numMuscles]);
        
        %save current unit table (checkMuscles only)
        kinEmgCondensedTable(n,:,:) = muscAvgNeurEval(n,:,checkMuscles);
        %want every 53rd row
        for x=1:100
            kinEmgCondensedNeurEval(n+53*(x-1),:,:) = muscNeurEval(n+53*(x-1),:,checkMuscles);
        end

        %save top muscles for unit
        for c=1:numel(checkMuscles)
            unitMuscles{c,n} = muscleArray{checkMuscles(1,c)};
        end

        %find mean of 5 best muscle group for each category
        meanKinAct = mean(bestMuscKinAct);
        meanEMGAct = mean(bestMuscEMGAct);
        meanKinPas = mean(bestMuscKinPas);
        meanEMGPas = mean(bestMuscEMGPas);

        %add this unit pR2 values to master lists
        kinAct(1,n) = meanKinAct;
        emgAct(1,n) = meanEMGAct;
        kinPas(1,n) = meanKinPas;
        emgPas(1,n) = meanEMGPas;
    end

    emgNeurCondensedTable = mean(kinEmgCondensedTable,3);
    emgCondensedNeurEval = mean(kinEmgCondensedNeurEval,3);
    %emgNeurCondensedTable = table(emgNeurCondensedTable,'VariableNames',columnNames);

    if doPlots
        %scatter plots
        low = -1;
        high = 1;
        step = 0.25;
        figure('Name','Kinematics/EMG Active/Passive comparison across units');
        sgtitle('Kinematics/EMG Active/Passive comparison across units');
        x = [low:high];
        y = x;

        subplot(2,2,1);
        scatter(kinAct,kinPas,'filled')
        title('Kinematics active vs. passive');
        xlabel('Kinematics Active pR2')
        ylabel('Kinematics Passive pR2')
        xlim([low high])
        set(gca,'XTick',low:step:high)
        ylim([low high])
        set(gca,'YTick',low:step:high)
        hold on
        plot(x,y)

        subplot(2,2,2);
        scatter(emgAct,emgPas,'filled')
        title('EMG active vs. passive');
        xlabel('EMG Active pR2')
        ylabel('EMG Passive pR2')
        xlim([low high])
        set(gca,'XTick',low:step:high)
        ylim([low high])
        set(gca,'YTick',low:step:high)
        hold on
        plot(x,y)

        subplot(2,2,3);
        scatter(kinAct,emgAct,'filled')
        title('Active Kinematics vs EMG');
        xlabel('Active Kinematics pR2')
        ylabel('Active EMG pR2')
        xlim([low high])
        set(gca,'XTick',low:step:high)
        ylim([low high])
        set(gca,'YTick',low:step:high)
        hold on
        plot(x,y)

        subplot(2,2,4);
        scatter(kinPas,emgPas,'filled')
        title('Passive Kinematics vs EMG');
        xlabel('Passive Kinematics pR2')
        ylabel('Passive EMG pR2')
        xlim([low high])
        set(gca,'XTick',low:step:high)
        ylim([low high])
        set(gca,'YTick',low:step:high)
        hold on
        plot(x,y)

        %alternative act/pas scatter plots
        %scatter plots
        figure('Name','Active/Passive comparison across units for Kinematics and EMG models');
        sgtitle('Active/Passive comparison across units for Kinematics and EMG models');
        x = [low:high];
        y = x;

        subplot(1,2,1);
        scatter(kinAct,emgPas,'filled')
        title('Kinematics active vs. EMG passive');
        xlabel('Kinematics Active pR2')
        ylabel('EMG Passive pR2')
        xlim([low high])
        set(gca,'XTick',low:step:high)
        ylim([low high])
        set(gca,'YTick',low:step:high)
        hold on
        plot(x,y)

        subplot(1,2,2);
        scatter(emgAct,kinPas,'filled')
        title('EMG active vs. Kinematics passive');
        xlabel('EMG Active pR2')
        ylabel('Kinematics Passive pR2')
        xlim([low high])
        set(gca,'XTick',low:step:high)
        ylim([low high])
        set(gca,'YTick',low:step:high)
        hold on
        plot(x,y)
    end

    clear maxKinAct maxKinActInd bestMuscKinAct
    clear maxEMGAct maxEMGActInd bestMuscEMGAct
    clear maxKinPas maxKinPasInd bestMuscKinPas
    clear maxEMGPas maxEMGPasInd bestMuscEMGPas
    clear maxKinActInd maxKinPasInd maxEMGActInd maxEMGPasInd
    clear bests j check indexes addition muscleLabels x
    clear n c
    clear meanKinAct meanKinPas meanEMGAct meanEMGPas
    clear checkMuscles muscleArray1