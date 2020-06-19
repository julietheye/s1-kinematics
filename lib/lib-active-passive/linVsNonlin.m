function linNonlinCondensedTable4d = linVsNonlin(multiExp4dArray2Models,columnNames,params)

    %default params
    fullMuscleArray = {'abd\_poll\_longus','anconeus','bicep\_lh','bicep\_sh','brachialis',...
    'brachioradialis','coracobrachialis','deltoid\_ant','deltoid\_med','deltoid\_pos',...
    'dorsoepitrochlearis','ext\_carpi\_rad\_longus','ext\_carp\_rad\_brevis','ext\_carpi\_ulnaris',...
    'ext\_digitorum','ext\_digiti','ext\_indicis','flex\_carpi\_radialis','flex\_carpi\_ulnaris',...
    'flex\_digit\_profundus','flex\_digit\_superficialis','flex\_poll\_longus','infraspinatus',...
    'lat\_dorsi\_sup','lat\_dorsi\_cen','lat\_dorsi\_inf','palmaris\_longus','pectoralis\_sup',...
    'pectoralis\_inf','pronator\_quad','pronator\_teres','subscapularis','supinator',...
    'supraspinatus','teres\_major','teres\_minor','tricep\_lat','tricep\_lon','tricep\_sho'};
    dontWant = []; %REMOVE non-emg muscles from allMuscAvgNeurEval
    numMuscles = 5; %find numMuscles best muscles pR2 for each unit
    doPlots = true;
    assignParams(who,params); %overwrite params
    
    %remove dontWant muscles
    fullMuscleArray(dontWant) = [];
    muscleArray = fullMuscleArray;

linNonlinCondensedTable4d = zeros(size(multiExp4dArray2Models,1),16,numMuscles,size(multiExp4dArray,4));

% hack-y way to get scatter plots from normalized models
for d=1:size(multiExp4dArray,4)
    
    muscAvgNeurEval = multiExp4dArray2Models(:,:,:,d);
    
    %pre-allocate arrays
    linAct = zeros(1,numel(muscAvgNeurEval(:,1,1)));
    nonlinAct = zeros(1,numel(muscAvgNeurEval(:,1,1)));
    linPas = zeros(1,numel(muscAvgNeurEval(:,1,1)));
    nonlinPas = zeros(1,numel(muscAvgNeurEval(:,1,1)));
    unitMuscles = strings(numMuscles,numel(muscAvgNeurEval(:,1,1)));
    linNonlinCondensedTable = zeros(numel(muscAvgNeurEval(:,1,1)),16,numMuscles);
    
    %loop through each unit
    for n=1:numel(muscAvgNeurEval(:,1,1))
        
        %find best muscles pR2 evaluated in act
        [maxLinAct, maxLinActInd] = maxk(muscAvgNeurEval(n,3,:),numMuscles);
        bestMuscLinAct = reshape(permute(maxLinAct,[3,2,1]),[1,numMuscles]);
        maxLinActInd = reshape(permute(maxLinActInd,[3,2,1]),[1,numMuscles]);
        [maxNonlinAct, maxNonlinActInd] = maxk(muscAvgNeurEval(n,4,:),numMuscles);
        bestMuscNonlinAct = reshape(permute(maxNonlinAct,[3,2,1]),[1,numMuscles]);
        maxNonlinActInd = reshape(permute(maxNonlinActInd,[3,2,1]),[1,numMuscles]);
        %find best muscles pR2 evaluated in pas
        [maxLinPas, maxLinPasInd] = maxk(muscAvgNeurEval(n,5,:),numMuscles);
        bestMuscLinPas = reshape(permute(maxLinPas,[3,2,1]),[1,numMuscles]);
        maxLinPasInd = reshape(permute(maxLinPasInd,[3,2,1]),[1,numMuscles]);
        [maxNonlinPas, maxNonlinPasInd] = maxk(muscAvgNeurEval(n,6,:),numMuscles);
        bestMuscNonlinPas = reshape(permute(maxNonlinPas,[3,2,1]),[1,numMuscles]);
        maxNonlinPasInd = reshape(permute(maxNonlinPasInd,[3,2,1]),[1,numMuscles]);
        
        %%FOR COMPARING MEAN PR2 OF SAME MUSCLES ACROSS 4 CATEGORIES
        %find best group of muscles
        meanLinAct = mean(bestMuscLinAct);
        meanNonlinAct = mean(bestMuscNonlinAct);
        meanLinPas = mean(bestMuscLinPas);
        meanNonlinPas = mean(bestMuscNonlinPas);
        [bestGroup,bestGroupIdx] = maxk([meanLinAct meanNonlinAct meanLinPas meanNonlinPas],1);
        if bestGroupIdx==1
            checkMuscles = maxLinActInd;
        elseif bestGroupIdx==2
            checkMuscles = maxNonlinActInd;
        elseif bestGroupIdx==3
            checkMuscles = maxLinPasInd;
        elseif bestGroupIdx==4
            checkMuscles = maxNonlinPasInd;
        end
        bestMuscLinAct = reshape(permute(muscAvgNeurEval(n,3,checkMuscles),[3,2,1]),[1,numMuscles]);
        bestMuscNonlinAct = reshape(permute(muscAvgNeurEval(n,4,checkMuscles),[3,2,1]),[1,numMuscles]);
        bestMuscLinPas = reshape(permute(muscAvgNeurEval(n,5,checkMuscles),[3,2,1]),[1,numMuscles]);
        bestMuscNonlinPas = reshape(permute(muscAvgNeurEval(n,6,checkMuscles),[3,2,1]),[1,numMuscles]);
        
        linNonlinCondensedTable(n,:,:) = muscAvgNeurEval(n,:,checkMuscles);
        
        %save top muscles for unit
        for c=1:numel(checkMuscles)
            unitMuscles{c,n} = muscleArray{checkMuscles(1,c)};
        end
        
        %find mean of 5 best muscle group for each category
        meanLinAct = mean(bestMuscLinAct);
        meanNonlinAct = mean(bestMuscNonlinAct);
        meanLinPas = mean(bestMuscLinPas);
        meanNonlinPas = mean(bestMuscNonlinPas);

        %add this unit pR2 values to master lists
        linAct(1,n) = meanLinAct;
        nonlinAct(1,n) = meanNonlinAct;
        linPas(1,n) = meanLinPas;
        nonlinPas(1,n) = meanNonlinPas;
    end
    
    NonlinCondensedTable = mean(linNonlinCondensedTable,3);
    NonlinCondensedTable = table(NonlinCondensedTable,'VariableNames',columnNames);
    
    if doPlots
        %scatter plots
        low = -1;
        high = 1;
        step = 0.25;
        fig1 = figure('Name','Linear/Nonlinear Active/Passive comparison across units');
        sgtitle('Linear/Nonlinear Active/Passive comparison across units');
        x = [low:high];
        y = x;

        subplot(2,2,1);
        scatter(linAct,linPas)
        title('Linear active vs. passive');
        xlabel('Linear Active pR2')
        ylabel('Linear Passive pR2')
        xlim([low high])
        set(gca,'XTick',low:step:high)
        ylim([low high])
        set(gca,'YTick',low:step:high)
        hold on
        plot(x,y)

        subplot(2,2,2);
        scatter(nonlinAct,nonlinPas)
        title('Nonlinear active vs. passive');
        xlabel('Nonlinear Active pR2')
        ylabel('Nonlinear Passive pR2')
        xlim([low high])
        set(gca,'XTick',low:step:high)
        ylim([low high])
        set(gca,'YTick',low:step:high)
        hold on
        plot(x,y)

        subplot(2,2,3);
        scatter(linAct,nonlinAct)
        title('Active linear vs nonlinear');
        xlabel('Active Linear pR2')
        ylabel('Active Nonlinear pR2')
        xlim([low high])
        set(gca,'XTick',low:step:high)
        ylim([low high])
        set(gca,'YTick',low:step:high)
        hold on
        plot(x,y)

        subplot(2,2,4);
        scatter(linPas,nonlinPas)
        title('Passive linear vs nonlinear');
        xlabel('Passive Linear pR2')
        ylabel('Passive Nonlinear pR2')
        xlim([low high])
        set(gca,'XTick',low:step:high)
        ylim([low high])
        set(gca,'YTick',low:step:high)
        hold on
        plot(x,y)
    end
    
    linNonlinCondensedTable4d(:,:,:,d) = NonlinCondensedTable;
end
    
    clear maxLinAct maxLinActInd bestMuscLinAct
    clear maxNonlinAct maxNonlinActInd bestMuscNonlinAct
    clear maxLinPas maxLinPasInd bestMuscLinPas
    clear maxNonlinPas maxNonlinPasInd bestMuscNonlinPas
    clear linActInd linPasInd nonlinActInd nonlinPasInd
    clear bests j check indexes addition muscleLabels x
    clear n c
    clear meanLinAct meanLinPas meanNonlinAct meanNonlinPas
    clear meanAct meanPas meanLin meanNonlin
    clear checkMuscles muscleArray1