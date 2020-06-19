function topMuscBarPlots(muscAvgNeurEval,params)

% bar plots of pR2 values for each muscle
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
    assignParams(who,params); %overwrite params
    
    %remove dontWant muscles
    fullMuscleArray(dontWant) = [];
    muscleArray = fullMuscleArray;
    
    %loop through each unit
    k=0;
    ind=0;
    fig1 = figure('Name','Top muscles pR2 for units 1-15');
    fig2 = figure('Name','Top muscles pR2 for units 16-30');
    fig3 = figure('Name','Top muscles pR2 for units 31-45');
    fig4 = figure('Name','Top muscles pR2 for units 46-53');
    for n=1:numel(muscAvgNeurEval(:,1,1))
        k=k+1;
        muscleLabels = strings(1,numMuscles*4);
        if k<=15
            set(0,'CurrentFigure',fig1)
            subplot(3,5,k);
        elseif k<=30
            ind = k-15;
            set(0,'CurrentFigure',fig2)
            subplot(3,5,ind);
        elseif k<=45
            ind = k-30;
            set(0,'CurrentFigure',fig3)
            subplot(3,5,ind);
        else 
            ind = k-45;
            set(0,'CurrentFigure',fig4)
            subplot(2,5,ind);
        end
        
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
        
        %%FOR PLOTTING PR2 OF 5 TOP MUSCLES FOR EACH 4 CATEGORIES
        %order from least to greatest pR2
        [bestMuscLinAct,linActInd] = sort(bestMuscLinAct,'ascend');
        maxLinActInd = maxLinActInd(linActInd);
        [bestMuscNonlinAct,nonlinActInd] = sort(bestMuscNonlinAct,'ascend');
        maxNonlinActInd = maxNonlinActInd(nonlinActInd);
        [bestMuscLinPas,linPasInd] = sort(bestMuscLinPas,'ascend');
        maxLinPasInd = maxLinPasInd(linPasInd);
        [bestMuscNonlinPas,nonlinPasInd] = sort(bestMuscNonlinPas,'ascend');
        maxNonlinPasInd = maxNonlinPasInd(nonlinPasInd);
        y = [bestMuscLinAct bestMuscNonlinAct bestMuscLinPas bestMuscNonlinPas];
        
        %get muscle names of best muscles
        bests = [maxLinActInd maxNonlinActInd maxLinPasInd maxNonlinPasInd];
        for j=1:numel(bests)
            check = strncmpi(muscleLabels,muscleArray{bests(j)},numel(muscleArray{bests(j)}));
            indexes = find(check);
            if any(check)
                addition = strcat(muscleLabels{indexes(numel(indexes))},'.');
            else
                addition = muscleArray{bests(j)};
            end
            muscleLabels{j} = addition;
        end
        x = categorical(muscleLabels);
        x = reordercats(x,muscleLabels);
        
        %plot active linear (dark black)
        b1 = barh(x(1:numMuscles),y(1:numMuscles),'FaceColor','k','facealpha',1);
        
        hold on
        % get the current tick labels
        ticklabels = get(gca,'YTickLabel');
        % prepend a color for each tick label
        ticklabels_new = cell(size(ticklabels));
        for i = 1:length(ticklabels)
            if i>numMuscles && i<=numMuscles*2
                ticklabels_new{i} = ['\color{red} ' ticklabels{i}];
            elseif i>numMuscles*3 && i<=numMuscles*4
                ticklabels_new{i} = ['\color{red} ' ticklabels{i}];
            else
                ticklabels_new{i} = ['\color{black} ' ticklabels{i}];
            end
        end
        % set the tick labels
        set(gca, 'YTickLabel', ticklabels_new);
        
        %plot passive linear (dark red)
        b2 = barh(x(numMuscles+1:2*numMuscles),y(numMuscles+1:2*numMuscles),'FaceColor','r','facealpha',1);
        %plot active nonlinear (light black)
        b3 = barh(x(1+2*numMuscles:3*numMuscles),y(1+2*numMuscles:3*numMuscles),'FaceColor','k','facealpha',0.3);
        %plot passive nonlinear (light red)
        b4 = barh(x(1+3*numMuscles:4*numMuscles),y(1+3*numMuscles:4*numMuscles),'FaceColor','r','facealpha',0.3);
        
        %axes and legend
        xlim([-.5 0.5])
        set(gca,'XTick',-.5:0.25:0.5)
        title(strcat('Unit ',string(n)));
        sgtitle('Top muscles pR2 per unit');
        if k==1 || ind==1
            legend([b4,b3,b2,b1],'Passive, Nonlinear','Active, Nonlinear','Passive, Linear','Active, Linear',...
                'Position',[0 0.8 0.12 0.1])
        end
        suplabel('pR2');
    end
    
    clear k ind  maxLinAct maxLinActInd bestMuscLinAct
    clear maxNonlinAct maxNonlinActInd bestMuscNonlinAct
    clear maxLinPas maxLinPasInd bestMuscLinPas
    clear maxNonlinPas maxNonlinPasInd bestMuscNonlinPas
    clear bests j check indexes addition muscleLabels x
    clear n
