% MVPA decoding script. Pretty much does what brainVoyager does but enables
% batch processing, enabling easier data access/manipulation. Made to work 
% with BV file outputs. Analysis is done % using the open source machine 
% learning library LIBSVM (https://www.csie.ntu.edu.tw/~cjlin/libsvm/). 
% Calls up VOI and leave one run out cross validation mvp files created
% from brainVoyager. So does require some GUI based work. Future code might
% change that..but this should work for now. Oh, also loops for total # of 
% participants
%
% edit [2016.10.19] - uses vmp files to analyze, making future analysis
% easier

%% Cleanup before working
clear all; clc;
disp('start svm')

%% Participant Setup
par         = {};
par{end+1}  = {'p01'}; % SHA
par{end+1}  = {'p02'}; % PSS
par{end+1}  = {'p03'}; % AJC
par{end+1}  = {'p04'}; % NKC
allPars     = size(par,2);
root_dir    = '/Volumes/macLab/01_brainVoyager/semDuo'; % root directory
runSize     = 48; % # of trials

% matrices of zeros for labeling trials
allLabels   = [];% zeros(48,3,8);
allData     = [];
%% run for all participants
for n = 1%:allPars
    %% grabs labels for conditions you want to decode for all runs
    runs = [];
    for runIndex = 1:8
        tempAllLabels   = [];
        SR_RGHT_10 = [];
        SR_RGHT_20 = [];
        SR_RGHT_30 = [];
        SR_RGHT_40 = [];
        runN    = strcat('r',num2str(runIndex));
        par_dir = fullfile(root_dir,par{n}{1});
        run_dir = fullfile(par_dir,runN);
        cd(run_dir)
        vmpFile = strcat(runN,'_SCCTBL_3DMCTS_THPGLMF2c_TAL_16-Conditions_GLM-2G_PreOn-1-PostCond-6_z-t_Trials.vmp');
        runs(runIndex).run = BVQXfile(vmpFile);
        
        %% checks label of trial and relabel in numbers
        % labels = [condition hem type]
        for trialIndex = 1:runSize;
            
            vmpName = runs(runIndex).run.Map(1,trialIndex).Name;
            vmpData = runs(runIndex).run.Map(1,trialIndex).VMPData;
            
            nameChk = vmpName(1:10);
            if nameChk == 'SR_RIGHT_1'
                labels = [1 1 1];
                SR_RGHT_10 = [SR_RGHT_10; trialIndex];
            elseif nameChk == 'SR_RIGHT_2'
                labels = [1 1 2];
                SR_RGHT_20 = [SR_RGHT_20; trialIndex];
            elseif nameChk == 'SR_RIGHT_3'
                labels = [1 1 3];
                SR_RGHT_30 = [SR_RGHT_30; trialIndex];
            elseif nameChk == 'SR_RIGHT_4'
                labels = [1 1 4];
                SR_RGHT_40 = [SR_RGHT_40; trialIndex];
            elseif nameChk == 'NR_RIGHT_1'
                labels = [2 1 1];
            elseif nameChk == 'NR_RIGHT_2'
                labels = [2 1 2];
            elseif nameChk == 'NR_RIGHT_3'
                labels = [2 1 2];
            elseif nameChk == 'NR_RIGHT_4'
                labels = [2 1 4];
                
            elseif nameChk == 'SR_LEFT_10'
                labels = [1 2 1];
            elseif nameChk == 'SR_LEFT_20'
                labels = [1 2 2];
            elseif nameChk == 'SR_LEFT_30'
                labels = [1 2 3];
            elseif nameChk == 'SR_LEFT_40'
                labels = [1 2 4];
            elseif nameChk == 'NR_LEFT_10'
                labels = [2 2 1];
            elseif nameChk == 'NR_LEFT_20'
                labels = [2 2 2];
            elseif nameChk == 'NR_LEFT_30'
                labels = [2 2 3];   
            elseif nameChk == 'NR_LEFT_40'
                labels = [2 2 4];
            end
            
            tempAllLabels = [tempAllLabels; labels];
        end
        allLabels(:,:,runIndex) = tempAllLabels;
    end
    
    %% read in VOI
    voi_dir = fullfile(par_dir,'rois'); % voi/stats directory
    cd(voi_dir);
    voi     = BVQXfile('allHem.voi'); % need to compile voi file with all VOIs order: left -> right
    myVOIs  = voi.VOI;      % 
    voiN    = voi.NrOfVOIs; % # of VOIs
    
    %% run for all ROIs
    for voiIndex = 1:voiN
        voiVoxel                        = myVOIs(voiIndex).Voxels;  % Get the x,y,z coordinates of all the voxels
        voiBetas.Clust(voiIndex).Name   = myVOIs(voiIndex).Name;    % Set the VOIBetas.Clust(j) name as the given voi name.
        [voiVoxelSize, columns]         = size(voiVoxel);           % # of voxels in the clust; x, y, z coordinates
        
        for pairIndex = 1%:4
            classValuesRun = [];
            
            for runIndex = 1:8
                
                for trialIndex = 1:size(allLabels,1) % runs through all trials
                    curLabels = allLabels(:,:,runIndex);
                    
                    if voiIndex < voiN/2 % for the leftHem ROIs
                        %% extract relevant labels and designate correct labels. This is where changes might have to be made for different decoding comparisons
                        % if just comparing 2, labels: 1, -1 [check: 2016.10.20]
                        % if multi, just in numbers
                        
                        if curLabels(trialIndex,:) == [1 1 pairIndex];      % SR_RIGHT_*0
                            classValues = 1;
                        elseif curLabels(trialIndex,:) == [2 1 pairIndex];  % NR_RIGHT_*0
                            classValues = -1;
                        end
                        
                    elseif voiIndex >= voiN/2 % for the rghtHem ROIs
                        if curLabels(trialIndex,:) == [1 2 pairIndex];      % SR_LEFT*0
                            classValues = 1;
                        elseif curLabels(trialIndex,:) == [2 2 pairIndex];  % NR_LEFT_*0
                            classValues = -1;
                        end
                    end
                    classValuesRun = [classValuesRun; classValues];
                end
                d
            end
        end
    end
end
d


% empty matrices for data
bothHemMean     = [];
allMeanROI      = [];
parMean         = [];
% totalAcc    = [];
totalInfo   = [];
hem         = {'right','left'};

  % total number of participants
allRuns     = 8;            % total number of runs
allHems     = size(hem,2);  % only 2 hemispheres

%% runs analysis for all participants
for n = 1:allPars
%     avgMean     = [];
    par_dir     = fullfile(root_dir,par{n}{1}); % participant directory
    voi_dir     = fullfile(par_dir,'stats');    % voi/stats directory
    d
    bothHemMean     = [];
    allMeanROI      = [];
    
    %% loop through hemisphere
    for h = 1:allHems;
        allMeanROI = [];
        cd(voi_dir); % cd to voi/stats dir
        % vind voi file for hemisphere
        if h == 1;
            voiFile     = dir('evc_rh_combDV.voi');       % find voi file
        elseif h == 2;
            voiFile     = dir('evc_lh_combDV.voi');       % find voi file
        end
        
        voi         = fullfile(voiFile.name);   % extract voi file name
        curVoi      = BVQXfile(voi);            % load voiFile using BVQXfile f(x)
        allVoi      = curVoi.NrOfVOIs;          % number of total VOIs in file
        
        %% loop through VOI
        for v = 1:allVoi;
            
            learn       = 'Learn';
            test        = 'Test';
            
            voiName     = curVoi(v).Name;
            avgMeanPairs     = [];
            %% loop through ROIs
            for pairs = 1:4
                pairNo = pairs*10;
                totalInfo   = [];
                allRunAcc    = [];
                for r = 1:allRuns;
                    
                    % get name of each run file. Change accordingly to fit your output
                    runName     = 'RunLevelSplit-'; % damn you bv output..
                    runNumber   = num2str(r);       % current run number (in string)
                    runSearch   = strcat(runName,runNumber,'_'); % the name of file to search
                    if h ==1
                        mvpa_dir    = fullfile(par_dir,'stats','/mvpa/EVC_combDV/SR_left');    % voi/stats directory
                        pairName    = 'SR_LEFT_';
                    elseif h == 2
                        mvpa_dir    = fullfile(par_dir,'stats','/mvpa/EVC_combDV/SR_right');    % voi/stats directory
                        pairName    = 'SR_RIGHT_';
                    end
                    cd(mvpa_dir)
                    % specify, locate, and extract train mvp file names
                    trainFileName       = strcat(runSearch,learn,'_',voiName,'_',pairName,num2str(pairNo),'*');     % name of mvp file used for training
                    mvpTrainFile        = dir(trainFileName);                       % find the mvp file
                    mvpTrainFileName    = fullfile(mvpTrainFile.name);              % the name of mvp file
                    bvTrainMVP          = BVQXfile(mvpTrainFileName);               % load mvpFile using BVQXfile f(x)
                    classTrain          = bvTrainMVP.ClassValues;                   % get class labels (conditions)
                    featTrain           = bvTrainMVP.FeatureValues;                 % get the value of features (voxels)
                    sfeat               = sparse(featTrain);                        % sparse apart matrix "help sparse" for more info
                    
                    % specify, locate, and extract test mvp file names
                    testFileName        = strcat(runSearch,test,'_',voiName,'_',pairName,num2str(pairNo),'*');   % name of mvp file used for training
                    mvpTestFile         = dir(testFileName);                        % find the mvp file
                    mvpTestFileName     = fullfile(mvpTestFile.name);               % the name of mvp file
                    bvTestMVP           = BVQXfile(mvpTestFileName);                % load mvpFile using BVQXfile f(x)
                    classTest           = bvTestMVP.ClassValues;                    % get class labels (conditions)
                    featTest            = bvTestMVP.FeatureValues;                  % get the value of features (voxels)
                    sfeatT              = sparse(featTest);                         % sparse apart matrix "help sparse" for more info
                    
                    % The SVM part of the code. Utilizes LIBSVM as does BV
                    model               = svmtrain(classTrain,sfeat,'-t 0 -m 10000');            % train model, kernel = linear
                    [plab,acc,dvpe]     = svmpredict(classTest,featTest,model,'-q');             % test model
                    
                    %% Print data into matrix for easy extraction
                    fullInfo            = strcat(par{n},'_',runSearch,voiName);
                    totalInfo           = [totalInfo, fullInfo];
                    allRunAcc            = [allRunAcc, acc(1)]; % accuracy for all runs (leave one run out & test on others)

                end
                avgMeanPairs = [avgMeanPairs mean(allRunAcc)]; % average for all 4 pairs
            end
            allMeanROI = [allMeanROI, avgMeanPairs]; % accuracy for all ROIs
        end
        
        bothHemMean = [bothHemMean allMeanROI]; % accuracy for all ROIs (left then right)
    end
    parMean = [parMean; bothHemMean]; % accuracy for all participants
end


mean(allRunAcc(1,:))
