% MVPA decoding script. Pretty much does what brainVoyager does but enables
% easier data access/manipulation. Made to work with BV file outputs. 
% Analysis is done using the open source machine learnig library LIBSVM
% (https://www.csie.ntu.edu.tw/~cjlin/libsvm/). 
% Calls up VOI and leave one run out cross validation mvp files created
% from brainVoyager. So does require some GUI based work.
%
% Not the most efficient code, finding how to utilize VMPs will be better

%% Cleanup before working
clear all; clc;

%% Participant Setup
par         = {};
par{end+1}  = {'p01'}; % SHA
par{end+1}  = {'p02'}; % PSS
par{end+1}  = {'p03'}; % AJC
par{end+1}  = {'p04'}; % NKC
disp('start svm')
root_dir    = '/Volumes/macLab/01_brainVoyager/semDuo';

% empty matrices for data
bothHemMean     = [];
allMeanROI      = [];
parMean         = [];
% totalAcc    = [];
totalInfo   = [];
hem         = {'left','rght'};

%% experiment variables
allPars     = size(par,2);  % total number of participants
allRuns     = 8;            % total number of runs
allHems     = size(hem,2);  % only 2 hemispheres
allPairs    = 2;            % total # of pairs is 2 because comparing SR NR 
learn       = 'Learn';
test        = 'Test';
runName     = 'RunLevelSplit-'; % damn you bv output..

%% runs analysis for total number of participants
for n = 1:allPars

    par_dir     = fullfile(root_dir,par{n}{1}); % participant directory
    voi_dir     = fullfile(par_dir,'rois');    % voi/stats directory
    
    bothHemMean     = [];
    allMeanROI      = [];
    
    %% loop through hemisphere
    for h = 1:allHems;
        allMeanROI = [];
        cd(voi_dir); % cd to voi/stats dir
        
        % vind voi file for hemisphere
        if h == 1;
            voiFile     = dir('p*_all_lh.voi');       % find voi file
        elseif h == 2;
            voiFile     = dir('p*_all_rh.voi');       % find voi file
        end
        
        voi         = fullfile(voiFile.name);   % extract voi file name
        curVoi      = BVQXfile(voi);            % load voiFile using BVQXfile f(x)
        allVoi      = curVoi.NrOfVOIs;          % number of total VOIs in file
        
        %% loop through VOI
        for v = 1:allVoi;

            voiName             = curVoi(v).Name;
            avgMeanPairs        = [];
            
            %% loop through ROIs
            for pairs = 1:allPairs
                pairNo      = '10';
                totalInfo   = [];
                allRunAcc   = [];
                
                for r = 1:allRuns;
                    
                    % get name of each run file. Change accordingly to fit your output
                    runNumber   = num2str(r);       % current run number (in string)
                    runSearch   = strcat(runName,runNumber,'_'); % the name of file to search
                    if h == 1 && pairs == 1
                        pairName    = 'SR_RGHT_';
                    elseif h == 1 && pairs == 2
                        pairName    = 'NR_RGHT_';
                    elseif h == 2 && pairs == 1
                        pairName    = 'SR_LEFT_';
                    elseif h == 2 && pairs == 2
                        pairName    = 'NR_LEFT_';
                    end
                    mvpa_dir    = fullfile(par_dir,'stats','/mvpa/multiClass/');
                    cd(mvpa_dir)
                    
                    % specify, locate, and extract train mvp file names
                    trainFileName       = strcat(runSearch,learn,'_',voiName,'_',pairName,num2str(pairNo),'*')     % name of mvp file used for training
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
%                     opt.kerntype = {'linear'};
%                     opt.svmtype = {'c'}';
                    model               = svmtrain(classTrain,sfeat,'-t 0 -m 10000');            % train model, kernel = linear
                    [plab,acc,dvpe]     = svmpredict(classTest,featTest,model,'-q');             % test model
                    
                    %% Print data into matrix for easy extraction
                    fullInfo            = strcat(par{n},'_',runSearch,voiName);
                    totalInfo           = [totalInfo, fullInfo];
                    allRunAcc           = [allRunAcc, acc(1)]; % accuracy for all runs (leave one run out & test on others)
                    
                end
                avgMeanPairs = [avgMeanPairs mean(allRunAcc)]; % average for all 4 pairs
            end
            allMeanROI = [allMeanROI, avgMeanPairs]; % accuracy for all ROIs
        end
        
        bothHemMean = [bothHemMean allMeanROI]; % accuracy for all ROIs (left then right)
    end
    parMean = [parMean; bothHemMean]; % accuracy for all participants
end

