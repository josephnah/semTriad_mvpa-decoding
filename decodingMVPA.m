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
% par{end+1}  = {'p02'}; % PSS
% par{end+1}  = {'p03'}; % AJC
% par{end+1}  = {'p04'}; % NKC
allPars     = size(par,2);
root_dir    = '/Volumes/macLab/01_brainVoyager/semDuo'; % root directory
vmpName     = '_SCCTBL_3DMCTS_THPGLMF2c_TAL_16-Conditions_GLM-2G_PreOn-1-PostCond-6_z-t_Trials.vmp';
trialN      = 48; % # of trials

objPairN    = {'10'; '20'; '30'; '40'}; 
pairSide    = {'RGHT_' 'LEFT_'};    % len: 5
pairCond    = {'SR_' 'NR_'};        % len: 3

runN        = 8;
pars = [];

vmpTrainData = [];
vmpTestData  = [];

trainTestRuns = [1 2 3 4 5 6 7 8];

%% load all VMPs for all participants and all runs
for n = 1:allPars
    for runIndex = 1:8 % loops through all runs
        runNumb = strcat('r',num2str(runIndex));
        par_dir = fullfile(root_dir,par{n}{1});
        run_dir = fullfile(par_dir,runNumb);
        cd(run_dir)
        vmpFile = strcat(runNumb,vmpName);
        pars(n).runs(runIndex).vmp = BVQXfile(vmpFile);
    end
end

%% run analysis for all participants
for n = 1:allPars
    for pairIndex = 1:size(objPairN,1) % loops through pairs (10, 20, 30, 40)
        
        classIndex   = [];
        for trialIndex = 1:trialN % loops through 48 trials in 1 run
            classLabl = pars(n).runs(1).vmp.Map(1,trialIndex).Name;
            % find appropriate trials based on hemisphere
            if pairIndex < size(objPairN,1)/2      % for Left Hem ROIs %%%%%%%%%%%%%%%% WRONG
                curLabl = strcat(pairSide(1),objPairN(pairIndex));
            elseif pairIndex >= size(objPairN,1)/2 % for right Hem ROIs
                curLabl = strcat(pairSide(2),objPairN(pairIndex));
            end
            
            % gets trial # of relevant classes for train/test matrix
            if strcmp(curLabl,classLabl(4:10))
                classIndex = [classIndex; trialIndex];
            end
        end
                        
        % loop through runs to grab train/test data
        for runIndex = 1:runN       % loops through runs (1 ~ 8), runIndex = test run
            
            testRun     = runIndex;
            trainRun    = setdiff(trainTestRuns,testRun);  % selects the testing runs: n = 7
            
            % loop through test runs and grab relevant data
            for testDataRun = 1:size(classIndex,1);
                vmpTestData = pars(n).runs(testRun).vmp.Map(1,classIndex).VMPData;
            end
            trainLablIndex = [];            
        end
        d
        %% read in VOMs
        vom_dir = fullfile(par_dir,'vom'); % voi/stats directory
        cd(vom_dir);
        
        vom     = dir('p0*.vom'); % finds all VOM files
        vomN    = size(vom,1); % # of VOIs (should be 14)
        
        %% run for all ROIs
        for vomIndex = 1:vomN
            curVOM                          = BVQXfile(vom(vomIndex).name); % grabs 1 VOM 
            vomVoxel                        = curVOM.VOM.Voxels;            % Get the x,y,z coordinates of all the voxels
            vomData.Clust(vomIndex).Name    = curVOM.VOM.Name;              % Set the VOIBetas.Clust(j) name as the given voi name.
            [voiVoxelSize, columns]         = size(vomVoxel);               % # of voxels in the clust; x, y, z coordinates
            talVoxel                        = [];                           % empty matrix for data
            nonTalVoxel                     = [];
            for voxel = 1:voiVoxelSize
                x = vomVoxel(voxel,1);
                y = vomVoxel(voxel,2);
                z = vomVoxel(voxel,3);
                
                % convert the coordinates of the voxels (x,y,z) to
                % matlab array indices (Mx, My, Mz)using Tal2Matlab
                [Mx,My,Mz]  = Tal2Matlab(x,y,z);
                nonTalVoxel = [nonTalVoxel; x y z];
                talVoxel    = [talVoxel; Mx My Mz];
            end
        end

        d
    end

    
        
        for pairIndex = 1%:4
            classValuesRun = [];
            
            for runIndex = 1:8
                
                for trialIndex = 1:size(allLabels,1) % runs through all trials
                    curLabels = allLabels(:,:,runIndex);
                    
                    if vomIndex < voiN/2 % for the leftHem ROIs
                        %% extract relevant labels and designate correct labels. This is where changes might have to be made for different decoding comparisons
                        % if just comparing 2, labels: 1, -1 [check: 2016.10.20]
                        % if multi, just in numbers
                        
                        if curLabels(trialIndex,:) == [1 1 pairIndex];      % SR_RIGHT_*0
                            classValues = 1;
                        elseif curLabels(trialIndex,:) == [2 1 pairIndex];  % NR_RIGHT_*0
                            classValues = -1;
                        end
                        
                    elseif vomIndex >= voiN/2 % for the rghtHem ROIs
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
    vom_dir     = fullfile(par_dir,'stats');    % voi/stats directory
    d
    bothHemMean     = [];
    allMeanROI      = [];
    
    %% loop through hemisphere
    for h = 1:allHems;
        allMeanROI = [];
        cd(vom_dir); % cd to voi/stats dir
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

% Input VOI and vmp files from trial estimation

rootDir = '/Users/joecool890/Dropbox/MATLAB/GWU/imaging/00_MVPA/02_SVM/fromFatma/';
cd(rootDir);

Experiment=[];
Experiment(1).run=BVQXfile('JL_r1_SCCAI2_3DMCTS_THPGLMF2c_TAL_Trials.vmp');
Experiment(2).run=BVQXfile('JL_r2_SCCAI2_3DMCTS_THPGLMF2c_TAL_Trials.vmp');
Experiment(3).run=BVQXfile('JL_r3_SCCAI2_3DMCTS_THPGLMF2c_TAL_Trials.vmp');
Experiment(4).run=BVQXfile('JL_r4_SCCAI2_3DMCTS_THPGLMF2c_TAL_Trials.vmp');


%now read in voi data
myVOI=BVQXfile('JL_LH_complete.voi');
myVOIs = myVOI.VOI;
nVOI=myVOI.NrOfVOIs

block=1;
L=1; % block and L are the same numberS)
for j = 1:1 %analysis for 1 ROI
    [A,B]=size(myVOIs(j).Voxels);
    coordinations= myVOIs(j).Voxels;
    %store coordinates, name and size for each voi
    Alldata=zeros(48,A);
    Alllabels=zeros(48,1);
    for r = 1:4
        %loops through all runs
        for b=1:12
            %loops thru all trials
            data=Experiment(r).run.Map(1,b).VMPData;
            %check label
            name=Experiment(r).run.Map(1,b).Name;
            k=name(1:6);
            if k=='LeftFa'
                Alllabels(L,1)=1;
            elseif k=='RightF'
                Alllabels(L,1)=2;
            elseif k=='LeftHo'
                Alllabels(L,1)=3;
            elseif k=='RightH'
                Alllabels(L,1)=4;
            end
            
            
            for v = 1:A %loops through all voxels in a VOI
                x= coordinations(v,1);
                y= coordinations(v,2);
                z= coordinations(v,3);
                
                [Mx,My,Mz] = Tal2Matlab(x,y,z);
                Alldata(block, v)=data(Mx,My,Mz);
                
            end
            block=block+1;
            L=L+1;
            
            
        end
    end
    
    clear coordinations
    %end of voi-clear for next
end

T=1;
for i=1:4
    for k=1:4
        if k~=i
            for z=1:12
                
                trainlabels(T,:)=Alllabels((k-1)*12+z,:);
                traindata(T,:)=Alldata((k-1)*12+z,:);
                T=T+1;
            end
        end
    end
    opt.kerntype={'linear'};
    opt.ncross=1;
    opt.svmtype='C';
    svm.cost=2.71;
    model=svmtrain(trainlabels, traindata,opt)
    
    TT=1;
    for z=1:12
        
        testdata(TT,:)=Alldata((i-1)*12+z,:);
        testlabels(TT,:)=Alllabels((i-1)*12+z,:);
        TT=TT+1;
    end
    opt.ncross=1;
    [plab,acc,dvpe]=svmpredict(testlabels,testdata,model,opt)
    
    SVM(i).trainmat=traindata;
    SVM(i).trainlabel=trainlabels;
    SVM(i).testmat=testdata;
    SVM(i).testlabel=testlabels;
    SVM(i).predictions=plab;
    SVM(i).accuracy=acc;
    SVM(i).decision=dvpe;
    SVM(i).Model=model;
    
    clear testdata
    clear testlabels
    clear trainlabels
    clear traindata
    
end
