%************************************************************************
%************************************************************************
%************************************************************************

% This work is licensed under a Creative Commons Attribution-NonCommercial
% -ShareAlike 3.0 United States License" 
% https://creativecommons.org/licenses/by-nc-sa/3.0/us/
%
%************************************************************************
%************************************************************************
%************************************************************************
clc
clear all
close all
fclose all;
commandwindow()
%
% If you have a GPU, you should use it
%
GPU_TF=true;%Use GPU
% GPU_TF=false;%Don't Use GPU

% This is where you enter the path to your data file.
%
% The format of each data file is:
%
% Subject Session F00001  F00002  F00003  F00004  F00005  F00006  F00007  F00008  F00009  F00010
%    1       1     0.07    0.34   -0.07    0.91    0.04    1.12   -0.92    0.31    0.22   -1.07
%    2       1     0.27   -0.83    0.81   -1.73    0.83   -0.58   -0.71    0.70    0.01   -0.92
%    3       1     0.31    1.38   -1.20    1.42   -1.98   -0.24    0.60    1.22    0.03   -0.64
%    4       1     0.24   -0.43    0.06   -1.13   -0.07   -1.21    0.32    0.24    0.46   -0.52
%    5       1    -1.00   -0.96   -0.41   -0.79   -0.45    0.26    0.57   -1.13   -0.59   -0.79
%    1       2     0.48    1.21    0.15    0.89    0.49    1.03   -0.02    0.15    0.73   -0.33
%    2       2    -0.22   -0.33    0.52   -1.97    0.65   -1.11   -1.53    0.99    0.30   -0.34
%    3       2     0.17    1.43   -0.78    1.57   -1.86    0.07    0.72    1.00   -0.53   -0.62
%    4       2    -0.29    0.79    0.07   -0.55   -0.48   -1.44    1.01   -0.35    1.50   -0.81
%    5       2    -0.50   -0.79   -0.36   -0.69   -0.36   -0.12    0.37   -0.36   -0.81   -0.17
%
%    This is a dummy set with only 5 subjects and 10 features. Typically,
%    you would have many more subjects. Session (Occasion) is the visit
%    number of the data in the row. This can be only 1 or 2.
%    This program assumes that you want all the subjects and all the
%    features in the dataset.  If you want to choose a subset of subjects
%    and a subset of features, you will need to do that, save that reduced
%    data set and use that file here.  Subjects do not have to be in
%    sequential order, but every subject in Session 1 must also occur in Session
%    2.  It is expected that the first row contains header information.
%
% This is the path to my file, in my data directory. The file has 10
% features and 10,000 subjects.
%
FileName='Band_6_NFeat_010_NumberOfSubjects_00010000.csv';%10,000 Subjects
%
% Also included are:
% FileName='Band_6_NFeat_010_NumberOfSubjects_00001000.csv';%1,000 Subjects
% FileName='Band_6_NFeat_010_NumberOfSubjects_00100000.csv';%100,000 Subjects
%

PathToFileName=[pwd '\InputData\' FileName];
fprintf('\nAbout to read in feature file %s\n',PathToFileName)
tic
DataArray=csvread(PathToFileName,1,0);
fprintf('\nFinished reading in feature file, time in minutes = %0.3f\n',toc/60) 
%
% Sort data by session and subject
%
fprintf('\nSorting Data by Session and Subject\n')
DataArray=sortrows(DataArray,[2 1]);
% 
% Get Some Input File Information and Extract feature vectors
%
[rows,cols]=size(DataArray);
NumberOfSubjects=rows/2;
NumberOfFeatures=cols-2;
features = DataArray(:,3:NumberOfFeatures+2);
%
% If you have a data set with less than about 10,000 subjects you can do a
% traditional ROC analysis as well as the Binary Search algoithm.  This
% allows one to compare a traditional EER produced by this data set to that
% produces by the binary search method.  The traditional approach is a
% particularly fast and efficient method using function fastAUC.m
% https://www.mathworks.com/matlabcentral/fileexchange/42860-fast-auc-calculator-and-roc-curve-plotter
% Note that for a traditional approach, FAR and FNR are evaluated only for actual
% similarity scores.  For the Binary Search method, the FAR and FNR do no
% depend or actual similarity scores, so there may be small differences.

if NumberOfSubjects <= 10000
    %
    % L_S_Array, returned by QuickDistance, is a 2 column array with column
    % 1 filled with simiarlity scores (N=NumberOfSubjects^2), and column 2
    % is either a 1 (for genuine scores) or 0 (for impostor scores).
    %
    fprintf('\nGetting similarity scores\n')
    L_S_Array=QuickDistance(features);
    
    fprintf('\nComputing EER using FastAUC\n')
    [~,fpr,tpr] = fastAUC(logical(L_S_Array(:,2)),L_S_Array(:,1));
    fnr=1-tpr;
    clear tpr
    figure(1)
    plot(fpr,fnr)
    axis square
    hold on
    xlabel('False positive rate') 
    ylabel('False negative rate')
    %
    % Find the absolute value of the differences between fpr and fnr.
    %
    mydiff=abs(fpr-fnr);
    %
    % Find the minimum differences between fpr and fnr
    %
    min_diff=min(mydiff);
    %
    % Find the index to the minimum difference
    %
    Index=find(mydiff==min_diff);
    % 
    % For EER, take the mean of fpr and fnr at the minimum index
    % EER is an error rate, not a percent at this point.
    %
    EER_FastAUC=mean([fpr(Index) fnr(Index)]);
    plot([EER_FastAUC EER_FastAUC],[0 EER_FastAUC],':r')
    plot([0 EER_FastAUC],[EER_FastAUC EER_FastAUC],':r')
    %
    % Convert EER to percent
    %
    EER_FastAUC=EER_FastAUC*100;
    title({['Number of Subjects = ' num2str(NumberOfSubjects)],['EER using FastAUC Function = ' num2str(EER_FastAUC,'%0.4f')]})
    fprintf('\nEER using FastAUC Function = %0.4f\n',EER_FastAUC)

end
%
% Compute the total number of Impostor scores and Genuine scores
%
NumberImpostorSimScores=(NumberOfSubjects^2)-NumberOfSubjects;
NumberGenuineSimScores = NumberOfSubjects;
% 
% Set up main data arrays
%
session1=[DataArray(1:NumberOfSubjects,1) DataArray(1:NumberOfSubjects,3:NumberOfFeatures+2)];
session2=[DataArray(NumberOfSubjects+1:NumberOfSubjects*2,1) DataArray(NumberOfSubjects+1:NumberOfSubjects*2,3:NumberOfFeatures+2)];
clear DataArray

features1 = session1(:,2:end);
features2 = session2(:,2:end);
clear session1 session2
%
% Keep track of time to complete each loop and FAR and FRR for each Binary
% Search Loop
%
LoopTimes=NaN(1000,1);
FAR_Loop=NaN(1000,1);
FRR_Loop=NaN(1000,1);
%
% Start a timer to determine how long each loop takes
%
BinSearchStartTime=tic;
LoopCount=1;
% 
% Initial Lower threshold = 0.9. Upper - 1.0
%
left=0.0;
right=1.0;
%
% Get the size of each batch. The number of similarity scores computed
% for each batch is n_batch^2
%
if NumberOfSubjects == 1e5
    n_batch=5000;
elseif NumberOfSubjects == 1e4
    n_batch=1000;
elseif NumberOfSubjects == 1e3
    n_batch=100;
end
%
% Main Binary Search Loop
%
fprintf('\nStarting Binary Search Loop\n')
while true
    
    %
    % Start a loop timer
    %
    LoopStartTime=tic;
    %
    % The actual threshold to be evaluated for this loop is middle
    % threshold (mid).
    %
    mid = (left + right) / 2;
    %
    % Calculate the number of false acceptances and false rejections
    %
    [fa_sum,fr_sum]=calc_roc(features1, features2, n_batch, mid, GPU_TF);
    %
    % compute the false acceptance rate and the false rejection rate
    %
    far=fa_sum/NumberImpostorSimScores;
    frr=fr_sum/NumberGenuineSimScores;
    FAR_Loop(LoopCount+1)=far;
    FRR_Loop(LoopCount+1)=frr;

    fprintf('LoopCount=%d, Left = %f, Right = %f, far=%0.6f,frr=%0.6f\n', LoopCount,left,right,far,frr)
    %
    % Stopping condition. User can change this to 1e-5, 1e-6, etc... for
    % greater precision of the EER estimate
    %
    if abs(far-frr)<1e-4
         LoopTimes(LoopCount)=toc(LoopStartTime);
         break
    end
    %
    % if not stopping, update lower and upper thresholds and middle
    % threshold (mid).
    %
     if far > frr
        left = mid;
     end
     if far < frr       
        right = mid;
     end
     
     LoopTimes(LoopCount)=toc(LoopStartTime);
     LoopCount=LoopCount+1;
end
%
% End of Binary Search Loop
%
fprintf('\nFinished with Binary Search\n')
fprintf('\n******************* Final Binary Search Loop LoopReport *******************\n')
StoppingCode='1e-4';
[EER_Binary_Search]=LoopReport(LoopCount,far,frr,StoppingCode,FileName,BinSearchStartTime);

figure(2)
LoopTimes(isnan(LoopTimes))=[];
NLoops=length(LoopTimes);
plot(1:NLoops,LoopTimes,'or','MarkerSize',6,'MarkerFaceColor','r')
xlabel('Loop Count')
ylabel('Seconds per Loop')
MaxLoopTime=max(LoopTimes);
ymax=(MaxLoopTime*1.1);
ylim([0 ymax])
xlim([0 LoopCount+1])
set(gca,'XTick',0:NLoops+1)
title({['Number of Subjects = ' num2str(NumberOfSubjects)],'Seconds Per Loop'})
grid(gca,'on')
% grid(gca,'minor')
ax=gca;
ax.GridAlpha = 0.5;
ax.GridColor=[0 0 0];

figure(3)
plot(1:NLoops,FRR_Loop(1:NLoops),'-b');
hold on
plot(1:NLoops,FAR_Loop(1:NLoops),'-g');
xlabel('Loop Count')
ylabel('Error Rate')
legend('FRR','FAR','location','northeast')
title({['Number of Subjects = ' num2str(NumberOfSubjects) ', Error Rates Per Loop'],['EER from Binary Search = ' num2str(EER_Binary_Search,'%0.4f')]})
%
%This function arranges all the figures so that all are visible on one
%monitor. https://www.mathworks.com/matlabcentral/fileexchange/48480-automatically-arrange-figure-windows
%
autoArrangeFigures();

function [fa_sum,fr_sum]=calc_roc(features1, features2, n_batch, mid, GPU_TF)

    fa_sum=[];
    fr_sum=[];
    %
    %  split each session into batches and update stats for one batch pair at a time
    %
    Iteration=1;
    if size(features1,1) == 100000
        LoopReportIterations_TF=true;
    else
        LoopReportIterations_TF=false;
    end

    for i = 1:n_batch:size(features1,1)
        for j = 1:n_batch:size(features2,1)
            
            if LoopReportIterations_TF && mod(Iteration,50)==1
                fprintf('Iteration = %d\n',Iteration)
            end

            if GPU_TF
                TheseFeatures1=gpuArray(features1(i:i+n_batch-1,:));
                TheseFeatures2=gpuArray(features2(j:j+n_batch-1,:));
            else
                TheseFeatures1=features1(i:i+n_batch-1,:);
                TheseFeatures2=features2(j:j+n_batch-1,:);
            end                
                  
            scores = cosine_similarity(TheseFeatures1,TheseFeatures2);
            if i==j
                gen_scores=diag(scores);
                scores(1:1+size(scores,1):end) = NaN;
                scores=scores(:);
                scores(isnan(scores))=[];
                imp_scores=scores;
            else
                gen_scores=[];
                imp_scores=scores;
            end           
            [fa_sum,fr_sum]= ...
                Process_Scores(gen_scores(:),imp_scores(:), mid,fa_sum,fr_sum);
                     
            Iteration=Iteration+1;
        end
    end

end

function [fa_sum,fr_sum]= Process_Scores(gen_scores,imp_scores, mid, fa_sum, fr_sum)

    fa_these_scores=length(find(imp_scores >= mid));
    fr_these_scores=length(find(gen_scores < mid));
    if isempty(fa_sum)
        fa_sum=fa_these_scores;
    else
        fa_sum=fa_sum+fa_these_scores;
    end
    if isempty(fr_sum)
        fr_sum=fr_these_scores;
    else
        fr_sum=fr_sum+fr_these_scores;
    end
end
        
function scores=cosine_similarity(X, Y)
    D=pdist2(X,Y,'cosine');
    MinD=min(min(D));
    MaxD=max(max(D));
    ScaledD=(D-MinD)/(MaxD-MinD);
    scores=1-ScaledD;
end

function [EER]=LoopReport(LoopCount,far,frr,stopping_code,FileName,BinSearchStartTime)
    DurationOfBinarySearchLoop=toc(BinSearchStartTime);
    DurationOfBinarySearchLoop=DurationOfBinarySearchLoop/60.0;
    fprintf('\n\nThere were %d Iterations to reach abs(FAR-FRR) < %s\n',LoopCount,stopping_code)
    fprintf('The duration of the binary search process was %0.3f min\n',DurationOfBinarySearchLoop)
    fprintf('The time per iteration is %0.4f min\n',DurationOfBinarySearchLoop/LoopCount)
    EER=mean([far frr])*100;
    fprintf('For abs(FAR-FRR) < %s, found EER = %0.20f\n',stopping_code,EER(1));

    WhereIsExtension=strfind(FileName,'.csv');
    OutFileName=[pwd '\OutputFiles\' FileName(1:WhereIsExtension-1) '_BinSearchStats.csv'];
    warning('off')
    delete(OutFileName)
    warning('on')
    fid=fopen(OutFileName,'a');
    fprintf(fid,'FileName,Loop,Code,DurLoop(min),DurIter(min),EER\n');
    fprintf(fid,'%s,%d,%s,%0.3f,%0.4f,%0.10f\n',FileName(1:WhereIsExtension-1),LoopCount,stopping_code,DurationOfBinarySearchLoop,DurationOfBinarySearchLoop/LoopCount,EER);
    fclose(fid);
end
