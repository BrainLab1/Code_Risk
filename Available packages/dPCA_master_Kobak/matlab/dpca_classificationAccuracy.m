% last update: 01.09.2015 -> varargout was added to the output to export distances between test and train data! (see bellow in the code)
%       bug fix: centeringTest was used for centering Xtest, previously it was centered using centeringTrain!   

function [accuracy, brier, varargout] = dpca_classificationAccuracy(Xfull, Xtrial, numOfTrials, varargin)

% accuracy = dpca_classificationAccuracy(X, Xtrial, numOfTrials) performs
% cross-validation to compute classification accuracy of the first several
% dPCA components. X is the data array. Xtrial is an array storing single
% trials. It has one extra dimension as compared with X and stores
% individual single trial firing rates, as opposed to the trial average.
% numOfTrials has one dimension fewer than X and for each neuron and
% combination of parameters (without time) specifies the number of
% available trials in X_trial. All entries have to be larger than 1.
%
% This function assumes that time parameter is stored in the last dimension of
% X. For datasets without time, some other cross-validation needs to be
% used.
%
% [accuracy, brier] = dpca_classificationAccuracy(X) returns Brier scores
% together with classification accuracy.
%
% [...] = dpca_classificationAccuracy(..., 'PARAM1',val1, 'PARAM2',val2, ...) 
% specifies optional parameter name/value pairs:
%
% 'numComps'        - number of components to test in each marginalization
%                     (default: 3)
%
% 'combinedParams'  - cell array of cell arrays specifying 
%                     which marginalizations should be added up together,
%                     e.g. for the three-parameter case with parameters
%                           1: stimulus
%                           2: decision
%                           3: time
%                     one could use the following value:
%                     {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}}.
%
%  'lambda'         - regularization parameter. It's going to be multiplied
%                     by the total variance of Xfull. Default value is
%                     zero. To use different lambdas for different
%                     marginalizations, provide an array instead of one
%                     number.
%
%  'numRep'         - number of cross-validation iterations. Default: 100
%
%  'decodingClasses'- specifies classes for each marginalization.
%                     E.g. for the three-parameter case with parameters
%                           1: stimulus  (3 values)
%                           2: decision  (2 values)
%                           3: time
%                     and combinedParams as specified above:
%                     {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}}
%                     one could use the following decodingClasses:
%                     {[1 1; 2 2; 3 3], [1 2; 1 2; 1 2], [], [1 2; 3 4; 5 6]}
%                     Default value is to use separate class for each
%                     condition, i.e.
%                     {[1 2; 3 4; 5 6], [1 2; 3 4; 5 6], [], [1 2; 3 4; 5 6]}
%
%  'verbose'        - If 'yes' (default) progress status is displayed. Set
%                     to 'no' to turn it off. Set to 'dots' for compact
%                     display.
%
%  'filename'       - If provided, accuracy and brier outputs will be saved
%                     in this file.
%
%  'timeSplits'     - an array of K integer numbers specifying time splits
%                     for time period splitting. All marginalizations will
%                     be additionally split into K+1 marginalizations,
%                     apart from the one corresponding to the last
%                     parameter (which is assumed to be time).
%
% 'timeParameter'   - is only used together with 'timeSplits', and must be
%                     provided. Specifies the time parameter. In the
%                     example above it is equal to 3.
%
% 'notToSplit'      - is only used together with 'timeSplits'. A cell array
%                     of cell arrays specifying which marginalizations
%                     should NOT be split. If not provided, all
%                     marginalizations will be split.

% ------------------- Bahareh 01.09.2015 ----------------------------------
% initialize the output
varargout = {[]};
% -------------------------------------------------------------------------

% default input parameters
options = struct('numComps',       3,                   ...  
                 'lambda',         0,                   ...
                 'numRep',         100,                 ...
                 'verbose',        'yes',               ...
                 'combinedParams', [],                  ...
                 'decodingClasses', [],                 ...
                 'timeSplits',      [],                 ...
                 'timeParameter',  [],                  ...
                 'notToSplit',     [],                  ...
                 'filename',       []);

% read input parameters
optionNames = fieldnames(options);
if mod(length(varargin),2) == 1
	error('Please provide propertyName/propertyValue pairs')
end
for pair = reshape(varargin,2,[])    % pair is {propName; propValue}
	if any(strcmp(pair{1}, optionNames))
        options.(pair{1}) = pair{2};
    else
        error('%s is not a recognized parameter name', pair{1})
	end
end

Xsum = bsxfun(@times, Xfull, numOfTrials);

% find time marginalization
timeComp = [];
if ~isempty(options.combinedParams)
    for k = 1:length(options.combinedParams)
        if options.combinedParams{k}{1} == length(size(Xfull))-1
            timeComp = k;
            break
        end
    end
else
    timeComp = ndims(Xfull)-1;
end

% set the number of components for the dPCA
if ~isempty(options.combinedParams)
    numCompsToUse = repmat(options.numComps, [1 length(options.combinedParams)]);
else
    numCompsToUse = repmat(options.numComps, [1 2^(ndims(Xfull)-1)-1]);
end
numCompsToUse(timeComp) = 0;

% set decoding classes
if isempty(options.decodingClasses)
    dims = size(Xfull);
    dims = dims(2:end-1);
    classes = reshape(1:prod(dims), dims);
    if ~isempty(options.combinedParams)
        numMarg = length(options.combinedParams);
    else
        numMarg = 2^(ndims(Xfull)-1)-1;
    end
    options.decodingClasses = {};
    for i=1:numMarg
        options.decodingClasses{i} = classes;
    end
end

    
dim = size(Xfull);
cln = {};
for i=1:length(dim)-2
    cln{i} = ':';
end

% compute dPCA on the full data
[~, Vfull, ~] = dpca(Xfull, numCompsToUse, ...
    'combinedParams', options.combinedParams, ...
    'lambda', options.lambda, ...
    'timeSplits', options.timeSplits, ...
    'timeParameter', options.timeParameter, ...
    'notToSplit', options.notToSplit, ...
    'order', 'no');

for rep = 1:options.numRep
    if strcmp(options.verbose, 'yes')
        display(['Repetition # ' num2str(rep) ' out of ' num2str(options.numRep)])
    elseif strcmp(options.verbose, 'dots')
        fprintf('.')
    end
    
    Xtest = dpca_getTestTrials(Xtrial, numOfTrials);
    Xtrain = bsxfun(@times, Xsum - Xtest, 1./(numOfTrials-1));
    
    [W,V,whichMarg] = dpca(Xtrain, numCompsToUse, ...
         'combinedParams', options.combinedParams, ...
         'lambda', options.lambda, ...
         'timeSplits', options.timeSplits, ...
         'timeParameter', options.timeParameter, ...
         'notToSplit', options.notToSplit, ...
         'order', 'no');
     
    % for debugging
    % dpca_plot(Xtrain, W, V, @dpca_plot_default, 'whichMarg', whichMarg);
    
    % rearranging the order of components so that fit mostly to the
    % components on the full data (probably not necessary at all)
    crosscorr = corr([V Vfull]);
    crosscorr = abs(crosscorr(size(V,2)+1:end, 1:size(V,2)));
    for marg = setdiff(unique(whichMarg), timeComp)
        d = find(whichMarg==marg);
        crc = crosscorr(d, d);
        order = [];
        for i=1:length(d)
            leftOver = setdiff(1:length(d), order);
            [~, num] = max(crc(leftOver,i));
            order(i) = leftOver(num);
        end
        [~, order] = max(crc);
        W(:,d) = W(:, d(order));
        V(:,d) = V(:, d(order));
    end
    
    Xtrain2D = Xtrain(:,:);
    centeringTrain = mean(Xtrain2D,2);
    Xtrain2D = bsxfun(@minus, Xtrain2D, centeringTrain);
    Xtest2D  = Xtest(:,:);
    centeringTest = mean(Xtest2D,2);
    Xtest2D  = bsxfun(@minus, Xtest2D,  centeringTest);
    
    % skip time marginalization
    for marg = setdiff(unique(whichMarg), timeComp)
        dd = find(whichMarg==marg);
        for d = 1:length(dd)
            % now we want to test linear decoder W(:,d) at each time point
            % t and see if it decodes classes in this marginalization
            % correctly, on the test pseudo-ensemble
            
            XtrainComp = reshape(W(:,dd(d))'*Xtrain2D, dim(2:end));
            XtestComp  = reshape(W(:,dd(d))'*Xtest2D,  dim(2:end));
            
            for t=1:dim(end)
                c = {cln{:}, t};
                trainClasses = XtrainComp(c{:});
                trainClasses = trainClasses(:);
                trainClassMeans = accumarray(options.decodingClasses{marg}(:), trainClasses, [], @mean);
            
                testValues = XtestComp(c{:});
                testValues = testValues(:);
                
                dist = bsxfun(@minus, testValues, trainClassMeans');
                [~, classification] = min(abs(dist), [], 2);
% % % % % ------------------- Bahareh 01.09.2015 ----------------------------------
% % % %                 dist2{marg,d,t,rep} = bsxfun(@minus, testValues, trainClassMeans');
% % % %                 [~, classification2{marg,d,t,rep}] = min(abs(dist), [], 2);
% % % % % -------------------------------------------------------------------------

                correctClassifications(marg, d, t, rep) = length(find(classification == options.decodingClasses{marg}(:)));
                
                prob = 1./abs(dist);
                prob = bsxfun(@times, prob, 1./sum(prob,2));
                actual = zeros(size(prob));
                actual(sub2ind(size(actual), (1:size(actual,1))', options.decodingClasses{marg}(:))) = 1;
                
                brier(marg, d, t, rep) = sum(sum((prob-actual).^2)) / length(unique(options.decodingClasses{marg}));
            end
        end
    end
end

% % % % % ------------------- Bahareh 01.09.2015 ----------------------------------
% % % % % varargout(1) = {dist2};
% % % % % varargout(2) = {classification2};
% % % % % clear dist2 classification2
% % % % for m=setdiff(unique(whichMarg), timeComp)
% % % %     for t=1:size(dist2,3)
% % % %         for rep=1:size(dist2,4)
% % % %             temp=zeros(size(dist2{m,1,t,rep}));
% % % %             for c=1:size(dist2,2)
% % % %                 temp = temp + (dist2{m,c,t,rep}.^2);
% % % %             end
% % % %             temp = sqrt(temp);
% % % %             [~,classResult{m,t,rep}] = min(abs(temp), [], 2);
% % % %             correctClassResult{m,t,rep} = length(find(classResult{m,t,rep} == options.decodingClasses{m}(:)));
% % % %             aa(marg, 1, t, rep) = correctClassResult{m,t,rep} ;
% % % %         end
% % % %     end
% % % % end
% % % % clear temp m t rep c
% % % % % -------------------------------------------------------------------------

for i=1:length(options.decodingClasses)
    numClasses(i) = length(unique(options.decodingClasses{i}));
end
numConditions = prod(dim(2:end-1));

accuracy = bsxfun(@times, sum(correctClassifications, 4), 1./(options.numRep*numConditions));
brier = bsxfun(@times, sum(brier,4), 1./(options.numRep*numConditions));


% % % % % ------------------- Bahareh 01.09.2015 ----------------------------------
% % % % accuracy2 = bsxfun(@times, sum(aa, 4), 1./(options.numRep*numConditions));
% % % % accuracy2(timeComp,:,:) = NaN;
% % % % varargout = {accuracy2};
% % % % % -------------------------------------------------------------------------

accuracy(timeComp,:,:) = NaN;
brier(timeComp,:,:) = NaN;

if strcmp(options.verbose, 'yes')
    rows = 1:size(accuracy,1);
    rows(timeComp) = NaN;
    rows(rows>timeComp) = rows(rows>timeComp) - 1;
    
    figure
    if 0  % bahareh 20160131
    for i=setdiff(1:size(accuracy,1), timeComp)
        for j=1:size(accuracy,2)
            subplot(length(rows)-1,3,(rows(i)-1)*3+j)
            title(['Marginalization #' num2str(i)])
            
            axis([1 size(accuracy,3) 0 1])
            hold on
            
            if ~isempty(options.filename) && exist(options.filename, 'file')
                vars = whos('-file', options.filename);
                if ismember('accuracyShuffle', {vars.name})
                    load(options.filename, 'accuracyShuffle')
                    maxSh = max(accuracyShuffle(i,:,:),[],3);
                    minSh = min(accuracyShuffle(i,:,:),[],3);
                    time = 1:length(maxSh);
                    h = patch([time fliplr(time)], [maxSh fliplr(minSh)], 'b');
                    set(h, 'FaceAlpha', 0.5)
                    set(h, 'EdgeColor', 'none')
                end
            end
            
            plot(xlim, 1/numClasses(i)*[1 1], 'k')
            plot(squeeze(accuracy(i,j,:)))
            
            plot(squeeze(brier(i,j,:)), 'r')        
        end
    end
    
    end %------ Bahareh
    
end

if ~isempty(options.filename)
    if exist(options.filename, 'file')
        save(options.filename, 'accuracy', 'brier', '-append')
    else
        save(options.filename, 'accuracy', 'brier')
    end
end