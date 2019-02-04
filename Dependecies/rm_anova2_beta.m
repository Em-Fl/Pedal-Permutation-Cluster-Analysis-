function stats = rm_anova2(Y,S,F1,F2,FACTNAMES)
%
% function stats = rm_anova2(Y,S,F1,F2,FACTNAMES)
%
% Two-factor, within-subject repeated measures ANOVA.
% For designs with two within-subject factors.
%
% Parameters:
%    Y          dependent variable (numeric) in a column vector
%    S          grouping variable for SUBJECT
%    F1         grouping variable for factor #1
%    F2         grouping variable for factor #2
%    FACTNAMES  a cell array w/ two char arrays: {'factor1', 'factor2'}
%
%    Y should be a 1-d column vector with all of your data (numeric).
%    The grouping variables should also be 1-d numeric, each with same
%    length as Y. Each entry in each of the grouping vectors indicates the
%    level # (or subject #) of the corresponding entry in Y.
%
% Returns:
%    stats is a cell array with the usual ANOVA table:
%      Source / ss / df / ms / F / p
%
% Notes:
%    Program does not do any input validation, so it is up to you to make
%    sure that you have passed in the parameters in the correct form:
%
%       Y, S, F1, and F2 must be numeric vectors all of the same length.
%
%       There must be at least one value in Y for each possible combination
%       of S, F1, and F2 (i.e. there must be at least one measurement per
%       subject per condition).
%
%       If there is more than one measurement per subject X condition, then
%       the program will take the mean of those measurements.
%
% Aaron Schurger (2005.02.04)
%   Derived from Keppel & Wickens (2004) "Design and Analysis" ch. 18
%

%
% Revision history...
%
% 11 December 2009 (Aaron Schurger)
% 
% Fixed error under "bracket terms"
% was: expY = sum(Y.^2);
% now: expY = sum(sum(sum(MEANS.^2)));
%

% AHORA,  Y es matriz abn x t

stats = cell(4,5);

F1_lvls = unique(F1);  %obtiene la lista de niveles de las variables
F2_lvls = unique(F2);
Subjs = unique(S);  

a = length(F1_lvls); % # of levels in factor 1
b = length(F2_lvls); % # of levels in factor 2
n = length(Subjs); % # of subjects
t = size(Y,2); % # de puntos de tiempo, no sé si lo uso para algo, pero ta


INDS = zeros(a,b,n); % this will hold arrays of indices
CELLS = zeros(a,b,n,t); % this will hold the data for each subject X condition
MEANS = zeros(a,b,n,t); % this will hold the means for each subj X condition

% Calculate means for each subject X condition.
% Keep data in CELLS, because in future we may want to allow options for
% how to compute the means (e.g. leaving out outliers > 3stdev, etc...).



for i=1:a % F1
    for j=1:b % F2
        for k=1:n % Subjs
            INDS(i,j,k) = find(F1==F1_lvls(i) & F2==F2_lvls(j) & S==Subjs(k));
            CELLS(i,j,k,:) = Y(INDS(i,j,k),:);
            MEANS(i,j,k,:) = mean(CELLS(i,j,k,:),2);
        end
    end
end

% make tables (see table 18.1, p. 402)
AB = reshape(sum(MEANS,3),a,b,t); % across subjects
AS = reshape(sum(MEANS,2),a,n,t); % across factor 2
BS = reshape(sum(MEANS,1),b,n,t); % across factor 1

A = reshape(sum(AB,2),a,t); % sum across columns, so result is axt column vector
B = reshape(sum(AB,1),b,t); % sum across rows, so result is bxt row vector
S = reshape(sum(AS,1),n,t); % sum across columns, so result is nxt row vector
T = sum(A); % could sum either A or B or S, choice is arbitrary

% degrees of freedom
dfA = a-1;
dfB = b-1;
dfAB = (a-1)*(b-1);
dfS = n-1;
dfAS = (a-1)*(n-1);
dfBS = (b-1)*(n-1);
dfABS = (a-1)*(b-1)*(n-1);

% bracket terms (expected value)
expA = sum(A.^2)./(b*n);
expB = sum(B.^2)./(a*n);
expAB = reshape(sum(sum(AB.^2))./n,1,t);
expS = sum(S.^2)./(a*b);
expAS = reshape(sum(sum(AS.^2))./b,1,t);
expBS = reshape(sum(sum(BS.^2))./a,1,t);
expY = reshape(sum(sum(sum(MEANS.^2))),1,t); %sum(Y.^2);
expT = T.^2 / (a*b*n);

% sums of squares
ssA = expA - expT;
ssB = expB - expT;
ssAB = expAB - expA - expB + expT;
ssS = expS - expT;
ssAS = expAS - expA - expS + expT;
ssBS = expBS - expB - expS + expT;
ssABS = expY - expAB - expAS - expBS + expA + expB + expS - expT;
ssTot = expY - expT;

% mean squares
msA = ssA / dfA;
msB = ssB / dfB;
msAB = ssAB / dfAB;
msS = ssS / dfS;
msAS = ssAS / dfAS;
msBS = ssBS / dfBS;
msABS = ssABS / dfABS;

% f statistic
fA = msA ./ msAS;
fB = msB ./ msBS;
fAB = msAB ./ msABS;

% p values
pA = 1-fcdf(fA,dfA,dfAS);
pB = 1-fcdf(fB,dfB,dfBS);
pAB = 1-fcdf(fAB,dfAB,dfABS);

% return values
stats = {'Source','SS','df','MS','F','p';...
         FACTNAMES{1}, ssA, dfA, msA, fA, pA;...
         FACTNAMES{2}, ssB, dfB, msB, fB, pB;...
         [FACTNAMES{1} ' x ' FACTNAMES{2}], ssAB, dfAB, msAB, fAB, pAB;...
         [FACTNAMES{1} ' x Subj'], ssAS, dfAS, msAS, [], [];...
         [FACTNAMES{2} ' x Subj'], ssBS, dfBS, msBS, [], [];...
         [FACTNAMES{1} ' x ' FACTNAMES{2} ' x Subj'], ssABS, dfABS, msABS, [], []};
 
 return