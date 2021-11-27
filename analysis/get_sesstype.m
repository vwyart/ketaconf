function [subj,sess] = get_sesstype(stype,sboth)
%  GET_SESSTYPE  Get session type for KETACONF study
%
%  Usage: [subj,sess] = GET_SESSTYPE(stype,sboth)
%
%  where stype corresponds to the session type (1=placebo, 2=ketamine)
%        sboth is set to true to keep only subjects with both sessions
%
%  The function outputs subj, corresponding to subject indices, and sess,
%  corresponding to the session number.
%
%  Valentin Wyart <valentin.wyart@inserm.fr>

% check input arguments
if nargin < 2
    sboth = false;
end
if nargin < 1
    error('Missing session type!');
end
if ~ismember(stype,[1,2])
    error('Invalid session type!');
end

% 24 subjects including 3 pilots (S101, S103, S104) => 21 test subjects
% ketamine session aborted prematurely for 2 subjects (S112, S115) => 19 remaining subjects
% session type => 0:missing 1:placebo 2:ketamine
sesstype = [ ...
    0 0; ... % S101 => pilot
    1 2; ... % S102
    0 0; ... % S103 => pilot
    0 0; ... % S104 => pilot
    2 1; ... % S105
    1 2; ... % S106
    1 2; ... % S107
    2 1; ... % S108
    1 2; ... % S109
    2 1; ... % S110
    2 1; ... % S111
    1 0; ... % S112 => ketamine session aborted prematurely
    1 2; ... % S113
    1 2; ... % S114
    0 1; ... % S115 => ketamine session aborted prematurely
    2 0; ... % S116 => EEG data missing for placebo session
    1 2; ... % S117
    1 2; ... % S118
    2 1; ... % S119
    2 1; ... % S120
    2 1; ... % S121
    2 1; ... % S122
    1 2; ... % S123
    1 2];    % S124

% get corresponding subject and session numbers
[subj,sess] = ind2sub(size(sesstype),find(sesstype(:) == stype));
[subj,i] = sort(subj);
sess = sess(i);

if sboth
    % exclude subjects for which the other session is missing
    isout = ismember(subj,find(any(sesstype == 0,2)));
    subj(isout) = [];
    sess(isout) = [];
end