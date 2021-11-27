%%
clear all
close all
clc

% get subject information
argindlg = inputdlg( ...
    {'Subject number (###)','Session (0/1/2)'}, ...
    'KETACONF',1,{'','','','',''});
if isempty(argindlg)
    error('Experiment cancelled!');
end
subj = str2num(argindlg{1});
sess = str2num(argindlg{2});

% run experiment
KETACONF_run_expe(subj,sess);