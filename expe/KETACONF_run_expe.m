function KETACONF_run_expe(subj,sess)

% clear command window
clc;

% set path
addpath('./Toolboxes/Rand');
addpath('./Toolboxes/IO');
addpath('./Toolboxes/Stimuli');

% initialize random number generator
seed = SetupRand;

% generate experiment structure
fprintf('Generating experiment structure... ');
if sess > 0, n = 6; else, n = 2; end
expe = KETACONF_gen_expe(n);
fprintf('Done.\n\n');

% store experiment information
expe.subj = subj; % subject number
expe.sess = sess; % session number
expe.date = datestr(now,'yyyymmdd-HHMM'); % date/time
expe.seed = seed; % random seed

% set stimulation parameters
syncflip = true; % synchronize video flips to screen refresh rate?
sendtrig = true; % send triggers to parallel port?
iscr     = max(Screen('Screens')); % screen index
maxres   = false; % use maximal screen resolution?
ppd      = 45; % screen pixels per degree of visual angle
fxtndmtr = deg2pix(24/60,2); % fixation point diameter
probdmtr = deg2pix(40/60,2); % response probe diameter
probwdth = deg2pix(4/60,1); % response probe width
pbarsize = deg2pix([16,2],1); % progress bar size (width/height)
pbarwdth = deg2pix(8/60,1); % progress bar width
lumibg   = 128/255; % background luminance
lumipbar = [112,160]/255; % progress bar luminance

% do not try to send triggers to parallel port on Mac!
if ismac
    sendtrig = false;
end

% set output folder name
foldname = sprintf('./Data/S%03d',subj);
if ~exist(foldname,'dir')
    mkdir(foldname);
end

% get experimental conditions
seqsig_lst = expe.cfg.seqsig_lst; % list of signal strengths
seqlen_lst = expe.cfg.seqlen_lst; % list of sequence lengths
lotval_lst = expe.cfg.lotval_lst; % list of lottery values

video   = [];
audio   = [];
aborted = false;
errmess = [];

try
    
    % hide cursor and stop spilling key presses into MATLAB windows
    HideCursor;
    FlushEvents;
    ListenChar(2);
    
    % check keyboard responsiveness before doing anything
    fprintf('Press any key to check keyboard responsiveness... ');
    if WaitKeyPress([],10) == 0
        fprintf('\n\n');
        error('No key press detected after 10 seconds.');
    else
        fprintf('Good.\n\n');
    end
    
    % define keys
    KbName('UnifyKeyNames');
    keywait = KbName('space');
    keyquit = KbName('ESCAPE');
    keyresp = KbName({'D','L'});
    
    % log Psychtoolbox startup information to disk
    flogname = sprintf('KETACONF_S%03d_sess%d_log.txt',expe.subj,expe.sess);
    diary(sprintf('%s/%s',foldname,flogname));
    
    % open main window
    % set screen synchronization properties
    % see 'help SyncTrouble',
    %     'help BeampositionQueries' or
    %     'help ConserveVRAMSettings' for more information
    if syncflip
        if maxres
            % set screen to highest/maximal resolution
            r = Screen('Resolutions',iscr);
            Screen('Resolution',iscr,max([r.width]),max([r.height]));
        end
        if ispc
            % soften synchronization test requirements
            Screen('Preference','SyncTestSettings',[],[],0.2,10);
            % enforce beamposition workaround for missing VBL interval
            Screen('Preference','ConserveVRAM',bitor(4096,Screen('Preference','ConserveVRAM')));
        end
        Screen('Preference','VisualDebuglevel',3);
    else
        % skip synchronization tests altogether, for debugging purposes only!
        Screen('Preference','SkipSyncTests',1);
        Screen('Preference','VisualDebuglevel',0);
        Screen('Preference','SuppressAllWarnings',1);
    end
    % set font properties
    if ismac
        txtfnt = 'Helvetica';
        txtsiz = deg2pix(1,1);
    elseif ispc
        txtfnt = 'Arial';
        txtsiz = deg2pix(0.5,1);
    end
    Screen('Preference','DefaultFontName',txtfnt);
    Screen('Preference','DefaultFontSize',txtsiz);
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask','General','UseFastOffscreenWindows');
    PsychImaging('AddTask','General','NormalizedHighresColorRange');
    video.i = iscr;
    video.res = Screen('Resolution',video.i);
    video.h = PsychImaging('OpenWindow',video.i,0);
    [video.x,video.y] = Screen('WindowSize',video.h);
    video.ifi = Screen('GetFlipInterval',video.h,100,50e-6,10);
    Screen('BlendFunction',video.h,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    Priority(MaxPriority(video.h));
    Screen('ColorRange',video.h,1);
    Screen('FillRect',video.h,lumibg);
    Screen('Flip',video.h);
    
    % store screen parameters
    expe.scr.ppd = ppd;
    expe.scr.res = [video.x,video.y];
    expe.scr.ifi = video.ifi;
    
    % open audio device
    InitializePsychSound(1);
    devices = PsychPortAudio('GetDevices');
    % select output device(s)
    isout = cellfun(@(n)n == 2,{devices.NrOutputChannels});
    if ismac
        % select non-AirPlay device(s)
        isval = cellfun(@(s)~strcmp(s,'AirPlay'),{devices.DeviceName});
    elseif ispc
        % select ASIO device(s)
        isval = cellfun(@(s)strcmp(s,'ASIO'),{devices.HostAudioAPIName});
    else
        error('OS not supported!');
    end
    i = find(isout & isval);
    if isempty(i)
        error('no valid audio device found!');
    elseif length(i) > 1
        % select device with lowest output latency
        [~,imin] = min([devices(i).HighOutputLatency]);
        i = i(imin);
    end
    audio.i = devices(i).DeviceIndex;
    audio.freq = devices(i).DefaultSampleRate;
    audio.h = PsychPortAudio('Open',audio.i,1,1,audio.freq,2);
    PsychPortAudio('RunMode',audio.h,1);
    PsychPortAudio('Volume',audio.h,0.9);
    
    % make feedback sounds
    [tonebuf,tonepos] = CreateAudioBuffer( ...
        CreateDoubleTone([660,880],[0.083,0.167],0,audio.freq), ...
        CreateDoubleTone([660,440],[0.083,0.167],0,audio.freq), ...
        CreateDoubleTone([440,440],[0.083,0.167],0,audio.freq));
    PsychPortAudio('FillBuffer',audio.h,tonebuf);
    
    if sendtrig
        config_io;
        write_lpt(0);
    end
    
    % stop logging for now
    diary('off');
    
    % make card frame texture
    cfg          = [];
    cfg.ppd      = ppd;
    cfg.cardtype = 0;
    img = make_card(cfg);
    framtex = Screen('MakeTexture',video.h,img,[],[],2);
    cardrec = CenterRectOnPoint(Screen('Rect',framtex),video.x/2,video.y/2);
    
    % create fixation point texture
    img = cat(3,zeros(fxtndmtr),CreateCircularAperture(fxtndmtr));
    fxtntex = Screen('MakeTexture',video.h,img,[],[],2);
    fxtnrec = CenterRectOnPoint(Screen('Rect',fxtntex),video.x/2,video.y/2);
    
    % create response probe texture
    img = cat(3,zeros(probdmtr),CreateCircle(probdmtr,probwdth));
    probtex = Screen('MakeTexture',video.h,img,[],[],2);
    probrec = CenterRectOnPoint(Screen('Rect',probtex),video.x/2,video.y/2);
    
    % create progress bar rectangle
    pbarrec = CenterRectOnPoint([0,0,pbarsize],video.x/2,video.y/2);
    
    t0 = Screen('Flip',video.h);
    
    % create result/event structure
    rslt = [];
    evnt = [];
    
    % get block structure
    blck = expe.blck;
    nblck = length(blck);
    
    % set prior for sd_cue
    sd_cue_hat = 0.5; % initial sd_cue estimate
    sd_val = 0.10:0.01:2.00; % sd_cue llf axis
    sd_llf = log(normpdf(log2(sd_val),log2(sd_cue_hat),1)); % prior sd_cue llf
    sd_llf = sd_llf-max(sd_llf); % normalize sd_cue llf
    
    update_ini = false;
    for iblck = 1:nblck
        
        % unpack block structure
        condtn = blck(iblck).condtn;
        seqcat = blck(iblck).seqcat;
        seqsig = blck(iblck).seqsig;
        seqlen = blck(iblck).seqlen;
        lotval = blck(iblck).lotval;
        
        % initialize event structure
        evnt(iblck).iseq = [];
        evnt(iblck).time = [];
        evnt(iblck).type = {};
        
        nseq = length(condtn); % number of sequences in block
        
        resp  = zeros(1,nseq); % categorization responses (1:orange or 2:blue)
        tresp = zeros(1,nseq); % categorization response times (s)
        conf  = zeros(1,nseq); % confirmation responses (1:confirm or 2:defect)
        tconf = zeros(1,nseq); % confirmation response times (s)
        rwrd  = zeros(1,nseq); % rewards

        if update_ini
            % update sd_cue estimate w.r.t. 2nd half of previous block
            cfg        = [];
            cfg.mu_cue = blck(iblck-1).mu_cue(37:end);
            cfg.seqcat = blck(iblck-1).seqcat(37:end);
            cfg.seqlen = blck(iblck-1).seqlen(37:end);
            cfg.resp   = rslt(iblck-1).resp(37:end);
            cfg.sd_val = sd_val; % prior sd_cue llf axis
            cfg.sd_llf = sd_llf; % prior sd_cue llf
            out = update_sd_cue(cfg);
            sd_cue_hat = out.sd_cue; % posterior sd_cue estimate
            sd_llf = out.sd_llf; % posterior sd_cue llf
        end
        
        % store titration step information
        titr(1).sd_cue = sd_cue_hat;
        titr(1).sd_val = sd_val;
        titr(1).sd_llf = sd_llf;
        
        % configure sequence generation
        cfg            = [];
        cfg.sd_cue     = sd_cue_hat;
        cfg.pcor       = mean(lotval_lst);
        cfg.seqsig_lst = seqsig_lst;
        cfg.seqlen_lst = seqlen_lst;
        cfg.seqcat     = seqcat(1:36);
        cfg.seqsig     = seqsig(1:36);
        cfg.seqlen     = seqlen(1:36);
        % generate sequences of angles for titration step
        out = gen_seqang(cfg);
        seqang = out.seqang;
        catdir = out.catdir;
        mu_cue = out.mu_cue;
        
        % short break before block start
        tonset = t0+roundfp(2.400,0.240);
        labeltxt = sprintf('Appuyez sur [espace] pour démarrer le bloc %d/%d.',iblck,nblck);
        labelrec = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,deg2pix(1,1));
        Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
        Screen('DrawingFinished',video.h);
        Screen('Flip',video.h,tonset);
        WaitKeyPress(keywait);
        Screen('DrawTexture',video.h,framtex,[],cardrec);
        Screen('DrawTexture',video.h,fxtntex,[],fxtnrec);
        Screen('DrawingFinished',video.h);
        t0 = Screen('Flip',video.h);
        if sendtrig
            write_lpt(255);
            WaitSecs(0.010);
            write_lpt(0);
        end
        tstart = t0; % block start time
        tonset = t0+roundfp(4.000,0.400);

        for iseq = 1:nseq
            
            % check if abort key is pressed
            if CheckKeyPress(keyquit)
                aborted = true;
                break
            end
            
            if iseq == 37
                
                % update sd_cue estimate w.r.t. 1st half of block
                cfg        = [];
                cfg.mu_cue = mu_cue(1:36);
                cfg.seqcat = seqcat(1:36);
                cfg.seqlen = seqlen(1:36);
                cfg.resp   = resp(1:36);
                cfg.sd_val = sd_val; % prior sd_cue llf axis
                cfg.sd_llf = sd_llf; % prior sd_cue llf
                out = update_sd_cue(cfg);
                sd_cue_hat = out.sd_cue; % posterior sd_cue estimate
                sd_llf = out.sd_llf; % posterior sd_cue llf
                
                % store titration step information
                titr(2).sd_cue = sd_cue_hat;
                titr(2).sd_val = sd_val;
                titr(2).sd_llf = sd_llf;
                
                % configure sequence generation
                cfg            = [];
                cfg.sd_cue     = sd_cue_hat;
                cfg.pcor       = mean(lotval_lst);
                cfg.seqsig_lst = seqsig_lst;
                cfg.seqlen_lst = seqlen_lst;
                cfg.seqcat     = seqcat(37:end);
                cfg.seqsig     = seqsig(37:end);
                cfg.seqlen     = seqlen(37:end);
                % generate sequences of angles for titration step
                out = gen_seqang(cfg);
                seqang = cat(2,seqang,out.seqang);
                catdir = cat(2,catdir,out.catdir);
                mu_cue = cat(2,mu_cue,out.mu_cue);
                
                if iblck == 1
                    % short break before block start
                    tonset = t0+roundfp(2.400,0.240);
                    labeltxt = 'Appuyez sur [espace] pour reprendre.';
                    labelrec = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,deg2pix(1,1));
                    Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
                    Screen('DrawingFinished',video.h);
                    Screen('Flip',video.h,tonset);
                    WaitKeyPress(keywait);
                    Screen('DrawTexture',video.h,framtex,[],cardrec);
                    Screen('DrawTexture',video.h,fxtntex,[],fxtnrec);
                    Screen('DrawingFinished',video.h);
                    t0 = Screen('Flip',video.h);
                    tonset = t0+roundfp(4.000,0.400);
                end
                
            end
            
            responded = false; % responded?
            confirmed = false; % confirmed?
            
            % make card textures
            cfg      = [];
            cfg.ppd  = ppd;
            cfg.axis = catdir(iseq);
            cardtex = cell(1,2);
            cfg.cardtype = 1;
            img = make_card(cfg);
            cardtex{1} = Screen('MakeTexture',video.h,img,[],[],2);
            ncards = seqlen(iseq);
            cardtex{2} = zeros(ncards,1);
            for i = 1:ncards
                cfg.cardtype = 3;
                cfg.angl = seqang{iseq}(i);
                img = make_card(cfg);
                cardtex{2}(i) = Screen('MakeTexture',video.h,img,[],[],2);
            end
            % make lottery texture
            cfg       = [];
            cfg.ppd   = ppd;
            cfg.pfill = lotval(iseq);
            img = make_pie(cfg);
            pietex = Screen('MakeTexture',video.h,img,[],[],2);
            
            % draw mapping cue
            Screen('DrawTexture',video.h,cardtex{1},[],cardrec);
            Screen('DrawTexture',video.h,fxtntex,[],fxtnrec);
            Screen('DrawingFinished',video.h);
            t = Screen('Flip',video.h,tonset);
            if sendtrig
                write_lpt(170);
                WaitSecs(0.010);
                write_lpt(0);
            end
            flag_evnt(t,'CUE');
            tonset = t+roundfp(0.800,0.080);
            
            % draw cards
            for i = 1:ncards
                Screen('DrawTexture',video.h,cardtex{2}(i),[],cardrec);
                Screen('DrawTexture',video.h,fxtntex,[],fxtnrec);
                Screen('DrawingFinished',video.h);
                while GetSecs < tonset
                    % check for early key presses
                    if CheckKeyPress(keyresp) > 0
                        flag_evnt(GetSecs,'RESP_EARLY');
                        break
                    end
                end
                t = Screen('Flip',video.h,tonset);
                if sendtrig
                    write_lpt(15);
                    WaitSecs(0.010);
                    write_lpt(0);
                end
                flag_evnt(t,'CARD');
                toffset = t+roundfp(0.100);
                tonset = t+roundfp(0.400,0.040);
                Screen('DrawTexture',video.h,cardtex{1},[],cardrec);
                Screen('DrawTexture',video.h,fxtntex,[],fxtnrec);
                Screen('DrawingFinished',video.h);
                while GetSecs < toffset
                    % check for early key presses
                    if CheckKeyPress(keyresp) > 0
                        flag_evnt(GetSecs,'RESP_EARLY');
                        break
                    end
                end
                Screen('Flip',video.h,toffset);
            end

            % draw response probe
            Screen('DrawTexture',video.h,cardtex{1},[],cardrec);
            Screen('DrawTexture',video.h,fxtntex,[],fxtnrec);
            Screen('DrawTexture',video.h,probtex,[],probrec);
            Screen('DrawingFinished',video.h);
            t = Screen('Flip',video.h,tonset);
            if sendtrig
                write_lpt(63);
                WaitSecs(0.010);
                write_lpt(0);
            end
            flag_evnt(t,'PROBE');
            tgo = t;
            toffset = t+roundfp(2.000);
            % wait for response
            while true
                [key,t] = CheckKeyPress(keyresp);
                if t > toffset
                    resp(iseq) = -1;
                    rt(iseq) = 0;
                    responded = false;
                    break
                end
                if key > 0
                    resp(iseq) = ceil(key);
                    rt(iseq) = t-tgo;
                    responded = true;
                    if sendtrig
                        write_lpt(42);
                        WaitSecs(0.010);
                        write_lpt(0);
                    end
                    flag_evnt(t,'RESP');
                    break
                end
            end
            Screen('DrawTexture',video.h,cardtex{1},[],cardrec);
            Screen('DrawTexture',video.h,fxtntex,[],fxtnrec);
            Screen('DrawingFinished',video.h);
            t = Screen('Flip',video.h);
            
            if condtn(iseq) == 0 || ~responded
                % play feedback
                tonset = t+roundfp(0.400,0.040);
                if ~responded
                    itone = 3;
                else
                    itone = 1+(resp(iseq) ~= seqcat(iseq));
                end
                PsychPortAudio('SetLoop',audio.h,tonepos(itone,1),tonepos(itone,2));
                PsychPortAudio('Start',audio.h,1,tonset,0);
                Screen('DrawTexture',video.h,framtex,[],cardrec);
                Screen('DrawTexture',video.h,fxtntex,[],fxtnrec);
                Screen('DrawingFinished',video.h);
                t = Screen('Flip',video.h,tonset);
                t0 = t; % end of trial
                tonset = t0+roundfp(1.200,0.120);
            else
                % wait for lottery
                if condtn(iseq) == 1 % short delay
                    tonset = t+roundfp(0.400,0.040);
                else % long delay
                    tonset = t+roundfp(4.000,0.400);
                end
                % draw lottery offer
                Screen('DrawTexture',video.h,pietex,[],cardrec);
                Screen('DrawTexture',video.h,fxtntex,[],fxtnrec);
                Screen('DrawingFinished',video.h);
                while GetSecs < tonset
                    % check for early key presses
                    if CheckKeyPress(keyresp) > 0
                        flag_evnt(GetSecs,'RESP_DELAY');
                        break
                    end
                end
                t = Screen('Flip',video.h,tonset);
                if sendtrig
                    write_lpt(47);
                    WaitSecs(0.010);
                    write_lpt(0);
                end
                flag_evnt(t,'OFFER');
                tgo = t;
                toffset = t+roundfp(2.000);
                % wait for response
                while true
                    [key,t] = CheckKeyPress(keyresp);
                    if t > toffset
                        conf(iseq) = -1;
                        tconf(iseq) = 0;
                        confirmed = false;
                        break
                    end
                    if key > 0
                        conf(iseq) = 1+(ceil(key) ~= resp(iseq)); 
                        tconf(iseq) = t-tgo;
                        confirmed = true;
                        if sendtrig
                            write_lpt(58);
                            WaitSecs(0.010);
                            write_lpt(0);
                        end
                        flag_evnt(t,'CONF');
                        break
                    end
                end
                Screen('DrawTexture',video.h,pietex,[],cardrec);
                Screen('DrawTexture',video.h,fxtntex,[],fxtnrec);
                Screen('DrawingFinished',video.h);
                t = Screen('Flip',video.h);
                tonset = t+roundfp(0.400,0.040);
                if ~confirmed
                    % play feedback
                    itone = 3;
                    PsychPortAudio('SetLoop',audio.h,tonepos(itone,1),tonepos(itone,2));
                    PsychPortAudio('Start',audio.h,1,tonset,0);
                end
                Screen('DrawTexture',video.h,framtex,[],cardrec);
                Screen('DrawTexture',video.h,fxtntex,[],fxtnrec);
                Screen('DrawingFinished',video.h);
                t = Screen('Flip',video.h,tonset);
                t0 = t; % end of trial
                tonset = t0+roundfp(1.200,0.120);
            end

            % close textures
            Screen('Close',cardtex{1});
            Screen('Close',cardtex{2});
            Screen('Close',pietex);
            
            % assign reward
            if resp(iseq) < 0 || conf(iseq) < 0 % timeout
                rwrd(iseq) = 0;
            elseif conf(iseq) == 1 % confirmed
                rwrd(iseq) = resp(iseq) == seqcat(iseq);
            else % defected
                rwrd(iseq) = lotval(iseq);
            end
            
        end
        
        % add titration-dependent information to block structure
        blck(iblck).seqang = seqang; % sequences of angles (rad)
        blck(iblck).catdir = catdir; % categorization directions (rad)
        blck(iblck).mu_cue = mu_cue; % mean evidence per cue
        blck(iblck).titr   = titr; % titration step information
        
        % store results
        rslt(iblck).resp  = resp;
        rslt(iblck).tresp = tresp;
        rslt(iblck).conf  = conf;
        rslt(iblck).tconf = tconf;
        rslt(iblck).rwrd  = rwrd;
        
        % save results
        fpath = foldname;
        fname = sprintf('KETACONF_S%03d_sess%d_blck%d_%s.mat',expe.subj,expe.sess,iblck,expe.date);
        save(sprintf('%s/%s',fpath,fname),'expe','blck','rslt','evnt');
        
        if aborted
            break
        end
        
        % draw end-of-block screen
        if iblck > 1
            % draw progress bar
            draw_progress(round(sum(rwrd)),nseq);
            labeltxt = 'votre score';
            labelrec = Screen('TextBounds',video.h,labeltxt);
            labelrec = CenterRectOnPoint(labelrec,video.x/2,video.y/2-2*ppd);
            labelrec = AlignRect(labelrec,pbarrec,'left');
            Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
            % draw proportion of confirmation among responded trials
            pconf = mean(conf(conf > 0) == 1);
            labeltxt = sprintf('%2.0f / %2.0f',100*pconf,100*(1-pconf));
            labelrec = Screen('TextBounds',video.h,labeltxt);
            labelrec = CenterRectOnPoint(labelrec,video.x/2,video.y/2+2*ppd);
            labelrec = AlignRect(labelrec,pbarrec,'left');
            Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
        end
        labeltxt = sprintf('Fin du bloc %d/%d. Appuyez sur [espace].',iblck,nblck);
        labelrec = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,deg2pix(1,1));
        Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
        Screen('DrawingFinished',video.h);
        t = Screen('Flip',video.h,tonset);
        WaitKeyPress(keywait);
        t0 = Screen('Flip',video.h);
        
        % update sd_cue at the beginning of next block
        update_ini = true;
        
    end
    
    % save results before exiting
    fsufx = '';
    if aborted
        fsufx = '_aborted';
    end
    fpath = foldname;
    fname = sprintf('KETACONF_S%03d_sess%d_%s%s.mat',expe.subj,expe.sess,expe.date,fsufx);
    save(sprintf('%s/%s',fpath,fname),'expe','blck','rslt','evnt');
    
    % resume logging
    diary('on');
    
    % close video
    Priority(0);
    Screen('CloseAll');
    FlushEvents;
    ListenChar(0);
    ShowCursor;
    
    % close audio
    PsychPortAudio('Stop',audio.h,1);
    PsychPortAudio('Close');
    
    % stop logging
    diary('off');
    
catch
    
    % save results before crashing
    if exist('rslt','var') && ~isempty(rslt)
        fsufx = '_crashed';
        fpath = foldname;
        fname = sprintf('KETACONF_S%03d_sess%d_%s%s.mat',expe.subj,expe.sess,expe.date,fsufx);
        save(sprintf('%s/%s',fpath,fname),'expe','blck','rslt','evnt');
    end
    
    % resume logging
    diary('on');
    
    % close video
    if ~isempty(video)
        Priority(0);
        Screen('CloseAll');
        FlushEvents;
        ListenChar(0);
        ShowCursor;
    end
    
    % close audio
    if ~isempty(audio)
        PsychPortAudio('Close');
    end
    
    % stop logging
    diary('off');

    % handle error
    if nargout > 1
        errmess = lasterror;
        errmess = rmfield(errmess,'stack');
    else
        rethrow(lasterror);
    end
    
end

    function [p] = deg2pix(d,b)
        p = d*ppd;
        if nargin > 1 && b > 0
            p = max(round(p/b),1)*b;
        end
    end

    function [t] = roundfp(t,dt)
        if nargin < 2
            dt = 0;
        end
        n = (t+dt*[-1,+1])/video.ifi;
        n = round(n);
        n = n(1)-1+ceil((diff(n)+1)*rand);
        t = (n-0.5)*video.ifi;
    end

    function flag_evnt(etime,etype)
        evnt(iblck).iseq(end+1) = iseq;
        evnt(iblck).time(end+1) = etime-tstart;
        evnt(iblck).type{end+1} = etype;
    end

    function draw_progress(ncur,nmax)
        Screen('FillRect',video.h,lumipbar(2),pbarrec);
        w = round(RectWidth(pbarrec)*ncur/nmax);
        h = RectHeight(pbarrec);
        fbarrec = [0,0,w,h];
        fbarrec = AlignRect(fbarrec,pbarrec,'left');
        fbarrec = AlignRect(fbarrec,pbarrec,'top');
        Screen('FillRect',video.h,lumipbar(1),fbarrec);
        Screen('FrameRect',video.h,0,pbarrec,pbarwdth);
        x = video.x/2+round(RectWidth(pbarrec)*0.25);
        y = video.y/2+h*[-1,+1];
        Screen('DrawLine',video.h,0,x,y(1),x,y(2),probwdth);
    end

end