% Initialize/calculate delays for discrete-delay reactions
% (Delays are manually set - if reactions are rearranged such that tau indicies change, this code 
% must be reset)

% Create a cache of previous species values for transcription/translation delays.
RelAn_tau6   = RelAn;
IkBat_tau8   = IkBat;
cReln_tau28   = cReln;
IkBet_tau30   = IkBet;
RelAn_tau49   = RelAn;
if v.PHASE ~= 1
    % Match current timepoint into t vector - append, if necessary
    idx = find(DELAY.t(1:DELAY.idx)>=t,1,'first');
    if isempty(idx)
        idx = DELAY.idx+1;
    end
    DELAY.t(idx) = t;
    DELAY.idx = idx;
    DELAY.rela(DELAY.idx) = RelAn;
    DELAY.ikbat(DELAY.idx) = IkBat;
    DELAY.crel(DELAY.idx) = cReln;
    DELAY.ikbet(DELAY.idx) = IkBet;
    if DELAY.idx > length(DELAY.t)
        error(['Delay index out of range: increase memory size (line 13 - ',...
            'sz = ',num2str(length(DELAY.t)),';). Reached t = ',num2str(DELAY.t(end))])
    end

    % 1) IkBa mRNA transcription/processing delay (RXN 6)
    tau = p(6,4);
    if tau > 0  % IkBa inducible txn delay
        if t > tau
            RelAn_tau6 = ...
                interp1(DELAY.t(1:DELAY.idx-1),...
                        DELAY.rela(1:DELAY.idx-1),t-tau);
        else
            RelAn_tau6 = ...
                DELAY.rela(1); % 1st index is the steady-state value
        end
    end
    
    % 2) IkBa protein translation/maturation delay (RXN 8)
    tau = p(8,2);
    if tau > 0  % IkBb inducible txn delay
        if t > tau
            IkBat_tau8 = ...
                interp1(DELAY.t(1:DELAY.idx-1),...
                        DELAY.ikbat(1:DELAY.idx-1),t-tau);
        else
            IkBat_tau8 = ...
                DELAY.ikbat(1); % 1st index is the steady-state value
        end
    end
    
    % 3) IkBe mRNA transcription/processing delay (RXN 28)
    tau = p(28,4);
    if tau > 0  % IkBa inducible txn delay
        if t > tau
            cReln_tau28 = ...
                interp1(DELAY.t(1:DELAY.idx-1),...
                        DELAY.crel(1:DELAY.idx-1),t-tau);
        else
            cReln_tau28 = ...
                DELAY.crel(1); % 1st index is the steady-state value
        end
    end
    
    % 4) IkBe protein translation/maturation delay (RXN 8)
    tau = p(30,2);
    if tau > 0  % IkBb inducible txn delay
        if t > tau
            IkBet_tau30 = ...
                interp1(DELAY.t(1:DELAY.idx-1),...
                        DELAY.ikbet(1:DELAY.idx-1),t-tau);
        else
            IkBet_tau30 = ...
                DELAY.ikbet(1); % 1st index is the steady-state value
        end
    end
    
    % 5) IkBa mRNA transcription/processing delay (RXN 49)
    tau = p(49,4);
    if tau > 0  % IkBa inducible txn delay
        if t > tau
            RelAn_tau49 = ...
                interp1(DELAY.t(1:DELAY.idx-1),...
                        DELAY.rela(1:DELAY.idx-1),t-tau);
        else
            RelAn_tau49 = ...
                DELAY.rela(1); % 1st index is the steady-state value
        end
    end
  
end

% Pulsed stimulus (ode_options.PULSE_TIME) - if time value is greater than pulse length, set external stimuli to zero
if isfield(v,'PULSE_TIME') && (t > v.PULSE_TIME)
    TNF = 0;
    LPS = 0;
    polyIC = 0;
    Pam3CSK = 0;
    CpG = 0;
    FLA = 0;
    R848 = 0;
    stim=0;
end
