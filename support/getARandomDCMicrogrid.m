function net = getARandomDCMicrogrid(N)
%GETARANDOMDCMICROGRID Build the evolving random DCMG used by main.mlx.
%   net = getARandomDCMicrogrid(N) creates N heterogeneous DGs, constructs
%   a connected proximity-based physical network, solves its equilibrium,
%   and assigns the current main-workflow initial condition.

    areaSide = 8;
    pos = areaSide*(rand(N,2)-0.5);            % random 2D positions
    
    % Reasonable DG parameter ranges
    Lvals = 1e-3 * (0.9 + 0.2*rand(N,1));      % 0.9–1.1 mH
    Cvals = 1e-3 * (0.9 + 0.2*rand(N,1));      % 0.9–1.1 mF
    Rfvals= 0.05 + 0.1*rand(N,1);              % 0.05–0.15 Ohm
    vref  = 48*ones(N,1);                      % 48 V nominal
    
    % --- DGs with per-node ratings ---
    %% ---------------- DGs with fixed Vr=48 and intentional ZI-load imbalance ----------------
    Vr_fixed = 48;                     % rated bus voltage for all DGs
    pctHeavy = 0.6;          % 0.6         % fraction of "heavy" nodes (forces network flows)
    idxHeavy = randperm(N, ceil(pctHeavy*N));
    
    % Power rating bands (W) for light/heavy converters (affects Irated)
    Pr_light = [100 200];              % 2.5–5 A @ 48V
    Pr_heavy = [200 300]; %500,600             % 5.4–8.8 A @ 48V
    % REVIEW PROPOSAL MAIN-02 (CHECK): Ir_k below is 2*Pr_k/Vr_k, so the
    % actual current bands are twice the comments above. Decide whether the
    % factor of two is a deliberate converter headroom factor.
    
    % Impedance safety bounds (Ohm)
    DGs = repmat(DG(1,struct()), 1, N);
    for k = 1:N
        Vr_k = Vr_fixed;                % keep 48 V for every DG
    
        % Converter rating (Irated) via light/heavy bands
        if ismember(k, idxHeavy)
            Pr_k = randi(Pr_heavy);     % heavier converter capability
        else
            Pr_k = randi(Pr_light);
        end
        Ir_k = 2*Pr_k / Vr_k;             % A (converter rated current)
    
        % --- Z/I split for the **load** at this node (not converter rating) ---
        % Start with a Z share alphaZ ~ N(0.6,0.2), clamp, then compute P_Z, P_I.
        alphaZ = 0.5;
        % REVIEW PROPOSAL MAIN-03 (ADD): use this bounded random alternative
        % if heterogeneous Z/I splits are desired without negative loads.
        % alphaZ = min(0.9,max(0.1,0.6+0.2*randn()));
        % Add per-node jitter to **load** power magnitude around converter rating to create mismatch
        Pload_k = Pr_k * (0.7 + 0.60*rand());   % ~Pr_k ±15%, biased +10% to induce imports
        PZ_k = alphaZ * Pload_k;              % Z-part power (W)
        PI_k = Pload_k - PZ_k;                 % I-part power (W)
        
        % % --- Z/I split for the **load** at this node (not converter rating) ---
        % % Start with a Z share alphaZ ~ N(0.6,0.2), clamp, then compute P_Z, P_I.
        % alphaZ = 0.6 + 0.2*randn();
        % % Add per-node jitter to **load** power magnitude around converter rating to create mismatch
        % Pload_k = Pr_k * (1.10 + 0.30*(2*rand()-1));   % ~Pr_k ±15%, biased +10% to induce imports
        % PZ_k = alphaZ * Pload_k;              % Z-part power (W)
        % PI_k = Pload_k - PZ_k;                 % I-part power (W)
    
        % Convert Z power to impedance and I power to constant current
        RL_k   = Vr_k^2 / PZ_k;   % Ohm
        Ibar_k = (PI_k / Vr_k);  % A
    
        % --- Build DG object ---
        params = struct( ...
            'L',Lvals(k), 'C',Cvals(k), 'Rf',Rfvals(k), ...
            'RL', RL_k, 'Ibar', Ibar_k, ...           % ZI load parameters (→ plant)
            'Vrated', Vr_k, 'Irated', Ir_k, ...       % converter ratings
            'u_s', Vr_k, ...                         % start input near rated voltage
            'pos', pos(k,:), ...
            'x', [randn(1)*Vr_k; randn(1)*Ir_k], ...                       % initial state [vC;iL]
            'iload_fun', @(t) 0, ...                % no time ripple now (add later if needed)
            'K', [0 0] );
        DGs(k) = DG(k, params);
    end
    
    %% ----------- Proximity-based graph with min-length & length-proportional R -----------
    lines = {};
    idc   = 1;
    
    % --- Tunables ---
    minLen      = 1.0;                 % [arb units] minimum allowable line length
    proxRadius  = 3.0;                 % initial proximity radius (will grow if not connected)
    rho_per_len = 0.1;   %0.3             % [Ohm / length-unit] base resistance per unit length
    randSpan    = 0.2;     %0.2           % ±20% multiplicative randomness on R
    
    % Pairwise distances
    D = squareform(pdist(pos));        % NxN Euclidean distances
    
    % Helper to (re)build edge list for a given radius
    buildEdges = @(R) find(triu((D >= minLen) & (D <= R), 1));   % upper-triangular linear indices
    
    % Grow radius until connected (or until it’s obviously enough)
    maxTries = 10;
    % REVIEW PROPOSAL MAIN-04 (ADD): prevents an undefined comps variable if
    % all attempts produce zero edges.
    % comps = 1:N;
    for attempt = 1:maxTries
        Eidx  = buildEdges(proxRadius);
        [ii, jj] = ind2sub([N N], Eidx);
        Gtmp = graph(ii, jj, D(Eidx), N);
        if numedges(Gtmp) == 0
            proxRadius = proxRadius * 1.5;    % nothing yet; expand search
            continue;
        end
        comps = conncomp(Gtmp);
        if numel(unique(comps)) == 1
            break;                             % connected!
        else
            proxRadius = proxRadius * 1.35;    % not connected; expand and retry
        end
    end
    
    % If still not connected, force connectivity with an MST over *eligible* pairs
    if numel(unique(comps)) ~= 1
        % Allow all pairs >= minLen for MST
        W = D;
        W(D < minLen) = inf;                   % disallow too-short edges
        Gfull = graph(W);
        T = minspantree(Gfull,'Method','dense');
        ii = [ii; T.Edges.EndNodes(:,1)];
        jj = [jj; T.Edges.EndNodes(:,2)];
    end
    
    % REVIEW PROPOSAL MAIN-05 (ADD): remove duplicate edges introduced when
    % MST edges overlap the existing proximity graph.
    % edgePairs = unique(sort([ii,jj],2),'rows','stable');
    % ii = edgePairs(:,1); jj = edgePairs(:,2);
    % Create TransmissionLine objects; R ∝ length with randomness
    for e = 1:numel(ii)
        i = ii(e); j = jj(e);
        len_ij = D(i,j);
        if isinf(len_ij) || len_ij < minLen, continue; end
        mult = 1 + randSpan*(2*rand()-1);      % in [1-randSpan, 1+randSpan]
        R_ij = rho_per_len * len_ij * mult;
        lines{end+1} = TransmissionLine(idc, i, j, R_ij); 
        idc = idc + 1;
    end
    
    Lines = [lines{:}];
    
    % Assemble microgrid
    net = DCMicrogrid(DGs, Lines);
    net.buildSystemMatrices();
    
    % ------ Steady state with epsilon matrices (equal current sharing) ------
    net.solveSteadyState();
    net.K = zeros(N, 2*N);
    
    x_s = net.x_s;
    mag = 0.2;
    % x0 = x_s + x_s.*(mag*(2*rand(size(x_s))-1));
    x0 = (0.5+0.5*rand(size(x_s))).*x_s;
    % REVIEW PROPOSAL MAIN-06 (ADD): uncomment to override the active initial
    % condition with the milder symmetric perturbation already shown above.
    % x0 = x_s + x_s.*(mag*(2*rand(size(x_s))-1));
    net.setStateVector(x0);

end
