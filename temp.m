%% Block-diagonal stabilization obstructions: generator + tests + YALMIP/MOSEK attempt
% Continuous-time: xdot = A x + B u,  u = K x
% N subsystems, each n states, each local input dimension m=1 (changeable)

clear; clc; close all;

%% User parameters
n = 2;          % states per subsystem
N = 4;          % number of subsystems
m = 1;          % inputs per subsystem (kept at 1 for simplicity)
tol = 1e-8;     % numerical tolerance

% YALMIP/MOSEK settings
yalmip('clear');
ops = sdpsettings('solver','mosek','verbose',1,'cachesolvers',1);

%% Run 3 cases:
% Case 1: Obstruction 1 holds (unstable A-invariant subspace inside ker(B))
% Case 2: Obstruction 2 holds (decentralized fixed-mode via left eigenvector condition)
% Case 3: Both obstructions invalid (should be stabilizable with block-diag K)

for case_id = 1:3
    fprintf('\n==================== CASE %d ====================\n', case_id);
    [A,B,desc] = generate_case(case_id,n,N,m);
    fprintf('%s\n', desc);

    % Check obstructions
    [obs1,info1] = check_obstruction1(A,B,tol);
    [obs2,info2] = check_obstruction2(A,B,n,N,m,tol);

    fprintf('Obstruction-1 (unstable invariant subspace in ker(B)) : %d\n', obs1);
    if obs1, disp(info1); end
    fprintf('Obstruction-2 (block-diag fixed mode via left eigenvector) : %d\n', obs2);
    if obs2, disp(info2); end

    % Try to find block-diagonal K via convex sufficient LMI:
    % Find block-diag P>0 and block-diag Y such that
    %   (A P + B Y) + (A P + B Y)' < -eps I
    % then K = Y P^{-1} (block-diagonal).
    [feas,K,P,Y,diagn] = design_blockdiag_K_LMI(A,B,n,N,m,ops);

    fprintf('LMI sufficient design feasible? %d\n', feas);
    if feas
        eigAcl = eig(A + B*K);
        fprintf('Closed-loop eigenvalues (real parts):\n');
        disp(real(eigAcl).');
        fprintf('Max real part = %.4g\n', max(real(eigAcl)));
    else
        fprintf('Infeasible / failed. YALMIP code: %d\n', diagn.problem);
        fprintf('YALMIP msg: %s\n', diagn.info);
    end

    % Show open-loop
    eigA = eig(A);
    fprintf('Open-loop max real part = %.4g\n', max(real(eigA)));
end

%% ======================= FUNCTIONS =======================

function [A,B,desc] = generate_case(case_id,n,N,m)
% Generates (A,B) with B block-diagonal and desired obstruction behavior.
% Dimensions:
%   A in R^{Nn x Nn}
%   B in R^{Nn x Nm}, block-diag with Bi in R^{n x m}
%
% We keep m=1 for clarity, but code supports any m.

    nx = N*n;
    nu = N*m;

    % Helper: block-diag B
    Bi = cell(N,1);

    switch case_id
        case 1
            % Obstruction 1:
            % Put an unstable mode fully inside ker(B) and make it A-invariant.
            %
            % Construction:
            % - Use m=1 with Bi = [0;1;0;...], i.e., actuator does NOT touch state 1.
            % - Make state 1 of subsystem 1 unstable and decoupled so span(e1) is invariant.
            %
            % This guarantees: exists S = span(e_{(sub1,state1)}) subset ker(B),
            % AS subset S, and eigenvalue >0 on S.

            for k=1:N
                b = zeros(n,m);
                b(min(2,n),1) = 1;       % actuate state 2 (or closest), not state 1
                Bi{k} = b;
            end
            B = blkdiag(Bi{:});

            A = zeros(nx,nx);

            % Subsystem 1: make x1 unstable, and decoupled from everything else
            A(1,1) = +1.2;               % unstable scalar mode on x_{1,1}
            if n>=2
                A(2,2) = -0.5;
                A(2,1) = 0.0;
                A(1,2) = 0.0;
            end

            % Other subsystems stable, mild couplings among them (not touching state1 of sub1)
            for i=2:N
                idx = (i-1)*n + (1:n);
                Ai = -0.7*eye(n);
                % a little internal coupling
                if n>=2
                    Ai(1,2) = 0.2;
                    Ai(2,1) = -0.1;
                end
                A(idx,idx) = Ai;
            end

            % Add interconnections among subsystems 2..N (but NOT involving state1 of subsystem1)
            for i=2:N
                for j=2:N
                    if i~=j
                        idx_i = (i-1)*n + (1:n);
                        idx_j = (j-1)*n + (1:n);
                        A(idx_i,idx_j) = 0.05*randn(n);
                    end
                end
            end

            desc = 'Case 1: Has unstable A-invariant subspace inside ker(B) (impossible for any K).';

        case 2
            % Obstruction 2 (decentralized fixed mode via left eigenvector):
            %
            % Want an unstable eigenvalue lambda>0 with a left eigenvector q partitioned as
            % q_i^T B_i = 0 for all i, so lambda is fixed under any block-diag K.
            %
            % Construction:
            % - Set Bi = e1 (actuate state 1 only). Then ker(Bi^T)=span(e2,...,en).
            % - Choose A block-diagonal with each local block having eigenvalue lambda on state 2,
            %   so left eigenvector can be q_i = e2 (orthogonal to Bi).
            %
            % Take lambda = +0.8 on each subsystem's state 2; keep other states stable.

            if n < 2
                error('Case 2 needs n>=2 to create q_i orthogonal to B_i.');
            end

            for k=1:N
                b = zeros(n,m);
                b(1,1) = 1;              % actuate state 1 only
                Bi{k} = b;
            end
            B = blkdiag(Bi{:});

            lambda = 0.8;                % unstable fixed mode candidate
            A = zeros(nx,nx);

            for i=1:N
                idx = (i-1)*n + (1:n);
                Ai = -0.6*eye(n);
                Ai(2,2) = lambda;        % put unstable eigenvalue on state 2
                % small coupling from x1->x2 etc (still keeps e2 as left eigenvector if Ai is diagonal)
                % keep Ai diagonal so e2 is both left/right eigenvector
                A(idx,idx) = Ai;
            end

            % Add off-diagonal couplings that DO NOT destroy the left-eigenvector q with support on e2's:
            % Easiest: keep A strictly block-diagonal so proof is crystal clear.
            % (You can add structured couplings later if you want.)
            desc = 'Case 2: Has unstable decentralized fixed mode (left eigenvector orthogonal to each B_i).';

        case 3
            % No obstructions (construct a system that IS stabilizable with block-diag K).
            %
            % Construction:
            % - Each Bi actuates enough (choose Bi = [1;0;0...] so local controllability is plausible)
            % - Choose A with unstable diagonal blocks but locally stabilizable, and modest interconnections.
            % - Expect LMI to find a block-diag K.

            for k=1:N
                b = zeros(n,m);
                b(1,1) = 1;              % actuate state 1 (local authority)
                Bi{k} = b;
            end
            B = blkdiag(Bi{:});

            A = zeros(nx,nx);

            for i=1:N
                idx = (i-1)*n + (1:n);

                % Make each local block slightly unstable but stabilizable via B_i (acting on state 1)
                Ai = -0.2*eye(n);
                Ai(1,1) = +0.4;          % local instability in actuated channel
                if n>=2
                    Ai(2,1) = 0.8;       % couple actuated state into others
                    Ai(1,2) = 0.2;
                end
                A(idx,idx) = Ai;
            end

            % Add mild interconnections (won't trigger obstructions by construction)
            rng(1);
            for i=1:N
                for j=1:N
                    if i~=j
                        idx_i = (i-1)*n + (1:n);
                        idx_j = (j-1)*n + (1:n);
                        A(idx_i,idx_j) = 0.03*randn(n);   % small coupling
                    end
                end
            end

            desc = 'Case 3: Neither obstruction holds (should admit block-diag stabilization).';

        otherwise
            error('Unknown case_id.');
    end
end

function [obs1,info] = check_obstruction1(A,B,tol)
% Obstruction 1:
% Exists nonzero subspace S such that S subset ker(B), AS subset S,
% and A|S has eigenvalue with Re>=0.
%
% We numerically test S = ker(B) as a sufficient witness:
% - If ker(B) is A-invariant and A restricted to ker(B) has unstable eigenvalue => obstruction 1 holds.
%
% Note: obstruction 1 may still hold for a smaller S even if ker(B) isn't invariant.
% For our generated case 1, ker(B) contains an A-invariant unstable subspace so this check will pass.

    Z = null(full(B));      % basis for ker(B)
    info = struct();

    if isempty(Z)
        obs1 = false;
        info.msg = 'ker(B) is trivial.';
        return;
    end

    % Check approximate invariance: (I - P) A Z ~ 0
    P = Z * pinv(Z);
    resid = norm((eye(size(A,1)) - P) * (A*Z), 'fro');

    info.dimKerB = size(Z,2);
    info.invarianceResidual = resid;

    if resid > 1e-6
        obs1 = false;
        info.msg = 'ker(B) not numerically A-invariant (this does not rule out a smaller invariant S).';
        return;
    end

    % Restrict A to ker(B): A_S = Z^+ A Z
    AS = pinv(Z) * A * Z;
    eigAS = eig(AS);
    info.eigRestricted = eigAS;

    obs1 = any(real(eigAS) >= -tol);
    if obs1
        info.msg = 'Found unstable eigenvalue on A restricted to ker(B).';
    end
end

function [obs2,info] = check_obstruction2(A,B,n,N,m,tol)
% Obstruction 2:
% There exists unstable eigenvalue lambda (Re>=0) with left eigenvector q such that
% q_i' * B_i = 0 for all i.
%
% We compute eigenpairs of A' (right eigenvectors of A') which are left eigenvectors of A.

    info = struct();
    [V,L] = eig(A.');                % A' v = lambda v  => v' A = lambda v'
    lambdas = diag(L);

    unstable_idx = find(real(lambdas) >= -tol);
    info.unstableLambdas = lambdas(unstable_idx);

    obs2 = false;
    witness = [];

    for k = unstable_idx.'
        q = V(:,k);                  % left eigenvector of A (possibly complex)
        ok = true;

        for i=1:N
            rows = (i-1)*n + (1:n);
            cols = (i-1)*m + (1:m);
            Bi = B(rows,cols);

            if norm((q(rows).')*Bi, 2) > 1e-7
                ok = false;
                break;
            end
        end

        if ok
            obs2 = true;
            witness.lambda = lambdas(k);
            witness.q = q;
            break;
        end
    end

    if obs2
        info.msg = 'Found unstable lambda with left eigenvector blockwise orthogonal to each B_i.';
        info.witness_lambda = witness.lambda;
        % Store norms per block as evidence
        bn = zeros(N,1);
        for i=1:N
            rows = (i-1)*n + (1:n);
            cols = (i-1)*m + (1:m);
            Bi = B(rows,cols);
            bn(i) = norm((witness.q(rows).')*Bi, 2);
        end
        info.blockOrthogonalityNorms = bn;
    else
        info.msg = 'No unstable left eigenvector found that is orthogonal to every B_i block.';
    end
end

function [feas,K,Pval,Yval,diagn] = design_blockdiag_K_LMI(A,B,n,N,m,ops)
% Convex sufficient LMI for block-diagonal static state feedback:
% Variables:
%   P = blkdiag(Pi), Pi = Pi' > 0
%   Y = blkdiag(Yi)
% Constraint:
%   (A P + B Y) + (A P + B Y)' < -eps I
% Then K = Y P^{-1} (block diagonal).

    nx = size(A,1);
    eps_ = 1e-6;

    P = cell(N,1);
    Y = cell(N,1);

    for i=1:N
        P{i} = sdpvar(n,n,'symmetric');
        Y{i} = sdpvar(m,n,'full');
    end

    Pblk = blkdiag(P{:});
    Yblk = blkdiag(Y{:});

    M = A*Pblk + B*Yblk;
    Cons = [];
    Cons = [Cons, Pblk >= eps_*eye(nx)];
    Cons = [Cons, M + M' <= -eps_*eye(nx)];

    diagn = optimize(Cons, [], ops);

    feas = (diagn.problem == 0);
    if ~feas
        K = [];
        Pval = [];
        Yval = [];
        return;
    end

    Pval = value(Pblk);
    Yval = value(Yblk);

    % Recover block-diagonal K = Y * P^{-1}
    K = Yval / Pval;  % safe because Pval is block diag PD
end
