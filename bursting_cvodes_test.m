clear;
%%
addpath FSP;

rng('default')          % use the same random seed to reproduce numerical results
%% Adjustable experimental constraints
n_total_cells = 100;
n_mle_test = 100; % number of data sets to generate
t_measure_min = 0.01; %minimum time to measure after sorting
t_measure_max = 600; %maximum time tto measure after sorting
t_long = 7*24*3600; %time to grow the cells from completely off state (we assume the distribution reaches very close to stationary at this point)
i_threshold_min = 30;
i_threshold_max = 80;
i_thresholds = i_threshold_min:i_threshold_max;
t_measurements = linspace(t_measure_min, t_measure_max, 1000);
%% Parameters for the full model
model_name = 'bursting';
ON = 1;
S = [1 0;-1 0;0 1;0 -1]'; % stoichiometry matrix
nmax = [1 10000];
states = fsp_get_states(nmax);

true_params = [0.05, 0.08, 100, 0.001]';
params_min = [1.0e-4, 1.0e-4, 1.0e-4, 1.0e-8]';
params_max = [1.0e3, 1.0e3, 1.0e+3, 1.0]';

params_dim = length(true_params);
ind_prop = @(gene, rna) [gene==0, gene==1, gene==1, rna]; % the parameter-independent part of the propensities
Acell = fsp_get_matrices(S, ind_prop, nmax); % parameter-independent terms of the FSP matrix

% gene state is not observable
unobserved_species = [1];
C_obs = sparse(kron(eye(nmax(2)+1), ones(1, nmax(1)+1)));
ichange = 1:3;

% propensity function at the true parameters
prop = @(x) true_params'.*ind_prop(x(1), x(2));
%% Compute the matrix and the partial derivatives wrt parameters
A = spalloc(size(Acell{1}, 1), size(Acell{1}, 2), size(Acell{1},1)*(size(S,2)+1));
for i = 1:params_dim
    A = A + true_params(i)*Acell{i};
end

matvec = @(t, v) A*v;
dmatvec  = cell(params_dim, 1);
for i = 1:params_dim
    dmatvec{i} = @(t, v) Acell{i}*v;
end

f_out = @(t, p) p;
p0 = zeros(size(A,1), 1);
dp0 = zeros(length(p0)*params_dim, 1);
p0(1) = 1;
tspan = linspace(0.1, 500, 100);
svecs0 = mat2cell(dp0, size(p0,1)*ones(4,1), 1);
stop_cond = struct('eps', 1.0e-4, 'mv', matvec);
jac = @(t,x) A;

%%
tic
[y_out, fsp_status] = FspCVodeMex(tspan, matvec, p0, f_out, @fsp_stop_eval, stop_cond);
toc
%%
tic
opt = odeset('Jacobian',jac);
[T,Y] = ode23s(matvec, tspan, p0, opt);
toc
%%
function [i_stop, stop_status] = fsp_stop_eval(t, p, stop_cond)
    i_stop = 0.0;
    stop_status = '';
end

function X = fsp_get_states( nmax )
% Given the max molecular counts in the FSP, generate all states in the
% hyper-cubic FSP.
%
% Input:
% ====
% nmax: vector of length N for max molecule populations.
%
% Output:
% =====
%
% X( 1:nst, 1:N ) : array of states, nst = number of states, N = number of
% species.

N = length( nmax );

X = cell(1,N);
Xvec = cell(1,N);
for i = 1:N
    Xvec{i} = (0:nmax(i));
end
[X{1:N}] = ndgrid(Xvec{1:N}); % Find all states in the N-d hyper-rectangle ...
for i = 1:N
    X{i} = reshape(X{i}, numel(X{i}), 1);
end
X = cell2mat(X);
end

function Aterms = fsp_get_matrices(S, ind_prop, nmax)
% Generate the parameter-independent terms to form the FSP matrix on a
% hyper-rectangle.
%
% Arguments:
% ---------
% S         : stoichiometry matrix of size N x M, where N is the number of species,
%             M the number of reactions.
%
% ind_prop  : function handle to compute the parameter-independent part of
%             the propensities.
%
% nmax      : vector of maximum numbers of molecules of all species.
%
%
%% Initialize variables
N = size(S,1);   % Number of species.
M = size(S,2);   % Number of reactions.

%% Compute Inf. Gen. Matrix. on the box defined by the first N constraints
Nst = prod(nmax+1);                     % Total number of states.

% Find the propensities of all reactions at all states in the
% hyper-rectangle ...
X = cell(1,N);
Xvec = cell(1,N);
for i = 1:N
    Xvec{i} = (0:nmax(i));
end
[X{1:N}] = ndgrid(Xvec{1:N}); % Find all states in the N-d hyper-rectangle ...
for i = 1:N
    X{i} = reshape(X{i}, numel(X{i}), 1);
end
props = ind_prop(X{1:N});
%keyboard
Aterms = cell(1,M);
for mu = 1:M
    % transform mu^th stoichiometry vector into 1D coordinates ...
    stoich_1D = S(1,mu);
    for i = 2:N
        stoich_1D = stoich_1D + S(i,mu)*prod(nmax(1:i-1)+1);
    end

    X_new = cell(1,N);
    for i = 1:N
        X_new{i} = X{i} + S(i,mu);
    end
    vaild_constraints = (cell2mat(X_new) >= zeros(Nst,N)) & (cell2mat(X_new) <= repmat(nmax, Nst, 1));
    props_keep = props(:,mu); props_keep(min(vaild_constraints,[],2)~=1) = 0;
    %keyboard
    % adding the matrix terms corresponding to the mu-th reaction
%     Aterms{mu} = spdiags(props_keep,-stoich_1D,Nst,Nst)-spdiags(props(:,mu),0,Nst,Nst);
    Aterms{mu} = spdiags(props_keep,-stoich_1D,Nst,Nst)-spdiags(props_keep,0,Nst,Nst);
end
end

