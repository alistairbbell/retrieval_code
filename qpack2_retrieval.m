% function L2 = qpack2_retrieval(M,R,QQ)
%
% This function does a Qpack2 retrieval, based on information
% about the measurement (M) and on the information about the
% retrieval (R). The QQ-Structure contains additional information
% for the Q-Structure. The information in the QQ-Structure will be
% replaced into the Q-Structure after the Q-Structure was created.
%
function L2 = qpack2_retrieval(M,R,QQ)

L2 = [];

%%% Load the measurements
disp(' ');
disp('Loading the measurements...')
msm = loadmsm(M);
if isempty(msm)
    disp('No measurements!!')
    return
elseif size(msm.tau,1)<20
    disp('Less than 20 measurements! Retrieval aborted!');
    return
end

%%% Create the Y-Structure required by Qpack2
disp('Creating the Y-structure...')
Y = msm2Y(msm);


%%% Create the Q-Structure required by Qpack2
disp('Creating the Q-structure...')
Q = createQ(Y,R);
if isempty(Q)
    disp('Could not create Q structure!')
    return
end
Y.Z_PLATFORM = Q.Z_SURFACE;
if exist('QQ') && ~isempty(QQ)
    QQnames = fieldnames(QQ);
    for i = 1:size(QQnames,1)
        eval(['Q.' QQnames{i} ' = QQ.' QQnames{i} ';']);
    end
end


%%% Check that all frequencies are OK
Y = rmfield(Y,'MIN_TIME');
Y = rmfield(Y,'MAX_TIME');
disp('Checking for consistency...')
if ~qp2_check_f( Q, Y, 1e3 );
    disp( 'Some mismatch between Q.F_BACKEND and frequencies of spectra.' );
    return
end


%%% Define OEM variables
O = qp2_l2( Q );  % This ensures that OEM returns the varibles needed
O.linear = true;


%%% Make inversion
disp('Doing inversion...')
tic
L2 = qpack2( Q, O, Y );
toc
if isempty(L2)
    disp('Could not do the inversion!')
    return
end

%%% Add additional information to the L2 structure:
L2.min_time = msm.min_time;
L2.max_time = msm.max_time;
L2.M = M;
L2.R = R;
L2.Q = Q;
L2.QQ = QQ;
L2.tint = msm.tint;



