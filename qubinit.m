% =========================================================================
% Initialization for both the algorithms. Here we initialize the number of
% qubit and the G matrix
% =========================================================================

%tot_qub=input('       Please enter the number of qubits: \n');
% tot_qub=5;
% G=[3 -1 0 0 -1;...
%     -1 3 -1 0 0;...
%     0 -1 3 -1 0;...
%     0 0 -1 3 -1;...
%     -1 0 0 -1 30];
tot_qub=4;




 % List of all the G matrices that work
 G=[ 0   1  0  0;...
     1   0  0  0;...
     0   0  0  0;...
     0   0  0  0];
 G=[ 0   0  1  0;...
     0   0  0  0;...
     1   0  0  0;...
     0   0  0  0];
% It seems safe to say that all the 2 by 2 matrice subspaces will work

 G=[ 0   0  0  0;...
     0   0  1  1;...
     0   1  0  1;...
     0   1  1  0];
% some weird case but still works
 G=[ 0   1  0  0;...
     1   0  1  1;...
     0   1  0  1;...
     0   1  1  0];
G= [ 0   1  1  0;...
     1   0  1  0;...
     1   1  0  0;...
     0   0  0  0];
 
 G=[ 0   1  0  1;...
     1   0  1  0;...
     0   1  0  1;...
     1   0  1  0];
G=[  0   1  0  1;...
     1   0  1  1;...
     0   1  0  1;...
     1   1  1  0];
% 
 
G=[  0   1  1  1;...
     1   0  1  1;...
     1   1  0  1;...
     1   1  1  0];
 
 G=[  0   1  0  1;...
     1   0  1  0;...
     0   1  0  1;...
     1   0  1  0];
G=[ 0    1   1   1;...
    1    0   0   0;...
    1    0   0   0;...
    1    0   0   0];

G=1/2*[ 0    1    1   1;...
    1    0    1   1;...
    1    1    0   1;...
    1    1    1   0];

G=[  0   1  1  1;...
    1   0  1  0;...
    1   1  0  0;...
    1   0  0  0];

G=[ 2.5  1  0  1;...
    1  0.5  1  0;...
    0  1  0.3  1;...
    1  0  1  0.2];
G=ones(4);
G(1,1)=2; G(4,4)=2;
G=rand(4);
G=triu(G)+triu(G,1)';
 
 G=[ 0   1  0  1;...
     1   0  1  0;...
     0   1  0  1;...
     1   0  1  0];

% G=[  0   1  1  1;...
%      1   0  1  1;...
%      1   1  0  1;...
%      1   1  1  0];
% G=[ 2.5  1  1  1;...
%     1  0.5  1  1;...
%     1  1  0.3  1;...
%     1  1  1  0.2];


%FAILED --> PASSED
% G=[ 0  1  0  0;...
%     1  0  0  1;...
%     0  0  3  1;...
%     0  1  1  0];

% %FAILED --> PASSED
% G=[ 0  1  1  0;...
%     1  0  1  0;...
%     1  1  3  1;...
%     0  0  1  0];

%FAILED --> PASSED
% G=[ 5  1  2  1;...
%     1  2  1  3;...
%     2  1  3  1;...
%     1  3  1  4];

%FAILED --> PASSED
% G=[ 0  1  0  0;...
%     1  0  1  1;...
%     0  1  3  1;...
%     0  1  1  0];

%FAILED --> PASSED
% G=[ 0  1  0  0;...
%     1  0  1  1;...
%     0  1  1  1;...
%     0  1  1  0];

%FAILED --> PASSED
% G=[ 0  1  1  0;...
%     1  0  1  1;...
%     1  1  1  1;...
%     0  1  1  0];


%FAILED --> PASSED
% G=[ 0  1  1  0;...
%     1  0  1  1;...
%     1  1  0  1;...
%     0  1  1  1];


%FAILED --> PASSED
% G=[  0   1  0  0;...
%      1   0  1  1;...
%      0   1  0  1;...
%      0   1  1  0];
% G=G+diag(rand(4,1));

%FAILED --> PASSED
% G=[  1   1  1  1;...
%      1   -100  1  1;...
%      1   1  -100  1;...
%      1   1  1   1];
% G=G+diag(rand(4,1));

% %random couplings, no potential
% G=rand(4);
% G=G+G';
% G=G-diag(diag(G));
% G(2,2)=-10;
% G(3,3)=-10;
% G(2,3)=0;
% G(3,2)=0;

% G=[        0    0.9238    0    0.5293;...
%     0.9238  -10.0000         1      0;...
%     0         1  -10.0000      0.5049;...
%     0.5293    0    0.5049         0];


% tot_qub=7;
% V=[ 0.9593;...
%     0.5472;...
%     0.1386;...
%     0.1493;...
%     0.2575;...
%     0.8407;...
%     0.2543];
% lambda=0.1;
%
% T=[ 0 -1  0  0  0  0 -1;...
%    -1  0 -1  0  0  0  0;...
%     0 -1  0 -1  0  0  0;...
%     0  0 -1  0 -1  0  0;...
%     0  0  0 -1  0 -1  0;...
%     0  0  0  0 -1  0 -1;...
%     -1  0  0  0  0 -1 0 ];
% G=T+lambda*diag(V);

% tot_qub=6;
% G=[ 0 -1  0  0  0 -1;
%    -1  3 -1  0  0  0;
%     0 -1  3 -1  0  0;
%     0  0 -1  3 -1  0;
%     0  0  0 -1  3 -1;
%    -1  0  0  0 -1  3 ];



% % G=G+5*eye(5);
% % Operator=[];
% tot_qub=3;
% G=[0 1 1;1 0 1; 1 1 0];
% tot_qub=3;
%     G=[0.01 1 1;1 0.06 1;1 1 0.03];
