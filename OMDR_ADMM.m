
 function[S,A,V,U,PI,Lamb,M,N]=OMDR_ADMM(X,S0,V0,U0,PI0,Lamb0,M0,N0,R,varargin)

%On-line Minimum Dispersion Regularization-Alternating Direction Method of
%Multipliers (OMDC-ADMM) 
 
% See
 
% L. Nus, S. Miron, D. Brie, "An ADMM-based algorithm with minimum dispersion
%regularization for on-line blind unmixing of hyperspectral images". Preprint in 
%Chemometrics and Intelligent Laboratory Systems, 2020. 

%--------------------------------------------------------------
% Inputs
%--------------------------------------------------------------

% X: one slice of the hyperspectral image
% S0: initial endmembers
% A0: initial abundances
% V0, U0: initial auxiliary variables
% PI0, Lamb0: initial dual variables
% M0, N0: initial matrices (past information)
% R: number of endmembers (rank)
% rho:  convergence rate parameter
% alpha: weighting coefficient (tracking capability)
% mu: regularization parameter of the minimum dispersion regularization
% N1, N2: outer and inner loop iterations

%--------------------------------------------------------------
% Outputs
%--------------------------------------------------------------

% S, V, U, PI, Lamb, M, N

%--------------------------------------------------------------
% Initialization
%--------------------------------------------------------------

S = S0;
V = V0; U = U0;
PI = PI0; Lamb = Lamb0;

D = eye(R,R)-(1/R)*ones(1,R)'*ones(1,R); % Matrix relative to the minimum dispersion regularization

%--------------------------------------------------------------
% Set the default parameters
%--------------------------------------------------------------

rho = 0.01;
alpha = 0.99;
mu = 10^-3;
N1 = 100;
N2 = 10;

%--------------------------------------------------------------
% Read the input parameters
%--------------------------------------------------------------

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'RHO'
                rho = varargin{i+1};
            case 'ALPHA'
                alpha = varargin{i+1};
            case 'MU'
                mu = varargin{i+1};
            case 'N1'
                N1 = varargin{i+1};
            case 'N2'
                N2 = varargin{i+1};
                
            otherwise
                
                error(['Unrecognized option: ''' varargin{i} '''']);
        end;
    end;
end

%---------------------------------------------
%   Main body
%---------------------------------------------  

for t1 = 1:N1
    
  %%%% Update the abundances 

E = [(1 - alpha)*(S'*S) + rho*eye(R,R), ones(R,1);ones(1,R),0];
inverse_E = inv(E);

    for t2 = 1:N2
    
    AA = inverse_E*[((1 - alpha)*(S'*X) + rho*(V - PI));ones(1,size(V,2))];
    A = AA(1:R,:);
    V = max(0, A + PI);
    PI = PI + A - V;

    end

 %%%% Update the endmembers 

N = alpha*N0 + (1 - alpha)*X*A';
M = alpha*M0 + (1 - alpha)*A*A';
S = (N + rho*(U - Lamb))*inv(M + rho*eye(R,R) + 2*mu*D*eye(R,R));
U = max(0, S + Lamb);
Lamb = Lamb + S - U;

end


