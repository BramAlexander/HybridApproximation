function [Time,Y]= pair_based_model(T,g,I0,TSPAN)
%Solves the pair-based model for an arbitrary transmission network.
%T is a sparse matrix where T_{ij} is the strength of contact from site j
%to site i.
%g is the removal rate.
%I0 is a vector of initial infected sites.
%TSPAN is the list of time points at which to evaluate the solution.
%Returned is a matrix of length(TSPAN) rows where each column specifies the
%probability time course for the susceptible/infectious states of each
%individual and the states of all IS SS and II pairs on the
%network.

T=sparse(T);%Ensure that transmission matrix is sparse
T=T'; %Change transmission network so that T_{ij} is a contact from 
%site i to site j. Note that this also means that the output for pair
%quantities IS is also transposed.

%Set initial state of infectious and susceptible populations.
N=length(T(:,1));
I_vec=zeros(N,1);I_vec(I0)=1;
S_vec=ones(N,1);S_vec(I0)=0;

%Generate sparse diagonal matrices for multiplication
I_mat=spdiags(I_vec,0,N,N);
S_mat=spdiags(S_vec,0,N,N);

G=spones(T);%Contact network.

%Symmetrise contact network.
G_A=G.*(1-G');
H=G+G_A';

%Obtain information on the presence of closed and open triples.
Q=min(H^2,1);Q=Q-diag(diag(Q));%Nodes joined by a path-length of 2.
H_c=Q.*H;%Nodes joined by a path-length of 1 and 2.
H_o=Q-H_c;%Nodes joined by a path length of 2 but not 1.

%Set initial conditions for pairs.
F_AB=H;
F_AA=tril(H);

IS=I_mat*F_AB*S_mat;
SS=S_mat*F_AA*S_mat;%Symmetry means that we only need half the matrix.
II=I_mat*F_AA*I_mat;

%Determine the locations of each pair.
W_AB=find(reshape(F_AB,N^2,1));
W_AA=find(reshape(F_AA,N^2,1));
d_AB=length(W_AB);
d_AA=length(W_AA);

Y0=[S_vec;I_vec;IS(W_AB);SS(W_AA);II(W_AA)];

options=odeset('abstol',1e-4); %Improves numerical stability for some networks.

[Time,Y]=ode45(@model_function,TSPAN,Y0,options,T,H_c,H_o,g,N,W_AA,W_AB,d_AA,d_AB);

%______________________________________________________
function [dY] = model_function(t,Y0,T,H_c,H_o,g,N,W_AA,W_AB,d_AA,d_AB)

Y0=max(Y0,0);%Ensures that numerical errors do not cause negative values.

IS=spalloc(N,N,d_AB);
SS=spalloc(N,N,d_AB);
II=spalloc(N,N,d_AB);

%Generate vectors of individuals and matrices of pairs.
S=Y0(1:N);
I=Y0(N+1:2*N);
IS(W_AB)=Y0(2*N+1:2*N+d_AB);
SS(W_AA)=Y0(2*N+d_AB+1:2*N+d_AB+d_AA);
II(W_AA)=Y0(2*N+d_AB+d_AA+1:2*N+d_AB+2*d_AA);
SS=SS+SS'; %Regenerate symmetric matrices.
II=II+II';

%Diagonal matrices of reciprocal values.
inv_S=spdiags(spfun(@inve,S),0,N,N);
inv_I=spdiags(spfun(@inve,I),0,N,N);

%_______________________________________________________
%Kirkwood-type closure for triples.
%Here r denotes "right" and l denoted "left" indicating the
%relevant direction of network contact.

R=T.*IS;

%IrSrS triple.
Q=inv_I*(IS.*H_c);
IrSrS=(R'*Q).*(inv_S*SS*inv_S);%Closed part.
IrSrS=IrSrS+R'*H_o.*(inv_S*SS);%Addition of open part.

%IrSlI triple.
Q=(II.*H_c)*inv_I;
IrSlI=(inv_I*IS*inv_S).*(Q*R);%Closed part.
IrSlI=IrSlI+IS*inv_S.*(H_o*R);%Addition of open part.

SrSlI=IrSrS';%Due to symmetrising of networks.
IrSrI=IrSlI';

%_________________________________________________________
%Pair-based model equations
dT=sum(R)';
dS=-dT;
dI=dT-g*I;
dSS=-IrSrS-SrSlI;
dIS=IrSrS-IrSlI-R-g*IS;
dII=IrSlI+IrSrI+R+R'-2*g*II;

dY=[dS;dI;dIS(W_AB);dSS(W_AA);dII(W_AA)];

%_________________________________________________________
function M_out=inve(M_in)
%Determines the reciprocal value.
M_out=M_in.^(-1);
