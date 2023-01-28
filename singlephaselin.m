function [V, x, Table,Volttable]= singlephaselin(Vs,Se,Area,end_area)

s1 = strcat('Area Data\Area',num2str(Area),'\linedata.txt');
load (s1);
branch=linedata;

s1 = strcat('Area Data\Area',num2str(Area),'\powerdata.txt');
load (s1);
powerdata = sortrows(powerdata,1);

load_mult = 1;
gen_mult = 1;

%% Graph Formation
fb = branch(:,1);
tb = branch(:,2);
G = graph(fb,tb);
tnb = length(fb);
nb = size(powerdata,1);

%% Line Data

bKVA = 1000;
bKV = 4.16/sqrt(3);
bZ = ((bKV)^2)*1000/bKVA;

resitance = ((branch(:,3))/bZ);
reactance = ((branch(:,4))/bZ);

PL = (powerdata(:,2).*load_mult)/bKVA;           %%% Pload of phase A
QL = (powerdata(:,3).*load_mult)/bKVA;           %%% Pload of phase B
QC = powerdata(:,4)/bKVA;

Pder = (powerdata(:,5).*gen_mult)/bKVA;
Sder = 1.2*powerdata(:,5)./bKVA;


if (~end_area)
    for j = 1:size(Se,1)
        PL(end-j+1) = real(Se(end-j+1));                 %in PU
        QL(end-j+1) = imag(Se(end-j+1));
        
    end
end


%% DER Configuration:
% find(Sder(:)~=0)--->   %DER connected BUS
S_DER = Sder(find(Sder(:)~=0));   %in PU
P_DER = Pder(find(Sder(:)~=0));   %in PU
l_bQ(:,1) = -sqrt((S_DER.^2)-(P_DER.^2));
u_bQ(:,1) = sqrt((S_DER.^2)-(P_DER.^2));

%%

T = dfsearch(G,1,'edgetonew');

R = zeros(nb);
X = zeros(nb);

% Matrix form of R and X in terms of graph
for i = 1:(nb-1)
    R(fb(i), tb(i)) = resitance(i);
    R(tb(i) ,fb(i)) = R(fb(i), tb(i)) ;
    X(fb(i), tb(i))= reactance(i);
    X(tb(i) ,fb(i)) = X(fb(i), tb(i)) ;
end

%% defining the unknowns for phaseA
Ap = 1:nb-1;
Aq = nb:2*(nb-1);
Av = 2*(nb):3*(nb-1)+1;
Table = [T(:,1) T(:,2) Ap'  Aq'   Av'];
Da = 2*(nb)-1:3*(nb-1)+1;
Volttable = Da';

%% Initialization-

CVR_P = 0;                %%% this will make the loads as constant power load
CVR_Q = 0;                %%% this will make the loads as constant power load
Aeq = zeros(3*(nb-1),Table(end));
beq = zeros(3*(nb-1),1);

%% A and b matrix formulation-

for i =2:nb
    i;
    % Which rows has the children buses for Node i:
    ChildB = find(i == T(:,1)) ;
    
    %Find the row where the unique Parent is located for Node i:
    ParentB = find(i == T(:,2)) ;
    
    %Find the row where the unique Parent of Node i, has all the child nodes:
    PoC = find(T(ParentB,1) == T(:,1));
    
    %% Aeq formulations
    %P equations
    Aeq(ParentB,Table(ParentB,3))= 1;
    Aeq(ParentB,Table(ParentB,5))= -(CVR_P/2)*PL(Table(ParentB,2));
    
    %Q equations
    Aeq(ParentB+(nb-1),Table(ParentB,4))= 1;
    Aeq(ParentB+(nb-1),(Table(ParentB,5)))= -(CVR_Q/2)*QL(Table(ParentB,2));
    
    % For nodes with child bus
    if ~isempty(ChildB)
        for j = 1:length(ChildB)
            Aeq(ParentB, Table(ChildB(j),3)) =   - 1;   % for P
            Aeq(ParentB+(nb-1),Table(ChildB(j),4)) =  -  1;   % for Q
        end
    end
    
    % V equations
    Aeq(ParentB+2*(nb-1),Table(ParentB,5))= 1;
    Aeq(ParentB+2*(nb-1),Volttable(PoC(1)))= -1;
    Aeq(ParentB+2*(nb-1),Table(ParentB,3))= 2*(R(T((ParentB),1),T((ParentB),2)));
    Aeq(ParentB+2*(nb-1),Table(ParentB,4))= 2*(X(T((ParentB),1),T((ParentB),2)));
    
    
    %% beq Formulation
    beq(ParentB)=(1-(CVR_P/2))*PL(Table(ParentB,2))-Pder(Table(ParentB,2));
    %     beq(ParentB+(nb-1)) =  (1-(CVR_Q/2))*QL(Table(ParentB,2));
    beq(ParentB+(nb-1)) =  (1-(CVR_Q/2))*QL(Table(ParentB,2))-QC(Table(ParentB,2));
    
end

%% subtaion voltage equation
Aeq(3*(nb-1)+1,Volttable(1)) = 1;
beq(3*(nb-1)+1) = Vs;

%% DER equation addition
DER_Bus = find(Sder(:)~=0);
for k22 = 1:size(DER_Bus,1)
    Aeq((Table(find(DER_Bus(k22)==Table(:,2)),4)),end+1) = 1;
end

%%

Tnvar = size(Aeq,2);         %% total number of variables

%% formation of objective function

f = zeros(Tnvar,1);
f(Table(1,3)) = 0;

for i=1:Tnvar             % total # of P,Q V
    ctype(i)='C';
end

lb(1) = 0;                  % this is to limit the power flow going reverse at the substation
lb(2:2*(nb-1),1)= (-1500*ones(2*(nb-1)-1,1));
lb(2*(nb-1)+1:3*(nb-1)+1,1)= ((0.95^2)*ones(nb,1));

ub(1:2*(nb-1),1)= (1500*ones(2*(nb-1),1));
ub(2*(nb-1)+1:3*(nb-1)+1,1)= (1.05^2*ones(nb,1));

lb = [lb ;l_bQ];
ub = [ub; u_bQ];

options = optimoptions('intlinprog','Display','off');
[x,fval,exitflag,output] = intlinprog(f,[],[],[],Aeq,beq,lb,ub);

for j = 1:size(Table,1)
    V2(Table(j,2),1)= x((end-nb+1-(size(DER_Bus,1)))+j);
end

V2(1) = 1.05^2;
V = sqrt(V2);
P1 = x(1);
Q1 = x(nb);
% end

