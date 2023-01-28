function [V,S1, S_all, mic_iter, Dec_Var, time_dist, R, T] = NL_OPF_dist(Vs, Se,Area, end_area, loss_min, CVR)

V_min = 0.95;
V_max = 1.05;
%% Import Data

s1 = strcat('Area Data\Area',num2str(Area),'\linedata.txt');
load (s1);
branch=linedata;

s1 = strcat('Area Data\Area',num2str(Area),'\powerdata.txt');
load (s1);

powerdata = sortrows(powerdata,1);

load_mult = 1;
gen_mult = 1;

%% Bus Identifier:

GFI_bus = powerdata(find(powerdata(:,6)==3),1);

%% Graph Formation
fb = branch(:,1);
tb = branch(:,2);
G = graph(fb,tb);

global nb;
nb = size(powerdata,1);     % total number of nodes
%% Line Data

bKVA = 1000;                                        % base kVA
bKV = 4.16/sqrt(3);                                 % base kV
bZ = ((bKV)^2)*1000/bKVA;                           % base Z

resitance = ((branch(:,3))/bZ);                     % R of edge
reactance = ((branch(:,4))/bZ);                     % X of edge

PL = (powerdata(:,2).*load_mult)/bKVA;           % Rated Pload of node
QL = (powerdata(:,3).*load_mult)/bKVA;           % Rated Qload of node
QC = powerdata(:,4)/bKVA;                           % Rated Q of Capacitor

Pder = (powerdata(:,5).*gen_mult)/bKVA;          % Rated  DER active power
Sder = 1.2*powerdata(:,5)./bKVA; 

% % Large Pder in some buses (to match the result with GFI buses)
% Pder(GFI_bus) = Pder(GFI_bus)*10;                       % Rated  DER Complex power
% Sder(GFI_bus) = Sder(GFI_bus)*10;

if (~end_area)
    for j = 1:size(Se,1)
        PL(end-j+1) = real(Se(end-j+1));                 %in PU
        QL(end-j+1) = imag(Se(end-j+1));
        
    end
end


%% DER Configuration:
% find(Sder(:)~=0)--->   %DER connected BUS
DER_Bus = find(Sder(:)~=0);
S_DER = Sder(DER_Bus);   %in PU
P_DER = Pder(DER_Bus);   %in PU
l_bQ(:,1) = -sqrt((S_DER.^2)-(P_DER.^2));
u_bQ(:,1) = sqrt((S_DER.^2)-(P_DER.^2));

%%
% global T;
T = dfsearch(G,1,'edgetonew');

% global R;
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
Ap=1:nb-1;                                  % defining the variabes for P
Aq=nb:2*(nb-1);                             % defining the variabes for Q
Ai=2*(nb)-1:3*(nb-1);                       % defining the variabes for l(I^2)
Av=3*(nb)-1:4*(nb-1)+1;                     % defining the variabes for V
global Table;
Table = [T(:,1) T(:,2) Ap'  Aq'  Ai' Av'];  % creating Table for variabes P, Q ,l, V
Da = 3*(nb)-2:4*(nb-1)+1;                   % voltage variables including substation
Volttable = Da';

%% Initialization-

CVR_P = CVR (1);
CVR_Q = CVR (2);


%% A and b matrix formulation-

for i =2:nb
    
    % Which rows has the children buses for Node i:
    ChildB = find(i == T(:,1)) ;
    
    %Find the row where the unique Parent is located for Node i:
    ParentB = find(i == T(:,2)) ;
    
    %Find the row where the unique Parent of Node i, has all the child nodes:
    PoC = find(T(ParentB,1) == T(:,1));
    
    %% Aeq formulations
    %P equations
    Aeq(ParentB,Table(ParentB,3))= 1;
    Aeq(ParentB,Table(ParentB,5))= -(R(T((ParentB),1),T((ParentB),2)));
    Aeq(ParentB,Table(ParentB,6))= -(CVR_P/2)*PL(Table(ParentB,2));
    
    %Q equations
    Aeq(ParentB+(nb-1),Table(ParentB,4))= 1;
    Aeq(ParentB+(nb-1),Table(ParentB,5))= -(X(T((ParentB),1),T((ParentB),2)));
    Aeq(ParentB+(nb-1),(Table(ParentB,6)))= -(CVR_Q/2)*QL(Table(ParentB,2));
    
    % For nodes with child bus
    if ~isempty(ChildB)
        for j = 1:length(ChildB)
            Aeq(ParentB, Table(ChildB(j),3)) =   - 1;   % for P
            Aeq(ParentB+(nb-1),Table(ChildB(j),4)) =  -  1;   % for Q
        end
    end
    
    % V equations
    Aeq(ParentB+2*(nb-1),Table(ParentB,6))= 1;
    Aeq(ParentB+2*(nb-1),Volttable(PoC(1)))= -1;
    Aeq(ParentB+2*(nb-1),Table(ParentB,3))= 2*(R(T((ParentB),1),T((ParentB),2)));
    Aeq(ParentB+2*(nb-1),Table(ParentB,4))= 2*(X(T((ParentB),1),T((ParentB),2)));
    Aeq(ParentB+2*(nb-1),Table(ParentB,5))= -((R(T((ParentB),1),T((ParentB),2)))^2 + (X(T((ParentB),1),T((ParentB),2)))^2) ;
    
    %% beq Formulation
    beq(ParentB)=(1-(CVR_P/2))*PL(Table(ParentB,2))-Pder(Table(ParentB,2));
    beq(ParentB+(nb-1)) =  (1-(CVR_Q/2))*QL(Table(ParentB,2))-QC(Table(ParentB,2));
    
end

%% subtaion voltage equation
Aeq(3*(nb-1)+1,Volttable(1)) = 1;
beq(3*(nb-1)+1) = Vs;

%% DER equation addition
Table_DER = DER_Bus;
for k22 = 1:size(DER_Bus,1)
    Aeq((Table(find(DER_Bus(k22)==Table(:,2)),4)),end+1) = 1;
    
    %setting other parameters for DGs:
    Table_DER (k22,2) = size(Aeq,2);
    
    % slope kq definiton:
    Table_DER (k22,3) = 2*u_bQ(k22)/(V_max-V_min); % Qmax at Vmin, and vice versa 
%     Table_DER (k22,3) = 0.07;   % Known slope
    
    % Q_ref, V_ref definition:
    Table_DER (k22,4) = 0.00;  %Qref
    Table_DER (k22,5) = 1.00;  %Vref
end

%%

Tnvar = size(Aeq,2);         %% total number of variables

%% calling linear solution for intial point
[Vlin, xlin,Tablelin,Volttablelin]= singlephaselin(Vs,Se,Area,end_area);

Pflow0 = xlin(1:(nb-1));
Qflow0 = xlin(nb:(2*(nb-1)));
V0 =  Vlin;
Qder0 = xlin((end-(size(DER_Bus,1))+1):end);

for i=2:nb
    ParentB = find(i == T(:,2));
    PoC = find(T(ParentB,1) == T(:,1));
    Iflow0(Tablelin(ParentB,3))=(xlin(Tablelin(ParentB,3))^2+xlin(Tablelin(ParentB,4))^2)/xlin(Volttablelin(PoC(1)));
end

x0 = [Pflow0;Qflow0;Iflow0';V0;Qder0];


%% formation of objective function
% strt = ((nb-1)/2)+1;
% strt_end = strt+(nb-1);
% fun = @(x)( sum( x(strt:strt_end) ) );
% func = @(x)(0);

%% Definig Limits

lb(1,1) = -1500;                                        % this is to limit the power flow going reverse at the substation
lb(2:2*(nb-1),1)= (-1500*ones(2*(nb-1)-1,1));       % P Q limit
% lb(2:(nb-1),1)= (-150*zeros((nb-1)-1,1));       % P limit, stop flowing reverse
% lb(nb:2*(nb-1),1)= (-150*zeros((nb-1),1));       % Q limit, stop flowing reverse
lb(2*(nb-1)+1:3*(nb-1),1)= zeros((nb-1),1);         % I limit
lb(3*(nb-1)+1:4*(nb-1)+1,1)= ((V_min^2)*ones(nb,1)); % V limit


ub(1:2*(nb-1),1)= (1500*ones(2*(nb-1),1));          % P Q limit
ub(2*(nb-1)+1:3*(nb-1),1)= (1500*ones((nb-1),1));   % I limit
ub(3*(nb-1)+1:4*(nb-1)+1,1)= ((V_max^2)*ones(nb,1));      % V limit

lb = [lb ;l_bQ];
ub = [ub; u_bQ];

%%  Optimization -
% 
options = optimoptions('fmincon','Display','off','MaxFunctionEvaluations',100000000,'Algorithm','sqp');
tic
    % For Loss minimization:
    [x,fval,exitflag,output] = fmincon(@(x)objfun(x, R ,T),x0,[],[],Aeq,beq,lb,ub, @(x)eqcons(x, T), options);

time_dist = toc;
% mic_iter = output.iterations;


%% Result
P1 = x(1);
Q1 = x(nb);
Pall = x(1:nb-1);
Qall = x(nb:2*(nb-1));
Dec_Var_Q = (x(end-((size(DER_Bus,1))-1):end));

Dec_Var = zeros(nb,1);
for j=1:length(DER_Bus)
    Dec_Var(DER_Bus(j),1) = Dec_Var_Q(j);
end

S1 = complex(P1,Q1);    % In Pu
Sall = complex(Pall,Qall);

for j = 1:size(Table,1)
    V(Table(j,2),1)= x((end-nb+1-(size(DER_Bus,1)))+j);
    S_all((Table(j,2))-1,1) = Sall(j);
end
V(1) = Vs;
mic_iter = P1-(sum(PL)-sum(Pder));
fclose ('all');
% P1
% Q1



