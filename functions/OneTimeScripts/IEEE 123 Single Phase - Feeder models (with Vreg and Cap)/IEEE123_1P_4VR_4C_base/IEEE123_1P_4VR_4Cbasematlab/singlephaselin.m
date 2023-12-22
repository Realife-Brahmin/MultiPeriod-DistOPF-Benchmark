% clc;
% clear all;
function [V, x, Table,Volttable,tap]= singlephaselin()
load linedata.txt        % load the line data, it has R and X value for the edges
branch=linedata;
load powerdata.txt;    % load the power data, node number, P load, Q load, Capacitor, DG

mult = 1;                                % load mupltiplier
mult1 = 1;                               % PV mupltiplier

fb = branch(:,1);                        % from bus
tb = branch(:,2);                        % to bus
G = graph(fb,tb);                        % generate the tree of the network,
nb = size((powerdata),1);               % total number of nodes 
bKVA = 1000;                             % base kVA
bKV = 4.16/sqrt(3);                      % base kV
bZ = ((bKV)^2)*1000/bKVA;                % base Z
resitance = ((branch(:,3))/bZ);          % resistance of the edge
reactance = ((branch(:,4))/bZ);          % reactance of the edge

PL = powerdata(:,2)/bKVA;              % rated active power of node
QL = powerdata(:,3)/bKVA;              % rated reactive power of node
QC = powerdata(:,4)/bKVA;              % rated Q of capacitor 

PL = PL.*mult;                          % rated active power multipiled with load multiplier                     
QL = QL.*mult;                           % rated active power multipiled with load multiplier                     



T=dfsearch(G,1,'edgetonew');                % finding Parent child combination using dfsearch

%% calculating R and X matrix
R = zeros(nb);
X  = zeros(nb);
for i = 1:(nb-1)
    R(fb(i), tb(i)) = resitance(i);
    R(tb(i) ,fb(i)) = R(fb(i), tb(i)) ;
    X(fb(i), tb(i))= reactance(i);
    X(tb(i) ,fb(i)) = X(fb(i), tb(i)) ;    
end

%% Finding number of unknowns

 Ap = 1:nb-1;                          % defining the variabes for P
 Aq = nb:2*(nb-1);                      % defining the variabes for Q
 Av = 2*(nb):3*(nb-1)+1;                % defining the variabes for V

Table = [T(:,1) T(:,2) Ap'  Aq'   Av'];     % creating Table for variabes P, Q , V

Da = 2*(nb)-1:3*(nb-1)+1;               
Volttable = Da';                        % voltage variables including substation


for i = 1:(nb)                       % number of child of a node 
tc(i)=size(find((i)==T(:,1)),1) ;
row = find(i == T(:,1));
end
tc;

    
%% Formation of equality constarints
Aeq = zeros(381,382);
beq = zeros(381,1);

Vs= 1;                                % substation voltage

 for i =2:nb
    k = tc(i) ;                       %% total # of child
    row = find(i == T(:,1));
 
  if isempty(row)                       % if there is no child for a node
    Parent = find(i == T(:,2));         % find parent of a node
    
    if isempty(Parent)  
        Aeq = Aeq;

    else    
        Poc = find(T(Parent,1) == T(:,1));
        Aeq(Parent,Table(Parent,3))= 1;                                               % P12 = PL2
        Aeq(Parent+(nb-1),Table(Parent,4))= 1;                                        % Q12 = QL2-QC2
        Aeq(Parent+2*(nb-1),Table(Parent,5))= 1;                                      % V2 -V1 + 2rP+2xP =0                                           
        Aeq(Parent+2*(nb-1),Volttable(Poc(1)))= -1;
        Aeq(Parent+2*(nb-1),Table(Parent,3))= 2*(R(T((Parent),1),T((Parent),2)));
        Aeq(Parent+2*(nb-1),Table(Parent,4))= 2*(X(T((Parent),1),T((Parent),2)));
        beq(Parent)= 1*PL(Table(Parent,2));
        beq(Parent+(nb-1)) = 1*QL(Table(Parent,2))-QC(Table(Parent,2));
   
    end
 
  else                                                      % if the node has child nodes 
    Parent = find(i == T(:,2));
    
    if isempty(Parent)  
        Aeq = Aeq;
    
    else    
        Poc = find(T(Parent,1) == T(:,1));
        Aeq(Parent,Table(Parent,3))= 1;                             % P12- summation(P23) = PL2 
        Aeq(Parent,Table(row(1),3))= -1;
            for j = 1:length(row)-1
                Aeq(Parent, Table(row(j+1),3)) =   - 1;
            end
            
        Aeq(Parent+(nb-1),Table(Parent,4))= 1;                      % Q12- summation(Q23) = QL2-QC2 
        Aeq(Parent+(nb-1),(nb-1)+ Table(row(1),3))= -1;
            for j = 1:length(row)-1
                Aeq(Parent+(nb-1),(nb-1)+Table(row(j+1),3)) =  -  1;
            end
   
       Aeq(Parent+2*(nb-1),Table(Parent,5))= 1;                                 % V2 -V1 + 2rP+2xP =0  
       Aeq(Parent+2*(nb-1),Volttable(Poc(1)))= -1;
       Aeq(Parent+2*(nb-1),Table(Parent,3))= 2*(R(T((Parent),1),T((Parent),2)));
       Aeq(Parent+2*(nb-1),Table(Parent,4))= 2*(X(T((Parent),1),T((Parent),2)));
       beq(Parent)= 1*PL(Table(Parent,2));
       beq(Parent+(nb-1)) =  1*QL(Table(Parent,2))-QC(Table(Parent,2));
     end
  end
  
%  Aeq(3*(nb-1)+1,Volttable(1)) = 1;                                  % equation for substaion voltage
%  beq(3*(nb-1)+1) = Vs;
 
 end
Aeq;


%% defining tap
tap = [];
t_0 = 0.9;

for i = 1:32
    k = (1.1-0.9)/32;
    t_(i) = t_0 + i*k;
    tap = [tap t_(i)];
end
tap = [t_0, tap];

% VR_indx = [265,283,330];

Aeq(382,256) = 1 ;
beq(382) = (tap(18)^2)*Vs ;   

Aeq(265,265) = -(tap(12)^2);
Aeq(283,279) = -(tap(21)^2);
Aeq(330,330) = -(tap(13)^2);

% Aeq(382,256) = 1 ;
% beq(382) = (tap(17)^2)*Vs ;   
% 
% Aeq(265,265) = -(tap(17)^2);
% Aeq(283,279) = -(tap(17)^2);
% Aeq(330,330) = -(tap(17)^2);
nvar = size(Aeq,2);                                             % totale number of variables

%% formation of objective function

f = zeros(nvar,1);

%% limits defination
lb(1:254,1)= -(1500*ones(254,1));
lb(255:382,1)= ((0.8^2)*ones(128,1));
% lb(383:415,1)= zeros(33,1);

ub(1:254,1)= (1500*ones(254,1));
ub(255:382,1)= ((1.1)*ones(128,1));
% ub(383:415,1)= ones(33,1);

%%
optimoptions(@intlinprog,'Display','iter');
[x,fval,exitflag,output] = intlinprog(f,[],[],[],Aeq,beq,[],[]);

 V = sqrt(x(255:382));                  % voltage 
 Ptotal = x(1);                         % substaion active power
 Qtotal = x(128);                       % substaion reactive power

