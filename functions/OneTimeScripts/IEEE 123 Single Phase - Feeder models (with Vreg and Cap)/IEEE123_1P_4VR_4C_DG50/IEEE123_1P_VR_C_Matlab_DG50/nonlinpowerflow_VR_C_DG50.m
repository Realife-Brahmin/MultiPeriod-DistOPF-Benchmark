clc;
clear all;
tic
load linedata.txt;              % load the line data, it has R and X value for the edges
branch=linedata;
load powerdata.txt;            % load the power data, node number, P load, Q load, Capacitor, DG
fb = branch(:,1);               % from bus
tb = branch(:,2);               % to bus
G = graph(fb,tb);               % generate the tree of the network,

nb = size((powerdata),1);      % total number of nodes
bKVA = 1000;                    % base kVA
bKV = 4.16/sqrt(3);                      % base kV
bZ = ((bKV)^2)*1000/bKVA;                % base Z
resitance = ((branch(:,3))/bZ);          % resistance of the edge
reactance = ((branch(:,4))/bZ);          % reactance of the edge

loadmult =1;  
PVmult = 1;
PL = loadmult*powerdata(:,2)/bKVA;         % rated active power of node
QL = loadmult*powerdata(:,3)/bKVA;                  % rated reactive power of node
QC = powerdata(:,4)/bKVA;                  % rated Q of capacitor 
PG = PVmult*powerdata(:,5)/bKVA;              % rated active power of node    

T=dfsearch(G,1,'edgetonew');             % finding Parent child combination using dfsearch

%% calculating R and X matrix
R = zeros(nb);
X  = zeros(nb);


for i = 1:(nb-1)
    R(fb(i), tb(i)) = resitance(i);
    R(tb(i) ,fb(i)) = R(fb(i), tb(i)) ;
    X(fb(i), tb(i))= reactance(i);
    X(tb(i) ,fb(i)) = X(fb(i), tb(i)) ;
    
end

%% number of child for a node         
for i = 1:(nb)                       
tc(i)=size(find((i)==T(:,1)),1) ;
row = find(i == T(:,1));
end
tc;

 %% defining the unknowns 
 Ap=1:nb-1;                                 % defining the variabes for P                   
 Aq=nb:2*(nb-1);                            % defining the variabes for Q
 Ai=2*(nb)-1:3*(nb-1);                      % defining the variabes for l(I^2)
 Av=3*(nb)-1:4*(nb-1)+1;                    % defining the variabes for V
Table = [T(:,1) T(:,2) Ap'  Aq'  Ai' Av'];  % creating Table for variabes P, Q ,l, V
Da = 3*(nb)-2:4*(nb-1)+1;                   % voltage variables including substation
Volttable = Da';

 
%% calling linear solution for intial point
[Vlin, xlin,Tablelin,Volttablelin,tap]= singlephaselin();

%% formation of equality matrix 

Vs= 1;                                % substation voltage 

CVR_P = 0;                %% CVR factor for P = 0.6
CVR_Q = 0;                  %%% CVR factor for Q = 3

 for i =2:nb
    k = tc(i) ;
    row = find(i == T(:,1));
 
  if isempty(row)
    Parent = find(i == T(:,2));
   
    if isempty(Parent)  
        Aeq = Aeq;
    else
        
        Poc = find(T(Parent,1) == T(:,1));
        rowi = find(i == T(:,2));
  
    %%% P flow
    Aeq(Parent,Table(Parent,3))= 1;                                         %P12 -lr = PL2
    Aeq(Parent,Table(Parent,5))= -(R(T((Parent),1),T((Parent),2)));
    
    %%% Q flow
    Aeq(Parent+(nb-1),Table(Parent,4))= 1;                                  %Q12 -lx = QL2-QC2
    Aeq(Parent+(nb-1),Table(Parent,5))= -(X(T((Parent),1),T((Parent),2)));

   %%% V flow
   Aeq(Parent+2*(nb-1),Table(Parent,6))= 1;                                 % V2-V1+2rP+2xQ-(r^2+x^2)l = 0
   Aeq(Parent+2*(nb-1),Volttable(Poc(1)))= -1;
   Aeq(Parent+2*(nb-1),Table(Parent,3))= 2*(R(T((Parent),1),T((Parent),2)));
   Aeq(Parent+2*(nb-1),Table(Parent,4))= 2*(X(T((Parent),1),T((Parent),2)));
   Aeq(Parent+2*(nb-1),Table(Parent,5))= -((R(T((Parent),1),T((Parent),2)))^2 + (X(T((Parent),1),T((Parent),2)))^2) ;
   beq(Parent)= 1*PL(Table(Parent,2))-PG(Table(Parent,2));
   beq(Parent+(nb-1)) =  1*QL(Table(Parent,2))-QC(Table(Parent,2));
   
     end
 
 else  
    Parent = find(i == T(:,2));    
    
    if isempty(Parent)  
         Aeq = Aeq;
  
    else
     Poc = find(T(Parent,1) == T(:,1));
    
    %%% P flow
        Aeq(Parent,Table(Parent,3))= 1;
        Aeq(Parent,Table(row(1),3))= -1;
            for j = 1:length(row)-1
                Aeq(Parent, Table(row(j+1),3)) =   - 1;
            end     
        Aeq(Parent,Table(Parent,5))= -(R(T((Parent),1),T((Parent),2)));
    
      %%% Q flow 
          Aeq(Parent+(nb-1),Table(Parent,4))= 1;
          Aeq(Parent+(nb-1),(nb-1)+ Table(row(1),3))= -1;    
            for j = 1:length(row)-1
                Aeq(Parent+(nb-1),(nb-1)+Table(row(j+1),3)) =  -1;
            end
          Aeq(Parent+(nb-1),Table(Parent,5))= -(X(T((Parent),1),T((Parent),2)));
    
      %%% V 
           Aeq(Parent+2*(nb-1),Table(Parent,6))= 1;
           Aeq(Parent+2*(nb-1),Volttable(Poc(1)))= -1;
           Aeq(Parent+2*(nb-1),Table(Parent,3))= 2*(R(T((Parent),1),T((Parent),2)));
           Aeq(Parent+2*(nb-1),Table(Parent,4))= 2*(X(T((Parent),1),T((Parent),2)));
           Aeq(Parent+2*(nb-1),Table(Parent,5))= -((R(T((Parent),1),T((Parent),2)))^2 + (X(T((Parent),1),T((Parent),2)))^2) ;
  
   beq(Parent)=1*PL(Table(Parent,2))-PG(Table(Parent,2));
   beq(Parent+(nb-1)) = 1*QL(Table(Parent,2))-QC(Table(Parent,2));

end
    end
  Aeq(3*(nb-1)+1,Volttable(1)) = 1;
 beq(3*(nb-1)+1) = (tap(17)^2)*Vs;      % tap VR1 =0
%  beq(3*(nb-1)+1) = (tap(18)^2)*Vs;      % tap VR1 =1
 
 end
Aeq;

% Aeq(265,392) = -(tap(12)^2);      % tap VR2 =-5
% Aeq(283,406) = -(tap(21)^2);      % tap VR3 =4
% Aeq(330,457) = -(tap(13)^2);      % tap VR4 =-4

Aeq(265,392) = -(tap(17)^2);        % tap VR2 =0
Aeq(283,406) = -(tap(17)^2);        % tap VR3 =0
Aeq(330,457) = -(tap(17)^2);        % tap VR4 =0

%% objective function
fun = @(x)(0);          


Pflow0 = xlin(1:127);
Qflow0 = xlin(128:254);
V0 =  Vlin;


for i=2:nb
     Parent = find(i == T(:,2));
     Poc = find(T(Parent,1) == T(:,1));
     Iflow0(Tablelin(Parent,3))=(xlin(Tablelin(Parent,3))^2+xlin(Tablelin(Parent,4))^2)/xlin(Volttablelin(Poc(1)));
end
% Iflow0 = zeros(12,1);
x0 = [Pflow0;Qflow0;Iflow0';V0]; 

%% defining limits

lb(1:254,1)= (zeros(254,1));
lb(255:381,1)= (zeros(127,1));
lb(382:509,1)= ((0.8^2)*ones(128,1));

ub(1:254,1)= (1500*ones(254,1));
ub(255:381,1)= (1500*ones(127,1));
ub(382:509,1)= ((1.1)*ones(128,1));

options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',100000000,'Algorithm','sqp');
[x,fval,exitflag,output] = fmincon(fun,x0,[],[],Aeq,beq,[],[], @eqcons, options)


%% Results
Vnlin = sqrt(x(382:509));
Ptotal = x(1);
Qtotal = x(128);
% Pflowlin = xlin(1:12);
% Qflowlin = xlin(13:24);
 Pflownonlin = x(1:127)*1000;
Qflownonlin = x(128:254)*1000;
PflownonlinS = x(1)*1000
QflownonlinS = x(128)*1000

%% Nonlinear volatge using BFM

Vnonlin =  sqrt(x(Volttable(1)));        %% substaion voltage
for i = 2:nb
indx = find(i==Table(:,2));
Vnon1 =  sqrt(x(Table(indx,6)));
Vnonlin = [Vnonlin Vnon1];
end
Vnonlin = [Vs Vnonlin];

toc