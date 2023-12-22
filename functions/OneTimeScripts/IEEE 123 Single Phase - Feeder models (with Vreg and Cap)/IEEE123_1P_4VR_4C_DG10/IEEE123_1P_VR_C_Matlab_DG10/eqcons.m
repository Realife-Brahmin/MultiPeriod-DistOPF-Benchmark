
function[c, ceq] = eqcons(x)


load linedata.txt;              %% load line data

fb = linedata(:,1);
tb = linedata(:,2);
G = graph(fb,tb);
tnb = length(fb);
nb = 128;


T=dfsearch(G,1,'edgetonew');

%% defining variables

 Ap=1:nb-1;                          % defining the unknowns for phaseA
 Aq=nb:2*(nb-1);
 Ai=2*(nb)-1:3*(nb-1);
 Av=3*(nb)-1:4*(nb-1)+1;
Table = [T(:,1) T(:,2) Ap'  Aq'  Ai' Av'];


Da = 3*(nb)-2:4*(nb-1)+1;
Volttable = Da';


%% current equation lv = P^2 + Q^2
c = [];
 
for i =2:nb
    Parent = find(i == T(:,2));
    if isempty(Parent)  
       ceq = ceq;
    else
       Poc = find(T(Parent,1) == T(:,1));
        ceq(Parent) = x(Table(Parent,5))*x(Volttable(Poc(1))) -  (x(Table(Parent,3))^2 +  x(Table(Parent,4))^2);
      end
end

    