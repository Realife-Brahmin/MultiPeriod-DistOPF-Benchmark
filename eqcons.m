%[c, ceq] - NonLinear
%for Relaxed use [-ceq,c]
function[c, ceq] = eqcons(x, T, nb)
% global nb;
% global T;

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
        PoC = find(T(Parent,1) == T(:,1));
        ceq(Parent) = x(Table(Parent,5))*x(Volttable(PoC(1))) -  (x(Table(Parent,3))^2 +  x(Table(Parent,4))^2);
    end
end
