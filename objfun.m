
function f = objfun(x, nb, T, R, Table)
% global Table;
% global nb;
% global T;
% global R;
%%
% load linedata.txt;              %% load line data
% branch=linedata;
% fb = linedata(:,1);
% tb = linedata(:,2);
% G = graph(fb,tb);
% tnb = length(fb);
% nb = 128;
% T=dfsearch(G,1,'edgetonew');
% 
% bKVA = 1000;                    % base kVA
% bKV = 4.16/sqrt(3);                      % base kV
% bZ = ((bKV)^2)*1000/bKVA;                % base Z
% resitance = ((branch(:,3))/bZ);          % resistance of the edge
% reactance = ((branch(:,4))/bZ);          % reactance of the edge
% 
% %% calculating R and X matrix
% R = zeros(nb);
% X  = zeros(nb);
% 
% for i = 1:(nb-1)
%     R(fb(i), tb(i)) = resitance(i);
%     R(tb(i) ,fb(i)) = R(fb(i), tb(i)) ;
%     X(fb(i), tb(i))= reactance(i);
%     X(tb(i) ,fb(i)) = X(fb(i), tb(i)) ;
%     
% end
% 
%  %% defining the unknowns 
%  Ap=1:nb-1;                                 % defining the variabes for P                   
%  Aq=nb:2*(nb-1);                            % defining the variabes for Q
%  Ai=2*(nb)-1:3*(nb-1);                      % defining the variabes for l(I^2)
%  Av=3*(nb)-1:4*(nb-1)+1;                    % defining the variabes for V
% Table = [T(:,1) T(:,2) Ap'  Aq'  Ai' Av'];  % creating Table for variabes P, Q ,l, V
% Da = 3*(nb)-2:4*(nb-1)+1;                   % voltage variables including substation
% Volttable = Da';

%% onjective function I^2*R

f=0;
% f = x(1);
for i =2:nb
    Parent = find(i == T(:,2));
    
    if isempty(Parent)  
       f = f;
    else
        f =f+x((Table(Parent,5)))*R(T(Parent,1),T(Parent,2));
    end
    
end

 
 
 
    
