function [A,b,c,FEASIBLE]=myPreprocess(A,b,c,lbounds,ubounds,BIG)

% MYPREPROCESS This funcion is used to prcprocess the input data.
%
% [A,b,c,FEASIBLE]=myPreprocess(A,b,c,lbounds,ubounds,BIG)
%
% Check if the input date if feasible or not. 
% Make sure that A is full rank.
% Change the form of the input LP to standard form if necessary.
% Modifie from the preprocessing.m by Yin Zhang, 
% Department of Mathematics and Statistics,
% University of Maryland  Baltimore County, 1995
% 
% $ Verson 0.2 $         $ Date 21 Feb, 2012 $
% $ Yiming Yan $

%fprintf('Preprocessing ...\n');
[m,~] = size(A); FEASIBLE = 1;

if any(lbounds > ubounds)
   %fprintf('\nPreprocessor: Lower bound exceeds upper bound\n');
   %fprintf(1,'%c',7);  % ring a bell
   FEASIBLE = 0; return;
end;

if ~issparse(A) 
    A = sparse(A);
end;
b = sparse(b); c = sparse(c);
lbounds = sparse(lbounds);
ubounds = sparse(ubounds);

%--------------------------------
% delete fixed variables
%--------------------------------
fixed = lbounds == ubounds;
Fixed_exist = any(fixed);
if (Fixed_exist)
   ifix = find(fixed); 
   infx = find(1 - fixed);
   xfix = lbounds(ifix);
   c    = c(infx);
   b    = b - A(:,ifix)*sparse(xfix);
   A    = A(:,infx);
   lbounds = lbounds(infx);
   ubounds = ubounds(infx);
end 

%--------------------------------
% delete zero rows
%--------------------------------
rnnzct = sum(spones(A'));
if any(rnnzct == 0)
   izrows = rnnzct == 0;
   if any(b(izrows) ~= 0) 
      %fprintf('\nPreprocessor: problem infeasible\n');
      %fprintf(1,'%c',7);  % ring a bell
      FEASIBLE = 0; return;
   end;
   inzrows = find(rnnzct > 0);
   A = A(inzrows,:); b = b(inzrows); 
   rnnzct = rnnzct(inzrows);
   %fprintf(' (%i 0-rows)', length(izrows));
end

%--------------------------------
% make A structurally "full rank"
%--------------------------------
sprk = sprank(A');
if (sprk < m)
   [dmp, ~] = dmperm(A);
   irow = dmp(1:sprk);
   A = A(irow,:); b = b(irow); 
   rnnzct = rnnzct(irow);
   %fprintf(' (%i dep-rows)', m-sprk);
end

%--------------------------------
% delete zero columns
%--------------------------------
zrcol = (max(abs(A)) == 0)';
if any(zrcol == 1)
   izrcol = find(zrcol);
   if any(c(izrcol) < 0 & ubounds(izrcol) > BIG-1)
      %fprintf('\nPreprocessor: problem unbounded below\n');
      fprintf(1,'%c',7);  % ring a bell
      FEASIBLE = 0; return;
   end
   inzcol = find(1 - zrcol);
   A = A(:,inzcol);
   c = c(inzcol);
   lbounds = lbounds(inzcol);
   ubounds = ubounds(inzcol);
   %fprintf(' (%i 0-columns)', nnz(zrcol));
end

%--------------------------------
% solve singleton rows
%--------------------------------
singleton = (rnnzct == 1);
nsgrows = nnz(singleton);
if nsgrows >= max(1, .01*size(A,1))
   isgrows = find(singleton);
   iothers = find(1 - singleton);
   %fprintf(' (%i singletons)',nsgrows);

   Atmp = A(isgrows,:); Atmp1 = spones(Atmp); btmp = b(isgrows);
   if nsgrows == 1 
      isolved  = find(Atmp1); 
      insolved = find(Atmp1 == 0); 
      xsolved  = b(isgrows)/Atmp(isolved);
   else
      colnnzct = sum(Atmp1);
      isolved  = find(colnnzct);
      insolved = find(colnnzct == 0);
      [ii, jj] = find(Atmp); 
      Atmp = Atmp(ii,jj); btmp = btmp(ii);
      xsolved  = btmp./diag(Atmp);
      if any(colnnzct >  1)
         repeat = diff([0; jj]) == 0;
         for i = 1:length(xsolved) - 1
             if repeat(i+1) && xsolved(i+1) ~= xsolved(i)
                %fprintf('\nPreprocessor: problem infeasible\n');
                fprintf(1,'%c',7);  % ring a bell
                FEASIBLE = 0; return;
             end;
         end;
         ii = find(~repeat); jj = ii;
         Atmp = Atmp(ii,jj); btmp = btmp(ii);
         xsolved  = btmp./diag(Atmp);
      end;
   end;

   if any(xsolved < lbounds(isolved)) || ...
      any(xsolved > ubounds(isolved))
      %fprintf('\nPreprocessor: problem infeasible\n');
      fprintf(1,'%c',7);  % ring a bell
      FEASIBLE = 0; return;
   end;

   b = b(iothers) - A(iothers,isolved)*xsolved;
   A = A(iothers, insolved);
   c = c(insolved);
   lbounds = lbounds(insolved);
   ubounds = ubounds(insolved);
end;

%-----------------------------------------
% change it to standard form if necessary
%-----------------------------------------
[m,n] = size(A);
standardForm = 1;

Lbounds_non0 = any(lbounds ~= 0);
if Lbounds_non0
    %fprintf('lbounds exisit\n');
    standardForm = 0;
end

U = spones(ubounds<BIG-1);
if any(U)
    %fprintf('ubounds exisit\n');
    standardForm = 0;
    indexu = find(U>0);
else    
    indexu = [];
end

if ~standardForm
    %fprintf('Change it to standard form...\n');
    
    tmp_n=length(indexu);
    tmp_A=A;
    M = sparse(tmp_n, n);
    M(:, indexu) = speye(tmp_n);
    
    A=[tmp_A sparse(m,tmp_n); M speye(tmp_n)];
    b=[b-tmp_A*lbounds; ubounds(indexu)-lbounds(indexu)];
    c=[c;sparse(tmp_n,1)];
end

end