function C = mldivide(A,B)
% MRDIVIDE C = A\B
%   C = A\B matrix B is left divided by matrix A which is factorized.
%
%   See also SYMBOLMF, MULTIFRONTAL.

%   Copyright 2016 Yingzhou Li, Stanford University

if isa(A,'Multifrontal')
    
    MF = A;
    C = B;
    
    if MF.symm == 1
        LeftDivSymmRecursionUp  (MF.symboltree,MF.Ltree);
        LeftDivSymmRecursionDiag(MF.symboltree,MF.Dtree);
        LeftDivSymmRecursionDown(MF.symboltree,MF.Ltree);
    elseif MF.symm == 2
        LeftDivSymmRecursionUp  (MF.symboltree,MF.Ltree);
        LeftDivSymmRecursionDown(MF.symboltree,MF.Utree);
    end
    
else
    error('Multifrontal as a numeriter has not been implemented yet');
end

%=====================================================================
    function LeftDivSymmRecursionUp(Stree,Ltree)
        
        if strcmpi(Stree.type,'node')
            LeftDivSymmRecursionUp(Stree.ltree,Ltree.ltree);
            LeftDivSymmRecursionUp(Stree.rtree,Ltree.rtree);
        end
        
        idx = Stree.idx;
        actidx = Stree.actidx;
        C(idx,:) = Ltree.Matinv*C(idx,:);
        C(actidx,:) = C(actidx,:) - Ltree.AMatinv*C(idx,:);
        
    end

    function LeftDivSymmRecursionDiag(Stree,Dtree)
        
        if strcmpi(Stree.type,'node')
            LeftDivSymmRecursionDiag(Stree.ltree,Dtree.ltree);
            LeftDivSymmRecursionDiag(Stree.rtree,Dtree.rtree);
        end
        
        idx = Stree.idx;
        C(idx,:) = Dtree.Matinv*C(idx,:);
        
    end

    function LeftDivSymmRecursionDown(Stree,Utree)
        
        idx = Stree.idx;
        actidx = Stree.actidx;
        C(idx,:) = C(idx,:) - Utree.AMatinv'*C(actidx,:);
        C(idx,:) = Utree.Matinv'*C(idx,:);
        
        if strcmpi(Stree.type,'node')
            LeftDivSymmRecursionDown(Stree.rtree,Utree.rtree);
            LeftDivSymmRecursionDown(Stree.ltree,Utree.ltree);
        end
        
    end

end