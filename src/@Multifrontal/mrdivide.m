function C = mrdivide(A,B)
% C=A/B

if isa(B,'Multifrontal')
    
    MF = B;
    C = A;
    
    if MF.symm == 1
        RightDivSymmRecursionUp  (MF.symboltree,MF.Ltree);
        RightDivSymmRecursionDiag(MF.symboltree,MF.Dtree);
        RightDivSymmRecursionDown(MF.symboltree,MF.Ltree);
    elseif MF.symm == 2
        RightDivSymmRecursionUp  (MF.symboltree,MF.Utree);
        RightDivSymmRecursionDown(MF.symboltree,MF.Ltree);
    end
    
else
    error('Multifrontal as a numeriter has not been implemented yet');
end

%=====================================================================
    function RightDivSymmRecursionUp(Stree,Utree)
        
        if strcmpi(Stree.type,'node')
            RightDivSymmRecursionUp(Stree.rtree,Utree.rtree);
            RightDivSymmRecursionUp(Stree.ltree,Utree.ltree);
        end
        
        idx = Stree.idx;
        actidx = Stree.actidx;
        C(:,idx) = C(:,idx)*Utree.Matinv';
        C(:,actidx) = C(:,actidx) - C(:,idx)*Utree.AMatinv';
        
    end

    function RightDivSymmRecursionDiag(Stree,Dtree)
        
        if strcmpi(Stree.type,'node')
            RightDivSymmRecursionDiag(Stree.ltree,Dtree.ltree);
            RightDivSymmRecursionDiag(Stree.rtree,Dtree.rtree);
        end
        
        idx = Stree.idx;
        C(:,idx) = C(:,idx)*Dtree.Matinv;
        
    end

    function RightDivSymmRecursionDown(Stree,Ltree)
        
        idx = Stree.idx;
        actidx = Stree.actidx;
        C(:,idx) = C(:,idx) - C(:,actidx)*Ltree.AMatinv;
        C(:,idx) = C(:,idx)*Ltree.Matinv;
        
        if strcmpi(Stree.type,'node')
            RightDivSymmRecursionDown(Stree.ltree,Ltree.ltree);
            RightDivSymmRecursionDown(Stree.rtree,Ltree.rtree);
        end
        
    end

end