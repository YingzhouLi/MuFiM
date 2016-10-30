function C = mldivide(A,B)
% C=A\B

if isa(A,'Multifrontal')
    
    MF = A;
    C = B;
    
    if MF.symm == 1
        LeftDivSymmRecursionUp  (MF.symboltree,MF.Ltree);
        LeftDivSymmRecursionDiag(MF.symboltree,MF.Dtree);
        LeftDivSymmRecursionDown(MF.symboltree,MF.Ltree);
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
        C(idx,:) = Ltree.Linv*C(idx,:);
        C(actidx,:) = C(actidx,:) - Ltree.ALDinv*C(idx,:);
        
    end

    function LeftDivSymmRecursionDiag(Stree,Dtree)
        
        if strcmpi(Stree.type,'node')
            LeftDivSymmRecursionDiag(Stree.ltree,Dtree.ltree);
            LeftDivSymmRecursionDiag(Stree.rtree,Dtree.rtree);
        end
        
        idx = Stree.idx;
        C(idx,:) = Dtree.Dinv*C(idx,:);
        
    end

    function LeftDivSymmRecursionDown(Stree,Ltree)
        
        idx = Stree.idx;
        actidx = Stree.actidx;
        C(idx,:) = C(idx,:) - Ltree.ALDinv'*C(actidx,:);
        C(idx,:) = Ltree.Linv'*C(idx,:);
        
        if strcmpi(Stree.type,'node')
            LeftDivSymmRecursionDown(Stree.rtree,Ltree.rtree);
            LeftDivSymmRecursionDown(Stree.ltree,Ltree.ltree);
        end
        
    end

end