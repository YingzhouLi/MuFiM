function C = mrdivide(A,B)
% C=A/B

if isa(B,'Multifrontal')
    
    MF = B;
    C = A;
    
    if MF.symm == 1
        RightDivSymmRecursionUp  (MF.symboltree,MF.Ltree);
        RightDivSymmRecursionDiag(MF.symboltree,MF.Dtree);
        RightDivSymmRecursionDown(MF.symboltree,MF.Ltree);
    end
    
else
    error('Multifrontal as a numeriter has not been implemented yet');
end

%=====================================================================
    function RightDivSymmRecursionUp(Stree,Ltree)
        
        if strcmpi(Stree.type,'node')
            RightDivSymmRecursionUp(Stree.rtree,Ltree.rtree);
            RightDivSymmRecursionUp(Stree.ltree,Ltree.ltree);
        end
        
        idx = Stree.idx;
        actidx = Stree.actidx;
        C(:,idx) = C(:,idx)*Ltree.Linv';
        C(:,actidx) = C(:,actidx) - C(:,idx)*Ltree.ALDinv';
        
    end

    function RightDivSymmRecursionDiag(Stree,Dtree)
        
        if strcmpi(Stree.type,'node')
            RightDivSymmRecursionDiag(Stree.ltree,Dtree.ltree);
            RightDivSymmRecursionDiag(Stree.rtree,Dtree.rtree);
        end
        
        idx = Stree.idx;
        C(:,idx) = C(:,idx)*Dtree.Dinv;
        
    end

    function RightDivSymmRecursionDown(Stree,Ltree)
        
        idx = Stree.idx;
        actidx = Stree.actidx;
        C(:,idx) = C(:,idx) - C(:,actidx)*Ltree.ALDinv;
        C(:,idx) = C(:,idx)*Ltree.Linv;
        
        if strcmpi(Stree.type,'node')
            RightDivSymmRecursionDown(Stree.ltree,Ltree.ltree);
            RightDivSymmRecursionDown(Stree.rtree,Ltree.rtree);
        end
        
    end

end