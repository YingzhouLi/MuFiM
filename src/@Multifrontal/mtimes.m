function C = mtimes(A,B)
% C=A*B

if isa(A,'Multifrontal')
    
    MF = A;
    C = B;
    
    if MF.symm == 1
        LeftMulSymmRecursionUp  (MF.symboltree,MF.Ltree);
        LeftMulSymmRecursionDiag(MF.symboltree,MF.Dtree);
        LeftMulSymmRecursionDown(MF.symboltree,MF.Ltree);
    end
    
elseif isa(B,'Multifrontal')
    
    MF = B;
    C = A;
    
    if MF.symm == 1
        RightMulSymmRecursionUp  (MF.symboltree,MF.Ltree);
        RightMulSymmRecursionDiag(MF.symboltree,MF.Dtree);
        RightMulSymmRecursionDown(MF.symboltree,MF.Ltree);
    end
    
end

%=====================================================================
    function LeftMulSymmRecursionUp(Stree,Ltree)
        
        if strcmpi(Stree.type,'node')
            LeftMulSymmRecursionUp(Stree.ltree,Ltree.ltree);
            LeftMulSymmRecursionUp(Stree.rtree,Ltree.rtree);
        end
        
        idx = Stree.idx;
        actidx = Stree.actidx;
        C(idx,:) = Ltree.L'*C(idx,:);
        C(idx,:) = C(idx,:) + Ltree.ALDinv'*C(actidx,:);
        
    end

    function LeftMulSymmRecursionDiag(Stree,Dtree)
        
        if strcmpi(Stree.type,'node')
            LeftMulSymmRecursionDiag(Stree.ltree,Dtree.ltree);
            LeftMulSymmRecursionDiag(Stree.rtree,Dtree.rtree);
        end
        
        idx = Stree.idx;
        C(idx,:) = Dtree.D*C(idx,:);
        
    end

    function LeftMulSymmRecursionDown(Stree,Ltree)
        
        idx = Stree.idx;
        actidx = Stree.actidx;
        C(actidx,:) = C(actidx,:) + Ltree.ALDinv*C(idx,:);
        C(idx,:) = Ltree.L*C(idx,:);
        
        if strcmpi(Stree.type,'node')
            LeftMulSymmRecursionDown(Stree.rtree,Ltree.rtree);
            LeftMulSymmRecursionDown(Stree.ltree,Ltree.ltree);
        end
        
    end

%=====================================================================
    function RightMulSymmRecursionUp(Stree,Ltree)
        
        if strcmpi(Stree.type,'node')
            RightMulSymmRecursionUp(Stree.rtree,Ltree.rtree);
            RightMulSymmRecursionUp(Stree.ltree,Ltree.ltree);
        end
        
        idx = Stree.idx;
        actidx = Stree.actidx;
        C(:,idx) = C(:,idx)*Ltree.L;
        C(:,idx) = C(:,idx) + C(:,actidx)*Ltree.ALDinv;
        
    end

    function RightMulSymmRecursionDiag(Stree,Dtree)
        
        if strcmpi(Stree.type,'node')
            RightMulSymmRecursionDiag(Stree.ltree,Dtree.ltree);
            RightMulSymmRecursionDiag(Stree.rtree,Dtree.rtree);
        end
        
        idx = Stree.idx;
        C(:,idx) = C(:,idx)*Dtree.D;
        
    end

    function RightMulSymmRecursionDown(Stree,Ltree)
        
        idx = Stree.idx;
        actidx = Stree.actidx;
        C(:,actidx) = C(:,actidx) + C(:,idx)*Ltree.ALDinv';
        C(:,idx) = C(:,idx)*Ltree.L';
        
        if strcmpi(Stree.type,'node')
            RightMulSymmRecursionDown(Stree.ltree,Ltree.ltree);
            RightMulSymmRecursionDown(Stree.rtree,Ltree.rtree);
        end
        
    end

end