function MF = Factorization(MF,A)
%Factorization

if MF.symm == 1
    [MF.Ltree,MF.Dtree] = ...
        SymmFactorizationRecursion(MF.symboltree);
elseif MF.symm == 2
    [MF.Ltree,MF.Utree] = ...
        PatSymmFactorizationRecursion(MF.symboltree);
end

%====================================================================
    function [Ltree,Dtree] = SymmFactorizationRecursion(Stree)
        
        if strcmpi(Stree.type,'node')
            [Ltree.ltree,Dtree.ltree] = ...
                SymmFactorizationRecursion(Stree.ltree);
            
            [Ltree.rtree,Dtree.rtree] = ...
                SymmFactorizationRecursion(Stree.rtree);
        end
        
        idx = Stree.idx;
        actidx = Stree.actidx;
        [L,D] = ldl(full(A(idx,idx)));
        Linv = inv(L);
        Dinv = inv(D);
        
        ALDinv = A(actidx,idx)*Linv'*Dinv;
        
        A(actidx,actidx) = A(actidx,actidx) ...
            - sparse(A(actidx,idx)*Linv'*Dinv*Linv*A(actidx,idx)');
        
        Ltree.Mat = L;
        Ltree.Matinv = Linv;
        Ltree.AMatinv = ALDinv;
        Dtree.Mat = D;
        Dtree.Matinv = Dinv;
        
    end

%====================================================================
    function [Ltree,Utree] = PatSymmFactorizationRecursion(Stree)
        
        if strcmpi(Stree.type,'node')
            [Ltree.ltree,Utree.ltree] = ...
                PatSymmFactorizationRecursion(Stree.ltree);
            
            [Ltree.rtree,Utree.rtree] = ...
                PatSymmFactorizationRecursion(Stree.rtree);
        end
        
        idx = Stree.idx;
        actidx = Stree.actidx;
        [L,U] = lu(full(A(idx,idx)));
        Linv = inv(L);
        Uinv = inv(U);
        
        AUinv = A(actidx,idx)*Uinv;
        ALinv = (Linv*A(idx,actidx))';
        
        A(actidx,actidx) = A(actidx,actidx) ...
            - sparse(A(actidx,idx)*Uinv*Linv*A(idx,actidx));
        
        Ltree.Mat = L;
        Ltree.Matinv = Linv;
        Ltree.AMatinv = AUinv;
        Utree.Mat = U';
        Utree.Matinv = Uinv';
        Utree.AMatinv = ALinv;
        
    end

end