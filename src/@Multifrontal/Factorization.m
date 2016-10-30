function MF = Factorization(MF,A)
%Factorization

if MF.symm == 1
    [MF.Ltree,MF.Dtree] = ...
        SymmFactorizationRecursion(MF.symboltree);
end


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
        
        Ltree.L = L;
        Ltree.Linv = Linv;
        Ltree.ALDinv = ALDinv;
        Dtree.D = D;
        Dtree.Dinv = Dinv;
        
    end

end