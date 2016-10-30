classdef Multifrontal < SymbolMF
    properties
        Ltree
        Utree
        Dtree
    end
    methods
        function MF = Multifrontal(A,SMF)
            
            if nargin < 2
                SMF = SymbolMF(A);
            end
            
            MF.M          = SMF.M;
            MF.N          = SMF.N;
            MF.symm       = SMF.symm;
            MF.symboltree = SMF.symboltree;
            
            MF = Factorization(MF,A);
            
        end
    end
end