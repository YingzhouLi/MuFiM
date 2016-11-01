classdef SymbolMF
    properties
        M
        N
        symm
        symboltree
    end
    methods
        function SMF = SymbolMF(A)
            
            if nargin < 1
                return;
            end
            
            if ~issparse(A)
                error('Only support sparse matrix');
            end
            
            [SMF.M,SMF.N] = size(A);
            
            if ishermitian(A)
                SMF.symm = 1;
            elseif issymmetric(spones(A))
                % symbolic symmetric, factorized as LDU
                SMF.symm = 2;
            else
                SMF.symm = 0;
                error('Not supported yet');
            end
            
            if SMF.symm > 0
                SMF.symboltree = BiTreePartition(A);
            end
        end
    end
end