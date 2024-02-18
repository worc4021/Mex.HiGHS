classdef highstest < matlab.unittest.TestCase

    methods (Test)
    
        function LP_test(testCase)
            c = [1;1];
            A = [0,1;
                 1,2;
                 3,2];
            b = [ -inf,   7;...
                     5,  15;...
                     6, inf];
            bnds = [0,   4;
                    1, inf];
            
            [xDense,fValDense,infoDense] = highs_lp(c,A,b,bnds);

            [xSparse,fValSparse,infoSparse] = highs_lp(c,sparse(A),b,bnds);

            testCase.verifyEqual(xDense,xSparse);

            testCase.verifyEqual(fValDense,fValSparse);

            testCase.verifyEqual(infoDense.column_status, infoSparse.column_status);

            testCase.verifyEqual(infoDense.row_status,infoSparse.row_status);
        end
    end

end