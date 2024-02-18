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

        function linprog(testCase)
            if license('test','optimization_toolbox')
                A = [1 1
                    1 1/4
                    1 -1
                    -1/4 -1
                    -1 -1
                    -1 1];
                
                b = [2; 1; 2; 1; -1; 2];
                Aeq = [1 1/4];
                beq = 1/2;
                f = [-1; -1/3];
    
    
                xMatlab = linprog(f,A,b,Aeq,beq);
                Ahighs = [A;Aeq];
                bhighs = [-inf(size(b)),b;beq,beq];
                xHighs = highs_lp(f,Ahighs,bhighs);
                testCase.verifyEqual(xMatlab,xHighs);
            else
                testCase.verifyTrue(true);
            end
        end

        function milp(testCase)
            weightIngots = [5,3,4,6];
            costIngots = weightIngots.*[350,330,310,280];
            costAlloys = [500,450,400];
            costScrap = 100;
            carbonIngots = [5,4,5,3]/100;
            carbonAlloys = [8,7,6]/100;
            carbonScrap = 3/100;
            molybIngots = [3,3,4,4]/100;
            molybAlloys = [6,7,8]/100;
            molybScrap = 9/100;
            
            xBnds = [0,1;...    ingots
                     0,1;...    ingots
                     0,1;...    ingots
                     0,1;...    ingots
                     0,inf;...  alloy
                     0,inf;...  alloy
                     0,inf;...  alloy
                     0,inf];%   scrap

            f = [costIngots(:);costAlloys(:);costScrap];
            A = [weightIngots,ones(1,3),1;...
                 weightIngots.*carbonIngots,carbonAlloys,carbonScrap;...
                 weightIngots.*molybIngots,molybAlloys,molybScrap];
            b = [25,25;...
                 1.25,1.25;...
                 1.25,1.25];
            bInt = [true(4,1);false(3+1,1)];
            [~,fVal,info] = highs_lp(f,A,b,xBnds,bInt);

            testCase.verifyEqual(info.status,"Optimal");

            if license('test','optimization_toolbox')
                steelprob = optimproblem;
                ingots = optimvar('ingots',4,'Type','integer','LowerBound',0,'UpperBound',1);
                alloys = optimvar('alloys',3,'LowerBound',0);
                scrap = optimvar('scrap','LowerBound',0);
                cost = costIngots*ingots + costAlloys*alloys + costScrap*scrap;
                steelprob.Objective = cost;
                totalWeight = weightIngots*ingots + sum(alloys) + scrap;
                totalCarbon = (weightIngots.*carbonIngots)*ingots + carbonAlloys*alloys + carbonScrap*scrap;
                totalMolyb = (weightIngots.*molybIngots)*ingots + molybAlloys*alloys + molybScrap*scrap;
                steelprob.Constraints.conswt = totalWeight == 25;
                steelprob.Constraints.conscarb = totalCarbon == 1.25;
                steelprob.Constraints.consmolyb = totalMolyb == 1.25;
                opt = steelprob.optimoptions('Display','none');
                [~,fval] = solve(steelprob,'Options',opt);
                testCase.verifyEqual(fVal,fval);
            end
        end
    end

end