classdef Test < matlab.unittest.TestCase
    
    % Test Method Block
    methods (Test)
        
        % Test Function
        function test_Qp_1d(testCase)      
            % Exercise function under test
            qp.H = 1;
            qp.f = -1;
            qp.lb = -1;
            qp.ub = 2;

            [x, fval, ~, ~, lam, t] = ip_qp(qp);
            
            testCase.assertEqual(x, 1, 'AbsTol', 1e-6);
            testCase.assertEmpty(lam.eqlin);
            testCase.assertEmpty(lam.ineqlin);
            testCase.assertEqual(lam.lower, 0, 'AbsTol', 1e-6);
            testCase.assertEqual(lam.upper, 0, 'AbsTol', 1e-6);
            testCase.assertEqual(t, [2; 1], 'AbsTol', 1e-6);
            testCase.assertEqual(fval, -1/2);
        end
        
        function test_Qp_1d_constraint_active(testCase)      
            % 
            qp.H = 1;
            qp.f = -1;
            qp.lb = 1.5;
            qp.ub = 2;
            qp.solver = 'quadprog';
            qp.options = optimoptions('quadprog');

            [x, fval, ~, ~, lam, t] = ip_qp(qp);
            [x1, fval1, ~, ~, lam1] = quadprog(qp);
            
            testCase.assertEqual(x, x1, 'AbsTol', 1e-6);
            testCase.assertEmpty(lam.eqlin);
            testCase.assertEmpty(lam.ineqlin);
            testCase.assertEqual(lam.lower, lam1.lower, 'AbsTol', 1e-6);
            testCase.assertEqual(lam.upper, lam1.upper, 'AbsTol', 1e-6);
            testCase.assertEqual(t, [0; 0.5], 'AbsTol', 1e-6);
            testCase.assertEqual(fval, fval1, 'AbsTol', 1e-6);
        end
        
        function test_Qp_2d_equality_constrained(testCase)      
            % 
            qp.H = [2, 0.3; 0.3, 4];
            qp.f = [-1; -2];
            qp.lb = [-1.5; -2.4];
            qp.ub = [3; 5];
            qp.Aeq = [2, -0.5];
            qp.beq = 3;
            qp.solver = 'quadprog';
            qp.options = optimoptions('quadprog');

            [x, fval, ~, ~, lam, t] = ip_qp(qp);
            [x1, fval1, ~, ~, lam1] = quadprog(qp);
            
            testCase.assertEqual(x, x1, 'AbsTol', 1e-6);
            testCase.assertEqual(lam.eqlin, lam1.eqlin);
            testCase.assertEmpty(lam.ineqlin);
            testCase.assertEqual(lam.lower, lam1.lower, 'AbsTol', 1e-6);
            testCase.assertEqual(lam.upper, lam1.upper, 'AbsTol', 1e-6);
%             testCase.assertEqual(t, [0; 0.5], 'AbsTol', 1e-6);
            testCase.assertEqual(fval, fval1, 'AbsTol', 1e-6);
        end
    end
end