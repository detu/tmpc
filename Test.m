classdef Test < matlab.unittest.TestCase
    
    % Test Method Block
    methods (Test)
        
        % Test Function
        function testQp1(testCase)      
            % Exercise function under test
            qp.H = 1;
            qp.f = -1;
%             qp.Aeq = [];
%             qp.beq = [];
%             qp.A = [];
%             qp.b = [];
            qp.lb = -1;
            qp.ub = 2;

            [x, lam, mu, t] = ip_qp(qp);
            
            % Verify using test qualification
            % exp = your expected value
            testCase.assertEqual(x, 1, 'AbsTol', 1e-6);
            testCase.assertEmpty(lam);
            testCase.assertEqual(mu, [0; 0], 'AbsTol', 1e-6);
            testCase.assertEqual(t, [2; 1], 'AbsTol', 1e-6);
        end
    end
end