classdef QpNode
    %QPNODE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Dependent)
        Q;
        R;
        S;
        q;
        r;
        C;
        D;
        ld;
        ud;
        lx;
        ux;
        lu;
        uu;
%         Zl;
%         Zu;
%         zl;
%         zu;
%         lls;
%         lus;
    end
    
    
    properties (Access=private)
        Q_;
        R_;
        S_;
        q_;
        r_;
        C_;
        D_;
        ld_;
        ud_;
        lx_;
        ux_;
        lu_;
        uu_;
%         Zl;
%         Zu;
%         zl;
%         zu;
%         lls;
%         lus;
    end
    
    
    methods
        function obj = QpNode(nx, nu, nc)
            if nargin < 1
                nx = 0;
            end
            
            if nargin < 2
                nu = 0;
            end
            
            if nargin < 3
                nc = 0;
            end
            
            obj.Q_ = zeros(nx, nx);
            obj.R_ = zeros(nu, nu);
            obj.S_ = zeros(nx, nu);
            obj.q_ = zeros(nx, 1);
            obj.r_ = zeros(nu, 1);
            obj.C_ = zeros(nc, nx);
            obj.D_ = zeros(nc, nu);
            obj.ld_ = zeros(nc, 1);
            obj.ud_ = zeros(nc, 1);
            obj.lx_ = zeros(nx, 1);
            obj.ux_ = zeros(nx, 1);
            obj.lu_ = zeros(nu, 1);
            obj.uu_ = zeros(nu, 1);
%             obj.Zl_ = zeros(ns, ns);
%             obj.Zu_ = zeros(ns, ns);
%             obj.zl_ = zeros(ns, 1);
%             obj.zu_ = zeros(ns, 1);
%             obj.lls_ = zeros(ns, 1);
%             obj.lus_ = zeros(ns, 1);
        end
       
        
        function obj = set.Q(obj, val)
            assert(isequal(size(val), size(obj.Q_)));
            assert(isequal(val, val.'));
            obj.Q_ = val;
        end
        
        
        function val = get.Q(obj)
            val = obj.Q_;
        end
        
        
        function obj = set.R(obj, val)
            assert(isequal(size(val), size(obj.R_)));
            assert(isequal(val, val.'));
            obj.R_ = val;
        end
        
        
        function val = get.R(obj)
            val = obj.R_;
        end
        
        
        function obj = set.S(obj, val)
            assert(isequal(size(val), size(obj.S_)));
            obj.S_ = val;
        end
        
        
        function val = get.S(obj)
            val = obj.S_;
        end
        
        
        function obj = set.q(obj, val)
            assert(isequal(size(val), size(obj.q_)));
            obj.q_ = val;
        end
        
        
        function val = get.q(obj)
            val = obj.q_;
        end
        
        
        function obj = set.r(obj, val)
            assert(isequal(size(val), size(obj.r_)));
            obj.r_ = val;
        end
        
        
        function val = get.r(obj)
            val = obj.r_;
        end
        
        
        function obj = set.C(obj, val)
            assert(isequal(size(val), size(obj.C_)));
            obj.C_ = val;
        end
        
        
        function val = get.C(obj)
            val = obj.C_;
        end
        
        
        function obj = set.D(obj, val)
            assert(isequal(size(val), size(obj.D_)));
            obj.D_ = val;
        end
        
        
        function val = get.D(obj)
            val = obj.D_;
        end
        
        
        function obj = set.ld(obj, val)
            assert(isequal(size(val), size(obj.ld_)));
            obj.ld_ = val;
        end
        
        
        function val = get.ld(obj)
            val = obj.ld_;
        end
        
        
        function obj = set.ud(obj, val)
            assert(isequal(size(val), size(obj.ud_)));
            obj.ud_ = val;
        end
        
        
        function val = get.ud(obj)
            val = obj.ud_;
        end
        
        
        function obj = set.lx(obj, val)
            assert(isequal(size(val), size(obj.lx_)));
            obj.lx_ = val;
        end
        
        
        function val = get.lx(obj)
            val = obj.lx_;
        end
        
        
        function obj = set.ux(obj, val)
            assert(isequal(size(val), size(obj.ux_)));
            obj.ux_ = val;
        end
        
        
        function val = get.ux(obj)
            val = obj.ux_;
        end
        
        
        function obj = set.lu(obj, val)
            assert(isequal(size(val), size(obj.lu_)));
            obj.lu_ = val;
        end
        
        
        function val = get.lu(obj)
            val = obj.lu_;
        end
        
        
        function obj = set.uu(obj, val)
            assert(isequal(size(val), size(obj.uu_)));
            obj.uu_ = val;
        end
        
        
        function val = get.uu(obj)
            val = obj.uu_;
        end
    end
end

