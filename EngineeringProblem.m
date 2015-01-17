classdef EngineeringProblem
   properties
      ControlBounds
      
      % Problem Parameters
      delta
      p
      q
      c1
      c3
      r
      k
      
      % Let solver know that this is a maximization problem (flips the sign on
      % the objective during reporting and final solution)
      MinMax = 'Max';
   end
   
   methods
      function obj = EngineeringProblem(params, ControlBounds)
         obj.ControlBounds = ControlBounds;
         obj.delta = params.delta;
         obj.p = params.p;
         obj.q = params.q;
         obj.c1 = params.c1;
         obj.c3 = params.c3;
         obj.r = params.r;
         obj.k = params.k;
      end
      
      function value = F(obj, t, x, u)
         N = x(1,:);
         E = x(2,:);
         
         c1 = obj.c1;
         c3 = obj.c3;
         delta = obj.delta;
         k = obj.k;
         p = obj.p;
         q = obj.q;
         r = obj.r;
         
         value  = [-E.*N.*q-N.*r.*(N./k-1.0);u;exp(-delta.*t).*(E.*c1+c3.*u.^2-E.*N.*p.*q)];
      end
      
      function value = dFdx_times_vec(obj, t, x, ~, v)
         N = x(1,:);
         E = x(2,:);
         v1 = v(1,:);
         v3 = v(3,:);
         
         c1 = obj.c1;
         delta = obj.delta;
         k = obj.k;
         p = obj.p;
         q = obj.q;
         r = obj.r;
         
         t2 = 1.0./k;
         t3 = exp(-delta.*t);
         value = [-v1.*(E.*q+r.*(N.*t2-1.0)+N.*r.*t2)-E.*p.*q.*t3.*v3;-N.*q.*v1+t3.*v3.*(c1-N.*p.*q);0.0];
      end
      
      function value = dFdu_times_vec(obj, t, ~, u, v)
         v2 = v(2,:);
         v3 = v(3,:);
         
         c3 = obj.c3;
         delta = obj.delta;
         
         value = v2+c3.*u.*v3.*exp(-delta.*t).*2.0;
      end
   end
end


