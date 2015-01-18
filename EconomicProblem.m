classdef EconomicProblem < OCProblem
   
   properties
      ControlBounds
      
      % Problem Parameters
      delta
      p
      q
      c1
      c2
      r % allowed to be a constant OR function of time
      k
      
      % Let solver know that this is a maximization problem (flips the sign on
      % the objective during reporting and final solution)
      MinMax = 'Max';
   end
   
   methods
      function obj = EconomicProblem(params, ControlBounds)
         obj.ControlBounds = ControlBounds;
         obj.delta = params.delta;
         obj.p = params.p;
         obj.q = params.q;
         obj.c1 = params.c1;
         obj.c2 = params.c2;
         obj.r = params.r;
         obj.k = params.k;
      end
      
      function value = F(obj, t, x, E)
         N = x(1,:);
         
         c1 = obj.c1;
         c2 = obj.c2;
         delta = obj.delta;
         k = obj.k;
         p = obj.p;
         q = obj.q;
         
         if isa(obj.r, 'function_handle')
             r = obj.r(t);
         else
             r = obj.r;
         end
         
         value = [-E.*N.*q-N.*r.*(N./k-1.0);exp(-delta.*t).*(E.*c1+E.^2.*c2-E.*N.*p.*q)];
      end
      
      function value = dFdx_times_vec(obj, t, x, E, v)
         N = x(1,:);
         v1 = v(1,:);
         v2 = v(2,:);
         
         delta = obj.delta;
         k = obj.k;
         p = obj.p;
         q = obj.q;
         
         if isa(obj.r, 'function_handle')
             r = obj.r(t);
         else
             r = obj.r;
         end
         
         t2 = 1.0./k;
         value = [-v1.*(E.*q+r.*(N.*t2-1.0)+N.*r.*t2)-E.*p.*q.*v2.*exp(-delta.*t);0.0];
      end
      
      function value = dFdu_times_vec(obj, t, x, E, v)
         N = x(1,:);
         v1 = v(1,:);
         v2 = v(2,:);
         
         c1 = obj.c1;
         c2 = obj.c2;
         delta = obj.delta;
         p = obj.p;
         q = obj.q;
         
         value = -N.*q.*v1+v2.*exp(-delta.*t).*(c1+E.*c2.*2.0-N.*p.*q);
      end
   end
end

