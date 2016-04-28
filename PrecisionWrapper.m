

% Top wrapper for aiding precision analysis
% Set precision switch to the following values
%   prec_switch = 0 - Default, all in DP
%   prec_switch = 1 - only smoother
%   prec_switch = 2 - only residual
%   prec_switch = 3 - only prolong_&_correct
%   prec_switch = 4 - smoother + residual
%   prec_switch = 5 - residual + prolong_&_correct
%   prec_switch = 6 - prolong_&_correct + smoother
%   prec_switch = 7 - smoother + residual + prolonbg_&_correct
%   prec_switch = 8 - Do All

%  set the function here
       % mu = @(x,y,z)(1e6 * sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z));
	   mu = @(x,y,z)(1 + 1e6 * (cos(2*pi*x)^2 + cos(2*pi*y)^2 + cos(2*pi*z)^2));
		
		prec_switch = 8;
		start_index = 1;
		end_index   = 1;
		if prec_switch==8
		     start_index=1; end_index=7;
		else
		     start_index=prec_switch ; end_index = prec_switch ;
		end
		
		disp('----- First Evaluating the default ---------')
		g = create_hexmesh_grids(3, mu, @homg.xform.shell, [1 2 4], [2 4 8]);
		disp('======= Solving using Multigrid as Solver ======');
		disp('----- Smoother : Jacobi with prec_switch = 0(Default) -----');
		[u_jac_0, rr_jac_0, iter_jac_0] = g.solve(150, 'jacobi', 3, 3, g.L, g.get_u0, 0) ;
		disp('----- End Jacobi ------');
		disp('----- Smoother : Chebyshev with prec switch = 0(Default) ----');
		[u_cheb_0, rr_cheb_0, iter_cheb_0] = g.solve(150, 'chebyshev', 3, 3, g.L, g.get_u0, 0);
		disp('----- End Chebyshev -------');
		disp('----- Smoother : ssor with prec switch = 0(Default) ----');
		[u_ssor_0, rr_ssor_0, iter_ssor_0] = g.solve(150, 'ssor', 2, 1, g.L, g.get_u0, 0);
		disp('----- End SSOR ------\n\n');
		
		if(start_index > 0)
		    for i=start_index:end_index
		        prec_switch = 0;
		    	disp(['Iteration-level - ' num2str(i) '-------']);
		    	disp('-----Creating Mesh Hierarchy-------');
		    	  %  g = create_hexmesh_grids(3, mu, @homg.xform.shell, [1 2 4], [2 4 8]);
		    	    disp('======= Solving using Multigrid as Solver =======');
		    		disp(['----- Smoother : Jacobi with prec_switch = ' num2str(i) '--------------']);
		    		if i==1		[u_jac_1, rr_jac_1, iter_jac_1] = g.solve(150, 'jacobi', 3, 3, g.L, g.get_u0 , i) ;
		    		elseif i==2		[u_jac_2, rr_jac_2, iter_jac_2] = g.solve(150, 'jacobi', 3, 3, g.L, g.get_u0 , i) ;
		    		elseif i==3		[u_jac_3, rr_jac_3, iter_jac_3] = g.solve(150, 'jacobi', 3, 3, g.L, g.get_u0 , i) ;
		    		elseif i==4		[u_jac_4, rr_jac_4, iter_jac_4] = g.solve(150, 'jacobi', 3, 3, g.L, g.get_u0 , i) ;
		    		elseif i==5		[u_jac_5, rr_jac_5, iter_jac_5] = g.solve(150, 'jacobi', 3, 3, g.L, g.get_u0 , i) ;
		    		elseif i==6		[u_jac_6, rr_jac_6, iter_jac_6] = g.solve(150, 'jacobi', 3, 3, g.L, g.get_u0 , i) ;
		    		elseif i==7		[u_jac_7, rr_jac_7, iter_jac_7] = g.solve(150, 'jacobi', 3, 3, g.L, g.get_u0 , i) ;
		    		end
		    		disp('------  End Jacobi -------------');
		    		
		    		disp(['------ Smoother : Chebyshev with prec_switch = ' num2str(i) '-------------']);
		    		if i==1		[u_cheb_1, rr_cheb_1, iter_cheb_1] = g.solve(150, 'chebyshev', 3, 3, g.L, g.get_u0 , i) ;
		    		elseif i==2		[u_cheb_2, rr_cheb_2, iter_cheb_2] = g.solve(150, 'chebyshev', 3, 3, g.L, g.get_u0 , i) ;
		    		elseif i==3		[u_cheb_3, rr_cheb_3, iter_cheb_3] = g.solve(150, 'chebyshev', 3, 3, g.L, g.get_u0 , i) ;
		    		elseif i==4		[u_cheb_4, rr_cheb_4, iter_cheb_4] = g.solve(150, 'chebyshev', 3, 3, g.L, g.get_u0 , i) ;
		    		elseif i==5		[u_cheb_5, rr_cheb_5, iter_cheb_5] = g.solve(150, 'chebyshev', 3, 3, g.L, g.get_u0 , i) ;
		    		elseif i==6		[u_cheb_6, rr_cheb_6, iter_cheb_6] = g.solve(150, 'chebyshev', 3, 3, g.L, g.get_u0 , i) ;
		    		elseif i==7		[u_cheb_7, rr_cheb_7, iter_cheb_7] = g.solve(150, 'chebyshev', 3, 3, g.L, g.get_u0 , i) ;
		    		end
		    		disp('------ End Chebyshev -----------');
		    		
		    		disp(['------ Smoother : SSOR with prec_switch = ' num2str(i) '--------------']);
		    		if i==1		[u_ssor_1, rr_ssor_1, iter_ssor_1] = g.solve(150, 'ssor', 2, 1, g.L, g.get_u0 , i) ;
		    		elseif i==2		[u_ssor_2, rr_ssor_2, iter_ssor_2] = g.solve(150, 'ssor', 2, 1, g.L, g.get_u0 , i) ;
		    		elseif i==3		[u_ssor_3, rr_ssor_3, iter_ssor_3] = g.solve(150, 'ssor', 2, 1, g.L, g.get_u0 , i) ;
		    		elseif i==4		[u_ssor_4, rr_ssor_4, iter_ssor_4] = g.solve(150, 'ssor', 2, 1, g.L, g.get_u0 , i) ;
		    		elseif i==5		[u_ssor_5, rr_ssor_5, iter_ssor_5] = g.solve(150, 'ssor', 2, 1, g.L, g.get_u0 , i) ;
		    		elseif i==6		[u_ssor_6, rr_ssor_6, iter_ssor_6] = g.solve(150, 'ssor', 2, 1, g.L, g.get_u0 , i) ;
		    		elseif i==7		[u_ssor_7, rr_ssor_7, iter_ssor_7] = g.solve(150, 'ssor', 2, 1, g.L, g.get_u0 , i) ;
		    		end
		    		
		    		
		    end
	  end
		    		
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
			   
			   
