# This is an example code of steady incompressible Navier-Stokes equation
# based on the online tutorial of DOLFIN fenics, numerical example refers to
# the 2D lid-driven cavity given by the PhD thesis of professor
# J.A. Evans and other online resources. The code here based on the tutorial code
# provided by Prof. Kamensky's tIGAr package.

# Solving this problem requires a compatible Bspline discretization(DIV-free) with
# iterative penalty method and generalized-alpha method.

# because of the singularity in that domain, we need good initial guess for the iterative penalty method
# (Newton's method embedded), here we try RE = 10 as the starting point, up to 10000, max iter = 50, relative error norm 1e-5


from tIGAr import *
from tIGAr.NURBS import *
from igakit.nurbs import NURBS as NURBS_ik
from igakit.io import PetIGA
from tIGAr.compatibleSplines import *
from tIGAr.timeIntegration import *
import numpy as np
import math

parameters['form_compiler']['cpp_optimize'] = True  # may acclerate the assembly process
parameters["ghost_mode"] = "shared_facet"

# strain-rate tensor
def nablas(u):
	return 0.5*(grad(u) + grad(u).T)


# problem set up
for p in [1,2,3]: # k'
	QUAD_DEG =  p + 2 # highest deg in RT: p+1
	mesh = [16,32,128]
	form_compiler_parameters = {"quadrature_degree":QUAD_DEG}  # minimum quadrature_degree to ensure exact Integration
	for resol in mesh:
		# start to define B-spline mesh
		kv = [ uniformKnots(p,0.0,1.0,resol,False) for i in range(0,2) ]  # uniform knots in each direction

		controlMesh = ExplicitBSplineControlMesh([p,p],kv)                # B-spline control mesh

		fieldList = generateFieldsCompat(controlMesh,"RT",[p,p])          # soultion field basis over the control mesh, rt type
		splineGenerator = FieldListSpline(controlMesh,fieldList)

		# Strong normal BC
		for field in [0,1]: # u,v
			scalarSpline = splineGenerator.getFieldSpline(field)
			for side in [0,1]: # 1d: points, 2d: edges, 3d: faces
				sideDofs = scalarSpline.getSideDofs(field,side)
				splineGenerator.addZeroDofs(field,sideDofs)

		spline = ExtractedBSplineRT(splineGenerator,QUAD_DEG)
		n      = FacetNormal(spline.mesh)
		h      = CellDiameter(spline.mesh)
		h_avg  = 0.5*(h('+') + h('-'))
		dS  =  dS(metadata = {"quadrature_degree":QUAD_DEG}) #interior facets
		dx  =  dx(metadata = {"quadrature_degree":QUAD_DEG}) #interior
		x = spline.spatialCoordinates()

		# Nitsche's for tangential BC, set-up
		class Right(SubDomain):
			def inside(self,x,on_boundary):
				return near(x[0],1.0)

		class Left(SubDomain):
			def inside(self,x,on_boundary):
				return near(x[0],0.0)

		class Top(SubDomain):
			def inside(self,x,on_boundary):
				return near(x[1],1.0)

		class Bottom(SubDomain):
			def inside(self,x,on_boundary):
				return near(x[1],0.0)

		# Declare sub-subdomains
		top    = Top()
		bottom = Bottom()
		right  = Right()
		left   = Left()

		# Initialize mesh function for boundary domains
		boundaries = MeshFunction("size_t", spline.mesh, 1, spline.mesh.domains()) # 1 for edge, 2 for face
		boundaries.set_all(0)
		top.mark(boundaries, 1)
		bottom.mark(boundaries, 2)
		right.mark(boundaries, 3)
		left.mark(boundaries, 4)

		# Nitsche's parameters
		Cb  = 5.0*(p+1)                                      # penalty for Nitsche :  by J.Evans's dissertation(2011)
		ds = ds(subdomain_data = boundaries)
		ds  =  ds(metadata = {"quadrature_degree":QUAD_DEG}) # boundary facets
		uD  = Constant((1.0,0.0))                            # Dirichelt BC

		spline.relativeTolerance = 1e-4 # use relativly larger tol due to strong nonlinearity
		# Start method of continuation
		gamma = 1*10**(-p-1)
		for case in [2]: # case1: no model, case2: stabilized
			if case == 2:
				print('case2, stabilized')
				for Re in [10,100,250,400,600,800,1000,1200,1500,2000,2600,3200,3800,4400,5000,5600,6200,6800,7500,8200,8500,9000,10000]:
					nu = 1.0/Re
					if Re == 10:
						print("Starting Re = " + str(Re) + '_k = '  + str(p) + '_m = ' + str(resol))
						# always redefine, in case of unwanted initialization	        
						u_hat = Function(spline.V)
						v_hat = TestFunction(spline.V)

						u = spline.pushforward(u_hat)
						v = spline.pushforward(v_hat)

						# weak forms
						I_adv  = inner(grad(u)*u,v)*dx
						I_diff = 2.0*nu*inner(nablas(u),nablas(v))*dx

						# top edge, enforce u[0] = 1
						ut = as_vector([u[0],0.0])
						vt = as_vector([v[0],0.0])

						It_c = -2.0*nu*inner(nablas(ut)*n,vt) * ds(1)
						It_s = -2.0*nu*inner(nablas(vt)*n,ut-uD) * ds(1)
						It_p = Cb*nu/h*inner(vt,ut-uD) * ds(1)

						# bottom edge, u[0] = 0
						ub = as_vector([u[0],0.0])
						vb = as_vector([v[0],0.0])

						Ib_c = -2.0*nu*inner(nablas(ub)*n,vb) * ds(2)
						Ib_s = -2.0*nu*inner(nablas(vb)*n,ub) * ds(2)
						Ib_p = Cb*nu/h*inner(vb,ub) * ds(2)

						# left/right, u[1] = 0
						ulr = as_vector([0.0,u[1]])
						vlr = as_vector([0.0,v[1]])

						Ilr_c = -2.0*nu*inner(nablas(ulr)*n,vlr) * ( ds(3) + ds(4) )
						Ilr_s = -2.0*nu*inner(nablas(vlr)*n,ulr) * ( ds(3) + ds(4) )
						Ilr_p = Cb*nu/h*inner(vlr,ulr) * ( ds(3) + ds(4) )


						I_c = It_c + Ib_c + Ilr_c # consistency
						I_s = It_s + Ib_s + Ilr_s # symmetry
						I_p = It_p + Ib_p + Ilr_p # penalty

						I_Nit = I_c + I_s + I_p

						# solve it
						R     = I_adv + I_diff + I_Nit # residual, not stabilization
						spline.iteratedDivFreeSolve(R,u_hat,v_hat,penalty=Constant(1e6))
					else:
						print("Starting Re = " + str(Re) + '_k = '  + str(p) + '_m = ' + str(resol) + '_gamma = ' + str(gamma))
						# div-project the last result as new ic
						u_IG  = as_vector([u[0],u[1]])
						u_hat = spline.divFreeProject(u_IG)
						v_hat = TestFunction(spline.V)

						u,v   = spline.pushforward(u_hat), spline.pushforward(v_hat)


						I_adv  = inner(grad(u)*u,v)*dx
						I_diff = 2.0*nu*inner(nablas(u),nablas(v))*dx

						# top edge, enforce u[0] = 1
						ut = as_vector([u[0],0.0])
						vt = as_vector([v[0],0.0])

						It_c = -2.0*nu*inner(nablas(ut)*n,vt) * ds(1)
						It_s = -2.0*nu*inner(nablas(vt)*n,ut-uD) * ds(1)
						It_p = Cb*nu/h*inner(vt,ut-uD) * ds(1)

						# bottom edge, u[0] = 0
						ub = as_vector([u[0],0.0])
						vb = as_vector([v[0],0.0])

						Ib_c = -2.0*nu*inner(nablas(ub)*n,vb) * ds(2)
						Ib_s = -2.0*nu*inner(nablas(vb)*n,ub) * ds(2)
						Ib_p = Cb*nu/h*inner(vb,ub) * ds(2)

						# left/right, u[1] = 0
						ulr = as_vector([0.0,u[1]])
						vlr = as_vector([0.0,v[1]])

						Ilr_c = -2.0*nu*inner(nablas(ulr)*n,vlr) * ( ds(3) + ds(4) )
						Ilr_s = -2.0*nu*inner(nablas(vlr)*n,ulr) * ( ds(3) + ds(4) )
						Ilr_p = Cb*nu/h*inner(vlr,ulr) * ( ds(3) + ds(4) )


						I_c = It_c + Ib_c + Ilr_c
						I_s = It_s + Ib_s + Ilr_s
						I_p = It_p + Ib_p + Ilr_p

						I_Nit = I_c + I_s + I_p

						if Re >= 400:  # start to add stabilization term after Re>400
							# Define skeleton-stabilization term
							def Min(a, b): return (a+b-abs(a-b))/Constant(2.0)
						

							alpha = p - 1    # regularity
							Re_e  = avg( (dot(u,u) )**0.5)*h_avg/nu  # local Reynolds number
							#Note: this time we dont need to add the very small number
							eta_n = gamma * Min( Constant(1.0), Re_e ) * h_avg**( 2*alpha + 2) *  abs( dot(u('+'),n('+')) )

							# be careful to use the built-in jump operator, especially with higher order jumps. This is due to the convention: n('+') = - n('-')
							if p == 1:

								jump_u_n  = dot( ( grad(u)('+') - grad(u)('-') ) , n('+'))
								jump_v_n  = dot( ( grad(v)('+') - grad(v)('-') ) , n('+'))


							if p == 2:

								jump_u_n  = dot( grad( dot(grad(u)('+'),n('+')) )  - grad( dot(grad(u)('-'),n('+')) ) , n('+') )
								jump_v_n  = dot( grad( dot(grad(v)('+'),n('+')) )  - grad( dot(grad(v)('-'),n('+')) ) , n('+') )



							if p == 3:
								jump_u_n  = dot ( grad( dot( grad( dot(grad(u)('+'),n('+')) ), n('+') ) ) - grad( dot( grad( dot(grad(u)('-'),n('+')) ), n('+') ) ) , n('+'))
								jump_v_n  = dot ( grad( dot( grad( dot(grad(v)('+'),n('+')) ), n('+') ) ) - grad( dot( grad( dot(grad(v)('-'),n('+')) ), n('+') ) ) , n('+'))

							J_n   = eta_n *  inner( jump_u_n , jump_v_n  )  * dS


							R_n = I_adv + I_diff + I_Nit + J_n

						else:

							R_n = I_adv + I_diff + I_Nit

						# solve the system
						spline.iteratedDivFreeSolve(R_n,u_hat,v_hat,penalty=Constant(1e6))
						u_sol = spline.projectScalarOntoLinears(u[0])
						v_sol = spline.projectScalarOntoLinears(u[1])

						if Re in [3200,5000,7500,10000]:
							ufile = File('SS_results/ss_Re='+str(Re)+'_p=' + str(p) + '_m=' + str(resol)  + '/u.pvd')
							vfile = File('SS_results/ss_Re='+str(Re)+'_p=' + str(p) + '_m=' + str(resol)  + '/v.pvd')

							u_sol.rename('u_sol','u')
							v_sol.rename('v_sol','v')

							ufile << u_sol
							vfile << v_sol
			if case == 1:
				print('case1: unstabilized')
				for Re in [10,100,250,400,600,800,1000,1200,1500,2000,2600,3200,3800,4400,5000,5600,6200,6800,7500,8200,8500,9000,10000]:
					nu = 1.0/Re
					if Re == 10:
						print("Starting Re = " + str(Re) + '_k = '  + str(p) + '_m = ' + str(resol))
						# always redefine, in case of unwanted initialization	        
						u_hat = Function(spline.V)
						v_hat = TestFunction(spline.V)

						u = spline.pushforward(u_hat)
						v = spline.pushforward(v_hat)

						# weak forms
						I_adv  = inner(grad(u)*u,v)*dx
						I_diff = 2.0*nu*inner(nablas(u),nablas(v))*dx

						# top edge, enforce u[0] = 1
						ut = as_vector([u[0],0.0])
						vt = as_vector([v[0],0.0])

						It_c = -2.0*nu*inner(nablas(ut)*n,vt) * ds(1)
						It_s = -2.0*nu*inner(nablas(vt)*n,ut-uD) * ds(1)
						It_p = Cb*nu/h*inner(vt,ut-uD) * ds(1)

						# bottom edge, u[0] = 0
						ub = as_vector([u[0],0.0])
						vb = as_vector([v[0],0.0])

						Ib_c = -2.0*nu*inner(nablas(ub)*n,vb) * ds(2)
						Ib_s = -2.0*nu*inner(nablas(vb)*n,ub) * ds(2)
						Ib_p = Cb*nu/h*inner(vb,ub) * ds(2)

						# left/right, u[1] = 0
						ulr = as_vector([0.0,u[1]])
						vlr = as_vector([0.0,v[1]])

						Ilr_c = -2.0*nu*inner(nablas(ulr)*n,vlr) * ( ds(3) + ds(4) )
						Ilr_s = -2.0*nu*inner(nablas(vlr)*n,ulr) * ( ds(3) + ds(4) )
						Ilr_p = Cb*nu/h*inner(vlr,ulr) * ( ds(3) + ds(4) )


						I_c = It_c + Ib_c + Ilr_c # consistency
						I_s = It_s + Ib_s + Ilr_s # symmetry
						I_p = It_p + Ib_p + Ilr_p # penalty

						I_Nit = I_c + I_s + I_p

						# solve it
						R     = I_adv + I_diff + I_Nit # residual, not stabilization
						spline.iteratedDivFreeSolve(R,u_hat,v_hat,penalty=Constant(1e6))
					else:
						print("Starting Re = " + str(Re) + '_k = '  + str(p) + '_m = ' + str(resol))
						# div-project the last result as new ic
						u_IG  = as_vector([u[0],u[1]])
						u_hat = spline.divFreeProject(u_IG)
						v_hat = TestFunction(spline.V)

						u,v   = spline.pushforward(u_hat), spline.pushforward(v_hat)


						I_adv  = inner(grad(u)*u,v)*dx
						I_diff = 2.0*nu*inner(nablas(u),nablas(v))*dx

						# top edge, enforce u[0] = 1
						ut = as_vector([u[0],0.0])
						vt = as_vector([v[0],0.0])

						It_c = -2.0*nu*inner(nablas(ut)*n,vt) * ds(1)
						It_s = -2.0*nu*inner(nablas(vt)*n,ut-uD) * ds(1)
						It_p = Cb*nu/h*inner(vt,ut-uD) * ds(1)

						# bottom edge, u[0] = 0
						ub = as_vector([u[0],0.0])
						vb = as_vector([v[0],0.0])

						Ib_c = -2.0*nu*inner(nablas(ub)*n,vb) * ds(2)
						Ib_s = -2.0*nu*inner(nablas(vb)*n,ub) * ds(2)
						Ib_p = Cb*nu/h*inner(vb,ub) * ds(2)

						# left/right, u[1] = 0
						ulr = as_vector([0.0,u[1]])
						vlr = as_vector([0.0,v[1]])

						Ilr_c = -2.0*nu*inner(nablas(ulr)*n,vlr) * ( ds(3) + ds(4) )
						Ilr_s = -2.0*nu*inner(nablas(vlr)*n,ulr) * ( ds(3) + ds(4) )
						Ilr_p = Cb*nu/h*inner(vlr,ulr) * ( ds(3) + ds(4) )


						I_c = It_c + Ib_c + Ilr_c
						I_s = It_s + Ib_s + Ilr_s
						I_p = It_p + Ib_p + Ilr_p

						I_Nit = I_c + I_s + I_p

						R_n = I_adv + I_diff + I_Nit

						# solve the system
						spline.iteratedDivFreeSolve(R_n,u_hat,v_hat,penalty=Constant(1e6))
						u_sol = spline.projectScalarOntoLinears(u[0])
						v_sol = spline.projectScalarOntoLinears(u[1])

						if Re in [3200,5000,7500,10000]:
							ufile = File('NM_results/nm_Re='+str(Re)+'_p=' + str(p) + '_m=' + str(resol)  + '/u.pvd')
							vfile = File('NM_results/nm_Re='+str(Re)+'_p=' + str(p) + '_m=' + str(resol)  + '/v.pvd')

							u_sol.rename('u_sol','u')
							v_sol.rename('v_sol','v')

							ufile << u_sol
							vfile << v_sol





















