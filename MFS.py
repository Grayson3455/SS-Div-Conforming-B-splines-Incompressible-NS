# This is an example code of steady incompressible Navier-Stokes equation
# for the convergence test of the skeleton-stabilizatiion solver, example chosen
# from the A Buffa, C de Falco, and R Va ́zquez. Isogeometric analysis: Stable elements for the 2D Stokes equation.
# International Journal for Numerical Methods in Fluids, 65:1407–1422,20–30, 2011.

#The code here based on the tutorial code provided by Prof. Kamensky's tIGAr package

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

for p in [1,2,3]:   # k'
    QUAD_DEG =  10  # overkill for convergence
    mode = 'w'
    mesh = [4,8,16,32,64,128]  # keep mesh size fixed
    form_compiler_parameters = {"quadrature_degree":QUAD_DEG}  # minimum quadrature_degree to ensure exact Integration

    for resol in mesh:
        # start to define B-spline mesh
        kv = [ uniformKnots(p,0.0,1.0,resol,False) for i in range(0,2) ]  # uniform knots in each direction

        controlMesh = ExplicitBSplineControlMesh([p,p],kv)                # B-spline control mesh

        fieldList = generateFieldsCompat(controlMesh,"RT",[p,p])          # soultion field basis over the control mesh, rt type
        splineGenerator = FieldListSpline(controlMesh,fieldList)

        # Strong normal BC only, use Nitsche for tangential
        for field in [0,1]: # u,v
            scalarSpline = splineGenerator.getFieldSpline(field)
            for side in [0,1]: # 1d: points, 2d: edges, 3d: faces
                sideDofs = scalarSpline.getSideDofs(field,side)
                splineGenerator.addZeroDofs(field,sideDofs)

        # extraction and parameters
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


        #exact solutions

        u_ex = 2.0*exp(x[0])*(x[0]-1.0)*(x[0]-1.0)*x[0]*x[0]*(x[1]*x[1] - x[1])*(2.0*x[1] - 1.0)
        v_ex = -exp(x[0])*(x[0]-1.0)*x[0]*(x[0]*(3.0 + x[0]) - 2.0)*(x[1]-1.0)*(x[1]-1.0)*x[1]*x[1]
        p_ex = (-424.0+156.0*exp(1.0)+(x[1]*x[1]-x[1]) * (-456.0 + exp(x[0])*(456.0 + x[0]*x[0]* ( \
        	   228.0 - 5.0*(x[1]*x[1] - x[1] ))  +2.0*x[0]* (-228.0 + (x[1]*x[1] - x[1]) ) + \
        	   2.0*x[0]*x[0]*x[0]* (-36.0 + (x[1]*x[1] - x[1])) + x[0]*x[0]*x[0]*x[0]*(12.0 + x[1]*x[1]-x[1]))))

        u_exact = as_vector((u_ex,v_ex))


        # use iterative solver instead
        # MAX_KSP_IT = 5000
        # spline.linearSolver = PETScKrylovSolver("gmres","jacobi")
        # spline.linearSolver.parameters["relative_tolerance"] = 1e-2
        # spline.linearSolver.parameters["error_on_nonconvergence"] = False
        # spline.linearSolver.parameters["maximum_iterations"] = MAX_KSP_IT
        # spline.linearSolver.ksp().setGMRESRestart(MAX_KSP_IT)
        spline.relativeTolerance = 1e-8
        spline.maxIters = 50
        div_pent = 1e4


        for Re in [10]:
            nu = 1.0/Re

            #outFile1 = open('nm_p=' + str(p) + '_Re='+str(Re)+ '_m=' + str(resol) +'.csv', mode)
            outFile2 = open('error/ss_p=' + str(p) + '_Re='+str(Re)+ '_m=' + str(resol) +'.csv', mode)

            # solving
            for case in [2]:
                if case == 1: # no model

                    #-------------------------------------------------------#
                    #         NO NITSCHE'S FOR UNSTABILIZED FOR NOW         #

                    # reinitialization, otherwise last results will be initialized
                    u_hat = Function(spline.V)
                    v_hat = TestFunction(spline.V)
                    u = spline.pushforward(u_hat)
                    v = spline.pushforward(v_hat)


                    # start building the form
                    I_adv  = inner(grad(u)*u,v)*dx
                    I_diff = 2.0*nu*inner(nablas(u),nablas(v))*dx
                    f = inner(grad(u_exact)*u_exact,v)*dx - 2.0*nu*inner(div(nablas(u_exact)),v)*dx + inner(grad(p_ex),v)*dx


                    print('case='+str(case)+ ', p=' + str(p) + ', m=' + str(resol)+ ', Re=' + str(Re))


                    R     = I_adv + I_diff - f

                    spline.iteratedDivFreeSolve(R,u_hat,v_hat,penalty=Constant(div_pent)) # iterative penalty solver

                    # save the solution
                    u_sol = spline.projectScalarOntoLinears(u[0])
                    v_sol = spline.projectScalarOntoLinears(u[1])

                    ufile = File('nm_results/Re='+str(Re)+'_p=' + str(p) + '_m=' + str(resol) + '/u.pvd')
                    vfile = File('nm_results/Re='+str(Re)+'_p=' + str(p) + '_m=' + str(resol) + '/v.pvd')

                    u_sol.rename('u_sol','u')
                    v_sol.rename('v_sol','v')

                    ufile << u_sol
                    vfile << v_sol

                    # error estimate
                    error = u - u_exact
                    L2_err      = sqrt(assemble(inner(error,error)*dx))
                    Semi_H1_err = sqrt(assemble(inner(grad(error), grad(error))*dx))

                    #save the data to files
                    outFile1.write(str(L2_err) + "," + str(Semi_H1_err)  + "\n")


                if case == 2: # skeleton stabilized

                    # reinitialization, otherwise last results will be initialized
                    u_hat = Function(spline.V)
                    v_hat = TestFunction(spline.V)
                    u = spline.pushforward(u_hat)
                    v = spline.pushforward(v_hat)


                    # start building the form
                    I_adv  = inner(grad(u)*u,v)*dx
                    I_diff = 2.0*nu*inner(nablas(u),nablas(v))*dx

                    f = inner(grad(u_exact)*u_exact,v)*dx - 2.0*nu*inner(div(nablas(u_exact)),v)*dx + inner(grad(p_ex),v)*dx

                    # start to build Nitsche term
                    # top edge
                    ut = as_vector([u[0],0.0])
                    vt = as_vector([v[0],0.0])

                    It_c = -2.0*nu*inner(nablas(ut)*n,vt) * ds(1)
                    It_s = -2.0*nu*inner(nablas(vt)*n,ut) * ds(1)
                    It_p = Cb*nu/h*inner(vt,ut) * ds(1)

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

                    # Skeleton stabilization term
                    def Min(a, b): return (a+b-abs(a-b))/Constant(2.0)
                    gamma = 1*10**(-p-1)
                    print('case='+str(case)+ ', p=' + str(p) + ', m=' + str(resol)+', gamma='+str(gamma)+ ', Re=' + str(Re))

                    alpha = p - 1    # regularity
                    Re_e  = avg( (dot(u,u) + 1e-32)**0.5)*h_avg/nu  # local Reynolds number
                    #Note: a potential bug when sqrt is taken, had to add a very small positive number to avoid that 
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


                    R     = I_adv + I_diff - f + J_n + I_Nit

                    spline.iteratedDivFreeSolve(R,u_hat,v_hat,penalty=Constant(div_pent)) # iterative penalty method


                    # save the solution
                    u_sol = spline.projectScalarOntoLinears(u[0])
                    v_sol = spline.projectScalarOntoLinears(u[1])

                    ufile = File('ss_results/Re='+str(Re)+'_p=' + str(p) + '_m=' + str(resol) + '/u.pvd')
                    vfile = File('ss_results/Re='+str(Re)+'_p=' + str(p) + '_m=' + str(resol) + '/v.pvd')

                    u_sol.rename('u_sol','u')
                    v_sol.rename('v_sol','v')

                    ufile << u_sol
                    vfile << v_sol

                    # error estimate
                    error = u - u_exact
                    L2_err      = sqrt(assemble(inner(error,error)*dx))
                    Semi_H1_err = sqrt(assemble(inner(grad(error), grad(error))*dx))

                    # save the data to files
                    outFile2.write(str(L2_err) + "," + str(Semi_H1_err) + "\n")


            #outFile1.close()
            outFile2.close()






