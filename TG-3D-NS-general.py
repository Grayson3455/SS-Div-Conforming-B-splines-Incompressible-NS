# this is an example code of unsteady incompressible Navier-Stokes equation
# based on the online tutorial of DOLFIN fenics, numerical example refers to
# the 3D T-G vortex flow, asympotic solution given by the PhD thesis of professor
# J.A. Evans and other online resources. The code here based on the tutorial code
# provided by Prof. Kamensky's tIGAr package.

# Solving this problem requires a compatible Bspline discretization(DIV-free) with
# iterative penalty method and generalized-alpha method.


# data to be stored, velocities, vorticities and Q-criterion


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

# unSymmetrizied gradient operator:
def nablaus(u):
    return 0.5*(grad(u) - grad(u).T)

# Trace operator
def Tr(A):
    return A[0,0] + A[1,1] + A[2,2]

# problem set up
for p in [1]: # k'
    QUAD_DEG =  p + 2 # highest deg in RT: p+1
    mesh = [64]
    form_compiler_parameters = {"quadrature_degree":QUAD_DEG}  # minimum quadrature_degree to ensure exact Integration
    for resol in mesh:
        # start to define B-spline mesh
        kv = [ uniformKnots(p,0.0,math.pi,resol,False) for i in range(0,3) ]  # uniform knots in each direction

        controlMesh = ExplicitBSplineControlMesh([p,p,p],kv)                # B-spline control mesh

        fieldList = generateFieldsCompat(controlMesh,"RT",[p,p,p])          # soultion field basis over the control mesh, rt type
        splineGenerator = FieldListSpline(controlMesh,fieldList)

        # Strong normal BC
        for field in [0,1,2]: # u,v,w
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


        # Time integration parameters
        T   = 10.0
        dT  = math.pi/resol*0.25
        Num = int(T/dT)

        # Initial conditions
        x  = spline.spatialCoordinates()
        u0 = sin(x[0])*cos(x[1])*cos(x[2])
        v0 = -cos(x[0])*sin(x[1])*cos(x[2])
        U0 = as_vector([u0,v0,Constant(0.0)])

        # For 3D computations, use an iterative solver. From tIGAr example
        spline.linearSolver = PETScKrylovSolver("gmres","jacobi")
        spline.linearSolver.parameters["relative_tolerance"] = 1e-1
        # Linear solver sometimes fails to converge, but convergence of nonlinear
        # iteration is still enforced to within spline.relativeTolerance.
        spline.linearSolver.parameters["error_on_nonconvergence"] = False
        spline.linearSolver.parameters["maximum_iterations"] = 1000
        spline.linearSolver.ksp().setGMRESRestart(1000)
        spline.linearSolver.ksp().setNormType(PETSc.KSP.NormType.UNPRECONDITIONED)
        spline.relativeTolerance = 1e-3

         # construction of the predictor n-level
        u_hat        = Function(spline.V)            # some placeholder for t_{n}
        v_hat        = TestFunction(spline.V)        # test function in para space
        u_old_hat    = spline.divFreeProject(U0)     # solution for t = 0, at RT space, will be replaced at each time step
        udot_old_hat = Function(spline.V)            # zeros as udot
        # why zero damping
        Int          = GeneralizedAlphaIntegrator(0.5,dT,u_hat,(u_old_hat, udot_old_hat)) # generalized-alpha integator

        # construction of the predictor \alpha-level
        u_alp    = spline.pushforward(Int.x_alpha())
        udot_alp = spline.pushforward(Int.xdot_alpha())
        v        = spline.pushforward(v_hat)

        # variational forms
        Re = 1600
        nu = 1.0/Re

        I_t =  inner(udot_alp,v)*dx                       # time derivative
        I_c =  inner(grad(u_alp)*u_alp,v)*dx              # advection
        I_d =  2.0*nu*inner(nablas(u_alp), nablas(v))*dx  # diffusion

        # Normal skeleton-stabilization

        if p == 1:
            gamma = 0.01
        elif p == 2:
            gamma = 0.001

        def Min(a, b): return (a+b-abs(a-b))/Constant(2.0)

        alpha = p - 1    # regularity

        Re_e  = avg(dot(u_alp,u_alp)**0.5)*h_avg/nu # local Reynolds number

        eta_n = gamma * Min( Constant(1.0), Re_e ) * h_avg**( 2*alpha + 2) * sqrt( ( dot(u_alp('+'),n('+')) )**2 )

        if p == 1:

            jump_u_n  = dot( ( grad(u_alp)('+') - grad(u_alp)('-') ) , n('+'))
            jump_v_n  = dot( ( grad(v)('+') - grad(v)('-') ) , n('+'))


        if p == 2:

            jump_u_n  = dot( grad( dot(grad(u_alp)('+'),n('+')) )  - grad( dot(grad(u_alp)('-'),n('+')) ) , n('+') )
            jump_v_n  = dot( grad( dot(grad(v)('+'),n('+')) )  - grad( dot(grad(v)('-'),n('+')) ) , n('+') )

        J_n   = eta_n *  inner( jump_u_n , jump_v_n  )  * dS

        R   = I_t + I_c + I_d + J_n     # residual

        #------------save the results-----------------#
        qfile  = File('results/NM-kprime='+str(p)+'-m='+str(resol)+'/q.pvd')     # Q-criterion

        ufile  = File('results/NM-kprime='+str(p)+'-m='+str(resol)+'/u.pvd')     # u - velocity
        vfile  = File('results/NM-kprime='+str(p)+'-m='+str(resol)+'/v.pvd')     # v - velocity
        wfile  = File('results/NM-kprime='+str(p)+'-m='+str(resol)+'/w.pvd')     # w - velocity


        xvotfile = File('results/NM-kprime='+str(p)+'-m='+str(resol)+'/X-vorticity.pvd')  # x-vorticity
        yvotfile = File('results/NM-kprime='+str(p)+'-m='+str(resol)+'/Y-vorticity.pvd')  # y-vorticity
        zvotfile = File('results/NM-kprime='+str(p)+'-m='+str(resol)+'/Z-vorticity.pvd')  # z-vorticity

        w = Function(spline.V) # placeholder for function in iterative penalty method

        # Time integration
        for i in range(0,Num):

            if (mpirank==0):
                print(" t = " + str(Int.t))
            # solve
            spline.iteratedDivFreeSolve(R,u_hat,v_hat,penalty=Constant(1e4), w = w) # solve the system, save the solution in u_alp

            vor = curl(u_alp) # compute vorticity

            u_sol = spline.projectScalarOntoLinears(u_alp[0])
            v_sol = spline.projectScalarOntoLinears(u_alp[1])
            w_sol = spline.projectScalarOntoLinears(u_alp[2])

            vorX_sol = spline.projectScalarOntoLinears(vor[0])
            vorY_sol = spline.projectScalarOntoLinears(vor[1])
            vorZ_sol = spline.projectScalarOntoLinears(vor[2])

            Omega  = nablaus(u_alp)
            Strain = nablas(u_alp)

            Q      = 0.5*( Tr(Omega*Omega.T) - Tr(Strain*Strain.T) )
            Q_sol  = spline.projectScalarOntoLinears(Q)

            if i%20==0:
            	u_sol.rename("u_sol", "u")
            	v_sol.rename("v_sol", "v")
            	w_sol.rename("w_sol", "w")

            	vorX_sol.rename("X-vorticity", "X-vorticity")
            	vorY_sol.rename("Y-vorticity", "Y-vorticity")
            	vorZ_sol.rename("Z-vorticity", "Z-vorticity")

            	Q_sol.rename("Q", "Q")
            	qfile << Q_sol

            	ufile << u_sol
            	vfile << v_sol
            	wfile << w_sol

            	xvotfile << vorX_sol
            	yvotfile << vorY_sol
            	zvotfile << vorZ_sol


            # calculate dissipation rate
            dissipationRate  = assemble( (2.0*nu/(math.pi)**3)*inner(nablas(u_alp),nablas(u_alp))*dx )
            TKE              = 1.0/2 * (1.0/(pi)**3) *assemble(inner(u_alp,u_alp)*dx)
            MdissipationRate = (1.0/(pi)**3) * assemble( eta_n *  inner( jump_u_n , jump_u_n  )  * dS )

            if (mpirank==0):
                mode = 'a'
                if i == 0:
                    mode = 'w'
                DisFile = open('NS-kprime='+str(p)+'-m='+str(resol)+'-gamma='+str(gamma)+'.csv',mode)
                DisFile.write(str(Int.t) + "," + str(dissipationRate) + "," + str(MdissipationRate) + "," + str(TKE) + "\n")
                DisFile.close()
            Int.advance() # next step
