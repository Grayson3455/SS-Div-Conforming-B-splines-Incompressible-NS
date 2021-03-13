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
Re = 10
nu = 1.0/Re

for p in [1]:   # k'
    QUAD_DEG =  10  # overkill for convergence
    mode = 'w'
    mesh = [16]  # keep mesh size fixed
    form_compiler_parameters = {"quadrature_degree":QUAD_DEG}  # minimum quadrature_degree to ensure exact Integration

    for resol in mesh:
        # start to define B-spline mesh
        kv = [ uniformKnots(p,0.0,1.0,resol,False) for i in range(0,2) ]  # uniform knots in each direction

        controlMesh = ExplicitBSplineControlMesh([p,p],kv)                # B-spline control mesh

        fieldList = generateFieldsCompat(controlMesh,"RT",[p,p])          # soultion field basis over the control mesh, rt type
        splineGenerator = FieldListSpline(controlMesh,fieldList)

        # Strong homogeneous dirichelt BC everywhere
        for field in range(0,2):
            scalarSpline = splineGenerator.getFieldSpline(field)
            for side in range(0,2):
                for direction in range(0,2):
                    sideDofs = scalarSpline.getSideDofs(direction,side)
                    splineGenerator.addZeroDofs(field,sideDofs)

        # extraction and parameters
        spline = ExtractedBSplineRT(splineGenerator,QUAD_DEG)
        n      = FacetNormal(spline.mesh)
        h      = CellDiameter(spline.mesh)
        h_avg  = 0.5*(h('+') + h('-'))
        dS  =  dS(metadata = {"quadrature_degree":QUAD_DEG}) #interior facets
        dx  =  dx(metadata = {"quadrature_degree":QUAD_DEG}) #interior
        x = spline.spatialCoordinates()

        #exact solutions

        u_ex = 2.0*exp(x[0])*(x[0]-1.0)*(x[0]-1.0)*x[0]*x[0]*(x[1]*x[1] - x[1])*(2.0*x[1] - 1.0)
        v_ex = -exp(x[0])*(x[0]-1.0)*x[0]*(x[0]*(3.0 + x[0]) - 2.0)*(x[1]-1.0)*(x[1]-1.0)*x[1]*x[1]
        p_ex1 = (-424.0+156.0*exp(1.0)+(x[1]*x[1]-x[1]) * (-456.0 + exp(x[0])*(456.0 + x[0]*x[0]* ( \
        	   228.0 - 5.0*(x[1]*x[1] - x[1] ))  +2.0*x[0]* (-228.0 + (x[1]*x[1] - x[1]) ) + \
        	   2.0*x[0]*x[0]*x[0]* (-36.0 + (x[1]*x[1] - x[1])) + x[0]*x[0]*x[0]*x[0]*(12.0 + x[1]*x[1]-x[1]))))
        p_ex2 = 10000*x[0]*x[1] 

        u_exact = as_vector((u_ex,v_ex))

        for p_case in [1,2]:
            if p_case == 1: # apply term

                outFile2 = open('t1_ss_p=' + str(p) + '_Re='+str(Re)+ '_m=' + str(resol) +'.csv', mode)

                # reinitialization, otherwise last results will be initialized
                u_hat = Function(spline.V)
                v_hat = TestFunction(spline.V)
                u = spline.pushforward(u_hat)
                v = spline.pushforward(v_hat)


                # start building the form
                I_adv  = inner(grad(u)*u,v)*dx
                I_diff = 2.0*nu*inner(nablas(u),nablas(v))*dx

                f = inner(grad(u_exact)*u_exact,v)*dx - 2.0*nu*inner(div(nablas(u_exact)),v)*dx + inner(grad(p_ex1),v)*dx 


                def Min(a, b): return (a+b-abs(a-b))/Constant(2.0)
                gamma = 1*10**(-p-1)
                print('case='+str(p_case)+ ', p=' + str(p) + ', m=' + str(resol)+', gamma='+str(gamma)+ ', Re=' + str(Re))

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


                R     = I_adv + I_diff - f + J_n 

                spline.iteratedDivFreeSolve(R,u_hat,v_hat,penalty=Constant(1e4)) # iterative penalty method


                # save the solution
                u_sol = spline.projectScalarOntoLinears(u[0])
                v_sol = spline.projectScalarOntoLinears(u[1])

                ufile = File('t1_ss_results/Re='+str(Re)+'_p=' + str(p) + '_m=' + str(resol) + '/u.pvd')
                vfile = File('t1_ss_results/Re='+str(Re)+'_p=' + str(p) + '_m=' + str(resol) + '/v.pvd')

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
                outFile2.close()


            else:  # apply term 2

                outFile3 = open('t2_ss_p=' + str(p) + '_Re='+str(Re)+ '_m=' + str(resol) +'.csv', mode)

                # reinitialization, otherwise last results will be initialized
                u_hat = Function(spline.V)
                v_hat = TestFunction(spline.V)
                u = spline.pushforward(u_hat)
                v = spline.pushforward(v_hat)


                # start building the form
                I_adv  = inner(grad(u)*u,v)*dx
                I_diff = 2.0*nu*inner(nablas(u),nablas(v))*dx

                f = inner(grad(u_exact)*u_exact,v)*dx - 2.0*nu*inner(div(nablas(u_exact)),v)*dx + inner(grad(p_ex2),v)*dx 


                def Min(a, b): return (a+b-abs(a-b))/Constant(2.0)
                gamma = 1*10**(-p-1)
                print('case='+str(p_case)+ ', p=' + str(p) + ', m=' + str(resol)+', gamma='+str(gamma)+ ', Re=' + str(Re))

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


                R     = I_adv + I_diff - f + J_n 

                spline.iteratedDivFreeSolve(R,u_hat,v_hat,penalty=Constant(1e4)) # iterative penalty method


                # save the solution
                u_sol = spline.projectScalarOntoLinears(u[0])
                v_sol = spline.projectScalarOntoLinears(u[1])

                ufile = File('t2_ss_results/Re='+str(Re)+'_p=' + str(p) + '_m=' + str(resol) + '/u.pvd')
                vfile = File('t2_ss_results/Re='+str(Re)+'_p=' + str(p) + '_m=' + str(resol) + '/v.pvd')

                u_sol.rename('u_sol','u')
                v_sol.rename('v_sol','v')

                ufile << u_sol
                vfile << v_sol

                # error estimate
                error = u - u_exact
                L2_err      = sqrt(assemble(inner(error,error)*dx))
                Semi_H1_err = sqrt(assemble(inner(grad(error), grad(error))*dx))

                # save the data to files
                outFile3.write(str(L2_err) + "," + str(Semi_H1_err) + "\n")
                outFile3.close()
                