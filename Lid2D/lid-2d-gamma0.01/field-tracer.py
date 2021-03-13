# trace generated using paraview version 5.9.0-RC3

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

for p in [1,2,3]:
	for m in [16,32,128]:
		for Re in [7500,10000]:
			# create a new 'PVD Reader'
			upvd = PVDReader(registrationName='u.pvd', FileName='/Users/gxt/Desktop/My_Computer/open-sources/SS-Div-Conforming-B-splines-Incompressible-NS/Lid2D/Lid-2d-gamma0.01/SS_results/ss_Re='+str(Re)+'_p='+str(p)+'_m='+str(m)+'/u.pvd')
			upvd.PointArrays = ['u_sol']

			# get active view
			renderView1 = GetActiveViewOrCreate('RenderView')

			# show data in view
			upvdDisplay = Show(upvd, renderView1, 'UnstructuredGridRepresentation')

			# get color transfer function/color map for 'u_sol'
			u_solLUT = GetColorTransferFunction('u_sol')

			# get opacity transfer function/opacity map for 'u_sol'
			u_solPWF = GetOpacityTransferFunction('u_sol')

			# trace defaults for the display properties.
			upvdDisplay.Representation = 'Surface'
			upvdDisplay.ColorArrayName = ['POINTS', 'u_sol']
			upvdDisplay.LookupTable = u_solLUT
			upvdDisplay.SelectTCoordArray = 'None'
			upvdDisplay.SelectNormalArray = 'None'
			upvdDisplay.SelectTangentArray = 'None'
			upvdDisplay.OSPRayScaleArray = 'u_sol'
			upvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
			upvdDisplay.SelectOrientationVectors = 'None'
			upvdDisplay.ScaleFactor = 0.1
			upvdDisplay.SelectScaleArray = 'u_sol'
			upvdDisplay.GlyphType = 'Arrow'
			upvdDisplay.GlyphTableIndexArray = 'u_sol'
			upvdDisplay.GaussianRadius = 0.005
			upvdDisplay.SetScaleArray = ['POINTS', 'u_sol']
			upvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
			upvdDisplay.OpacityArray = ['POINTS', 'u_sol']
			upvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
			upvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
			upvdDisplay.PolarAxes = 'PolarAxesRepresentation'
			upvdDisplay.ScalarOpacityFunction = u_solPWF
			upvdDisplay.ScalarOpacityUnitDistance = 0.22272467953508482
			upvdDisplay.OpacityArrayName = ['POINTS', 'u_sol']

			# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
			upvdDisplay.ScaleTransferFunction.Points = [-0.3683701736346532, 0.0, 0.5, 0.0, 0.8195239102446503, 1.0, 0.5, 0.0]

			# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
			upvdDisplay.OpacityTransferFunction.Points = [-0.3683701736346532, 0.0, 0.5, 0.0, 0.8195239102446503, 1.0, 0.5, 0.0]

			# reset view to fit data
			renderView1.ResetCamera()

			#changing interaction mode based on data extents
			renderView1.InteractionMode = '2D'
			renderView1.CameraPosition = [0.5, 0.5, 10000.0]
			renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]

			# get the material library
			materialLibrary1 = GetMaterialLibrary()

			# show color bar/color legend
			upvdDisplay.SetScalarBarVisibility(renderView1, True)

			# update the view to ensure updated data information
			renderView1.Update()

			# save data
			SaveData('/Users/gxt/Desktop/My_Computer/open-sources/SS-Div-Conforming-B-splines-Incompressible-NS/Lid2D/Lid-2d-gamma0.01/SS_results/ss_Re='+str(Re)+'_p='+str(p)+'_m='+str(m)+'/u-xy.csv', proxy=upvd, PointDataArrays=['u_sol'])

			# create a new 'PVD Reader'
			vpvd = PVDReader(registrationName='v.pvd', FileName='/Users/gxt/Desktop/My_Computer/open-sources/SS-Div-Conforming-B-splines-Incompressible-NS/Lid2D/Lid-2d-gamma0.01/SS_results/ss_Re='+str(Re)+'_p='+str(p)+'_m='+str(m)+'/v.pvd')
			vpvd.PointArrays = ['v_sol']

			# show data in view
			vpvdDisplay = Show(vpvd, renderView1, 'UnstructuredGridRepresentation')

			# get color transfer function/color map for 'v_sol'
			v_solLUT = GetColorTransferFunction('v_sol')

			# get opacity transfer function/opacity map for 'v_sol'
			v_solPWF = GetOpacityTransferFunction('v_sol')

			# trace defaults for the display properties.
			vpvdDisplay.Representation = 'Surface'
			vpvdDisplay.ColorArrayName = ['POINTS', 'v_sol']
			vpvdDisplay.LookupTable = v_solLUT
			vpvdDisplay.SelectTCoordArray = 'None'
			vpvdDisplay.SelectNormalArray = 'None'
			vpvdDisplay.SelectTangentArray = 'None'
			vpvdDisplay.OSPRayScaleArray = 'v_sol'
			vpvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
			vpvdDisplay.SelectOrientationVectors = 'None'
			vpvdDisplay.ScaleFactor = 0.1
			vpvdDisplay.SelectScaleArray = 'v_sol'
			vpvdDisplay.GlyphType = 'Arrow'
			vpvdDisplay.GlyphTableIndexArray = 'v_sol'
			vpvdDisplay.GaussianRadius = 0.005
			vpvdDisplay.SetScaleArray = ['POINTS', 'v_sol']
			vpvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
			vpvdDisplay.OpacityArray = ['POINTS', 'v_sol']
			vpvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
			vpvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
			vpvdDisplay.PolarAxes = 'PolarAxesRepresentation'
			vpvdDisplay.ScalarOpacityFunction = v_solPWF
			vpvdDisplay.ScalarOpacityUnitDistance = 0.22272467953508482
			vpvdDisplay.OpacityArrayName = ['POINTS', 'v_sol']

			# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
			vpvdDisplay.ScaleTransferFunction.Points = [-0.7448443833059003, 0.0, 0.5, 0.0, 0.3825274162646751, 1.0, 0.5, 0.0]

			# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
			vpvdDisplay.OpacityTransferFunction.Points = [-0.7448443833059003, 0.0, 0.5, 0.0, 0.3825274162646751, 1.0, 0.5, 0.0]

			# show color bar/color legend
			vpvdDisplay.SetScalarBarVisibility(renderView1, True)

			# update the view to ensure updated data information
			renderView1.Update()

			# hide data in view
			Hide(upvd, renderView1)

			# save data
			SaveData('/Users/gxt/Desktop/My_Computer/open-sources/SS-Div-Conforming-B-splines-Incompressible-NS/Lid2D/Lid-2d-gamma0.01/SS_results/ss_Re='+str(Re)+'_p='+str(p)+'_m='+str(m)+'/v-xy.csv', proxy=vpvd, PointDataArrays=['v_sol'])

			ResetSession()

			# get active view
			renderView1_1 = GetActiveViewOrCreate('RenderView')

			#================================================================
			# addendum: following script captures some of the application
			# state to faithfully reproduce the visualization during playback
			#================================================================

			# get layout
			layout1 = GetLayout()

			#--------------------------------
			# saving layout sizes for layouts

			# layout/tab size in pixels
			layout1.SetSize(2594, 1608)

			#--------------------------------------------
			# uncomment the following to render all views
			# RenderAllViews()
			# alternatively, if you want to write images, you can use SaveScreenshot(...).