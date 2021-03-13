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

			# create a new 'Plot Over Line'
			plotOverLine1 = PlotOverLine(registrationName='PlotOverLine1', Input=upvd,
			    Source='Line')

			# init the 'Line' selected for 'Source'
			plotOverLine1.Source.Point1 = [0.5, 0.0, 0.0]
			plotOverLine1.Source.Point2 = [0.5, 1.0, 0.0]

			# show data in view
			plotOverLine1Display = Show(plotOverLine1, renderView1, 'GeometryRepresentation')

			# trace defaults for the display properties.
			plotOverLine1Display.Representation = 'Surface'
			plotOverLine1Display.ColorArrayName = ['POINTS', 'u_sol']
			plotOverLine1Display.LookupTable = u_solLUT
			plotOverLine1Display.SelectTCoordArray = 'None'
			plotOverLine1Display.SelectNormalArray = 'None'
			plotOverLine1Display.SelectTangentArray = 'None'
			plotOverLine1Display.OSPRayScaleArray = 'u_sol'
			plotOverLine1Display.OSPRayScaleFunction = 'PiecewiseFunction'
			plotOverLine1Display.SelectOrientationVectors = 'None'
			plotOverLine1Display.ScaleFactor = 0.1
			plotOverLine1Display.SelectScaleArray = 'u_sol'
			plotOverLine1Display.GlyphType = 'Arrow'
			plotOverLine1Display.GlyphTableIndexArray = 'u_sol'
			plotOverLine1Display.GaussianRadius = 0.005
			plotOverLine1Display.SetScaleArray = ['POINTS', 'u_sol']
			plotOverLine1Display.ScaleTransferFunction = 'PiecewiseFunction'
			plotOverLine1Display.OpacityArray = ['POINTS', 'u_sol']
			plotOverLine1Display.OpacityTransferFunction = 'PiecewiseFunction'
			plotOverLine1Display.DataAxesGrid = 'GridAxesRepresentation'
			plotOverLine1Display.PolarAxes = 'PolarAxesRepresentation'

			# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
			plotOverLine1Display.ScaleTransferFunction.Points = [-0.36260967817893625, 0.0, 0.5, 0.0, 0.738922618749873, 1.0, 0.5, 0.0]

			# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
			plotOverLine1Display.OpacityTransferFunction.Points = [-0.36260967817893625, 0.0, 0.5, 0.0, 0.738922618749873, 1.0, 0.5, 0.0]

			# Create a new 'Line Chart View'
			lineChartView1 = CreateView('XYChartView')

			# show data in view
			plotOverLine1Display_1 = Show(plotOverLine1, lineChartView1, 'XYChartRepresentation')

			# trace defaults for the display properties.
			plotOverLine1Display_1.CompositeDataSetIndex = [0]
			plotOverLine1Display_1.UseIndexForXAxis = 0
			plotOverLine1Display_1.XArrayName = 'arc_length'
			plotOverLine1Display_1.SeriesVisibility = ['u_sol']
			plotOverLine1Display_1.SeriesLabel = ['arc_length', 'arc_length', 'u_sol', 'u_sol', 'vtkValidPointMask', 'vtkValidPointMask', 'Points_X', 'Points_X', 'Points_Y', 'Points_Y', 'Points_Z', 'Points_Z', 'Points_Magnitude', 'Points_Magnitude']
			plotOverLine1Display_1.SeriesColor = ['arc_length', '0', '0', '0', 'u_sol', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845', 'vtkValidPointMask', '0.220004577706569', '0.4899977111467155', '0.7199969481956207', 'Points_X', '0.30000762951094834', '0.6899977111467155', '0.2899977111467155', 'Points_Y', '0.6', '0.3100022888532845', '0.6399938963912413', 'Points_Z', '1', '0.5000076295109483', '0', 'Points_Magnitude', '0.6500038147554742', '0.3400015259021897', '0.16000610360875867']
			plotOverLine1Display_1.SeriesPlotCorner = ['arc_length', '0', 'u_sol', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
			plotOverLine1Display_1.SeriesLabelPrefix = ''
			plotOverLine1Display_1.SeriesLineStyle = ['arc_length', '1', 'u_sol', '1', 'vtkValidPointMask', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'Points_Magnitude', '1']
			plotOverLine1Display_1.SeriesLineThickness = ['arc_length', '2', 'u_sol', '2', 'vtkValidPointMask', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'Points_Magnitude', '2']
			plotOverLine1Display_1.SeriesMarkerStyle = ['arc_length', '0', 'u_sol', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
			plotOverLine1Display_1.SeriesMarkerSize = ['arc_length', '4', 'u_sol', '4', 'vtkValidPointMask', '4', 'Points_X', '4', 'Points_Y', '4', 'Points_Z', '4', 'Points_Magnitude', '4']

			# get layout
			layout1 = GetLayoutByName("Layout #1")

			# add view to a layout so it's visible in UI
			AssignViewToLayout(view=lineChartView1, layout=layout1, hint=0)

			# Properties modified on plotOverLine1Display_1
			plotOverLine1Display_1.SeriesPlotCorner = ['Points_Magnitude', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'arc_length', '0', 'u_sol', '0', 'vtkValidPointMask', '0']
			plotOverLine1Display_1.SeriesLineStyle = ['Points_Magnitude', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'arc_length', '1', 'u_sol', '1', 'vtkValidPointMask', '1']
			plotOverLine1Display_1.SeriesLineThickness = ['Points_Magnitude', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'arc_length', '2', 'u_sol', '2', 'vtkValidPointMask', '2']
			plotOverLine1Display_1.SeriesMarkerStyle = ['Points_Magnitude', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'arc_length', '0', 'u_sol', '0', 'vtkValidPointMask', '0']
			plotOverLine1Display_1.SeriesMarkerSize = ['Points_Magnitude', '4', 'Points_X', '4', 'Points_Y', '4', 'Points_Z', '4', 'arc_length', '4', 'u_sol', '4', 'vtkValidPointMask', '4']

			# set active source
			SetActiveSource(upvd)

			# toggle 3D widget visibility (only when running from the GUI)
			Hide3DWidgets(proxy=plotOverLine1.Source)

			# set active source
			SetActiveSource(plotOverLine1)

			# save data
			SaveData('/Users/gxt/Desktop/My_Computer/open-sources/SS-Div-Conforming-B-splines-Incompressible-NS/Lid2D/Lid-2d-gamma0.01/SS_results/ss_Re='+str(Re)+'_p='+str(p)+'_m='+str(m)+'/u_y.csv', proxy=plotOverLine1, PointDataArrays=['arc_length', 'u_sol', 'vtkValidPointMask'])

			ResetSession()

			# create a new 'PVD Reader'
			vpvd = PVDReader(registrationName='v.pvd', FileName='/Users/gxt/Desktop/My_Computer/open-sources/SS-Div-Conforming-B-splines-Incompressible-NS/Lid2D/Lid-2d-gamma0.01/SS_results/ss_Re='+str(Re)+'_p='+str(p)+'_m='+str(m)+'/v.pvd')
			vpvd.PointArrays = ['v_sol']

			# get active view
			renderView1_1 = GetActiveViewOrCreate('RenderView')

			# show data in view
			vpvdDisplay = Show(vpvd, renderView1_1, 'UnstructuredGridRepresentation')

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

			# reset view to fit data
			renderView1_1.ResetCamera()

			#changing interaction mode based on data extents
			renderView1_1.InteractionMode = '2D'
			renderView1_1.CameraPosition = [0.5, 0.5, 10000.0]
			renderView1_1.CameraFocalPoint = [0.5, 0.5, 0.0]

			# get the material library
			materialLibrary1_1 = GetMaterialLibrary()

			# show color bar/color legend
			vpvdDisplay.SetScalarBarVisibility(renderView1_1, True)

			# update the view to ensure updated data information
			renderView1_1.Update()

			# create a new 'Plot Over Line'
			plotOverLine1_1 = PlotOverLine(registrationName='PlotOverLine1', Input=vpvd,
			    Source='Line')

			# init the 'Line' selected for 'Source'
			plotOverLine1_1.Source.Point1 = [0.0, 0.5, 0.0]
			plotOverLine1_1.Source.Point2 = [1.0, 0.5, 0.0]

			# show data in view
			plotOverLine1_1Display = Show(plotOverLine1_1, renderView1_1, 'GeometryRepresentation')

			# trace defaults for the display properties.
			plotOverLine1_1Display.Representation = 'Surface'
			plotOverLine1_1Display.ColorArrayName = ['POINTS', 'v_sol']
			plotOverLine1_1Display.LookupTable = v_solLUT
			plotOverLine1_1Display.SelectTCoordArray = 'None'
			plotOverLine1_1Display.SelectNormalArray = 'None'
			plotOverLine1_1Display.SelectTangentArray = 'None'
			plotOverLine1_1Display.OSPRayScaleArray = 'v_sol'
			plotOverLine1_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
			plotOverLine1_1Display.SelectOrientationVectors = 'None'
			plotOverLine1_1Display.ScaleFactor = 0.1
			plotOverLine1_1Display.SelectScaleArray = 'v_sol'
			plotOverLine1_1Display.GlyphType = 'Arrow'
			plotOverLine1_1Display.GlyphTableIndexArray = 'v_sol'
			plotOverLine1_1Display.GaussianRadius = 0.005
			plotOverLine1_1Display.SetScaleArray = ['POINTS', 'v_sol']
			plotOverLine1_1Display.ScaleTransferFunction = 'PiecewiseFunction'
			plotOverLine1_1Display.OpacityArray = ['POINTS', 'v_sol']
			plotOverLine1_1Display.OpacityTransferFunction = 'PiecewiseFunction'
			plotOverLine1_1Display.DataAxesGrid = 'GridAxesRepresentation'
			plotOverLine1_1Display.PolarAxes = 'PolarAxesRepresentation'

			# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
			plotOverLine1_1Display.ScaleTransferFunction.Points = [-0.4629252605300211, 0.0, 0.5, 0.0, 0.37680873240806095, 1.0, 0.5, 0.0]

			# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
			plotOverLine1_1Display.OpacityTransferFunction.Points = [-0.4629252605300211, 0.0, 0.5, 0.0, 0.37680873240806095, 1.0, 0.5, 0.0]

			# Create a new 'Line Chart View'
			lineChartView1_1 = CreateView('XYChartView')

			# show data in view
			plotOverLine1_1Display_1 = Show(plotOverLine1_1, lineChartView1_1, 'XYChartRepresentation')

			# trace defaults for the display properties.
			plotOverLine1_1Display_1.CompositeDataSetIndex = [0]
			plotOverLine1_1Display_1.UseIndexForXAxis = 0
			plotOverLine1_1Display_1.XArrayName = 'arc_length'
			plotOverLine1_1Display_1.SeriesVisibility = ['v_sol']
			plotOverLine1_1Display_1.SeriesLabel = ['arc_length', 'arc_length', 'v_sol', 'v_sol', 'vtkValidPointMask', 'vtkValidPointMask', 'Points_X', 'Points_X', 'Points_Y', 'Points_Y', 'Points_Z', 'Points_Z', 'Points_Magnitude', 'Points_Magnitude']
			plotOverLine1_1Display_1.SeriesColor = ['arc_length', '0', '0', '0', 'v_sol', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845', 'vtkValidPointMask', '0.220004577706569', '0.4899977111467155', '0.7199969481956207', 'Points_X', '0.30000762951094834', '0.6899977111467155', '0.2899977111467155', 'Points_Y', '0.6', '0.3100022888532845', '0.6399938963912413', 'Points_Z', '1', '0.5000076295109483', '0', 'Points_Magnitude', '0.6500038147554742', '0.3400015259021897', '0.16000610360875867']
			plotOverLine1_1Display_1.SeriesPlotCorner = ['arc_length', '0', 'v_sol', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
			plotOverLine1_1Display_1.SeriesLabelPrefix = ''
			plotOverLine1_1Display_1.SeriesLineStyle = ['arc_length', '1', 'v_sol', '1', 'vtkValidPointMask', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'Points_Magnitude', '1']
			plotOverLine1_1Display_1.SeriesLineThickness = ['arc_length', '2', 'v_sol', '2', 'vtkValidPointMask', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'Points_Magnitude', '2']
			plotOverLine1_1Display_1.SeriesMarkerStyle = ['arc_length', '0', 'v_sol', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
			plotOverLine1_1Display_1.SeriesMarkerSize = ['arc_length', '4', 'v_sol', '4', 'vtkValidPointMask', '4', 'Points_X', '4', 'Points_Y', '4', 'Points_Z', '4', 'Points_Magnitude', '4']

			# get layout
			layout1_1 = GetLayoutByName("Layout #1")

			# add view to a layout so it's visible in UI
			AssignViewToLayout(view=lineChartView1_1, layout=layout1_1, hint=0)

			# Properties modified on plotOverLine1_1Display_1
			plotOverLine1_1Display_1.SeriesPlotCorner = ['Points_Magnitude', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'arc_length', '0', 'v_sol', '0', 'vtkValidPointMask', '0']
			plotOverLine1_1Display_1.SeriesLineStyle = ['Points_Magnitude', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'arc_length', '1', 'v_sol', '1', 'vtkValidPointMask', '1']
			plotOverLine1_1Display_1.SeriesLineThickness = ['Points_Magnitude', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'arc_length', '2', 'v_sol', '2', 'vtkValidPointMask', '2']
			plotOverLine1_1Display_1.SeriesMarkerStyle = ['Points_Magnitude', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'arc_length', '0', 'v_sol', '0', 'vtkValidPointMask', '0']
			plotOverLine1_1Display_1.SeriesMarkerSize = ['Points_Magnitude', '4', 'Points_X', '4', 'Points_Y', '4', 'Points_Z', '4', 'arc_length', '4', 'v_sol', '4', 'vtkValidPointMask', '4']

			# save data
			SaveData('/Users/gxt/Desktop/My_Computer/open-sources/SS-Div-Conforming-B-splines-Incompressible-NS/Lid2D/Lid-2d-gamma0.01/SS_results/ss_Re='+str(Re)+'_p='+str(p)+'_m='+str(m)+'/v_x.csv', proxy=plotOverLine1_1, PointDataArrays=['arc_length', 'v_sol', 'vtkValidPointMask'])

			ResetSession()

			# get active view
			renderView1_2 = GetActiveViewOrCreate('RenderView')

			#================================================================
			# addendum: following script captures some of the application
			# state to faithfully reproduce the visualization during playback
			#================================================================

			# get layout
			layout1_2 = GetLayout()

			#--------------------------------
			# saving layout sizes for layouts

			# layout/tab size in pixels
			layout1_2.SetSize(2594, 1608)

			#--------------------------------------------
			# uncomment the following to render all views
			# RenderAllViews()
			# alternatively, if you want to write images, you can use SaveScreenshot(...).