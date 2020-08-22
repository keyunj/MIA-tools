import vtk
import os
from glob import glob
import msvcrt


class keyCallback(vtk.vtkCommand):
    def __init__(self, iren):
        vtk.vtkCommand.__init__()
        self.iren = iren

    def Execute():
        key = iren.GetKeySym()
        if key == 'w':
            # for each image
            img_list = glob('{}/*.mha'.format('.'))
            img_list.sort()
            
            for img_path in img_list:
                reader = vtk.vtkMetaImageReader()
                reader.SetFileName(img_path)
                reader.Update()

                volMapper.SetInputData(reader.GetOutput())

                renWin.Render()


reader = vtk.vtkMetaImageReader()
reader.SetFileName('./036.mha')
reader.Update()

volMapper = vtk.vtkSmartVolumeMapper()
volMapper.SetInputData(reader.GetOutput())
volMapper.SetBlendModeToComposite()

comp = vtk.vtkColorTransferFunction()
comp.AddRGBPoint(-1000.0, 0.0, 0.0, 0.0);
comp.AddRGBPoint(0.0, 0.5, 0.5, 0.5);
comp.AddRGBPoint(1000.0, 1.0, 1.0, 1.0);

volProp = vtk.vtkVolumeProperty()
volProp.ShadeOff()
volProp.SetInterpolationType(vtk.VTK_LINEAR_INTERPOLATION)
volProp.SetColor(comp)

vol = vtk.vtkVolume()
vol.SetMapper(volMapper)
vol.SetProperty(volProp)

ren = vtk.vtkRenderer()
ren.AddViewProp(vol)
ren.SetBackground(0.5,0.5,0.5)

renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)

irenstyle = vtk.vtkInteractorStyleTrackballCamera()

iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)
iren.SetInteractorStyle(irenstyle)

keycb = keyCallback(iren)
iren.AddObserver(vtk.vtkCommand.CharEvent, keycb)

iren.Initialize()
renWin.Render()
iren.Start()



