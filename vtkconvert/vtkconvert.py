import os
import sys
import vtk
from vtk.numpy_interface import algorithms as algs
from vtk.numpy_interface import dataset_adapter as dsa
import numpy as np
import scipy
import h5py
#import multiprocessing 


#input vars

global POSTPROCESS

InitialFlow = int(sys.argv[1])
EndFlow = int(sys.argv[2])
POSTPROCESS = int(sys.argv[3])
numprocs = int(sys.argv[4])

if POSTPROCESS == 1:
    POSTPROCESS = True
if POSTPROCESS == 0:
    POSTPROCESS = False



files_path = ''    
output_path = ''

#clear = lambda: os.system('clear')  # On Windows System
#clear()

   

def vtkconvert(data,mesh,outputfilepath,outputfilepath_postporcessing):
    #` Load mesh data
    print('Convert mesh')
    x = np.array(mesh['X'][0],dtype=np.float64)
    y = np.array(mesh['Y'][0],dtype=np.float64)
    z = np.array(mesh['Z'][0],dtype=np.float64)    

    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    wall = np.array(mesh['wall'])

    print('Convert data')

    U = np.array(data['U']).transpose((2,1,0))
    V = np.array(data['V']).transpose((2,1,0))
    W = np.array(data['W']).transpose((2,1,0))
    R = np.array(data['E']).transpose((2,1,0))
    E = np.array(data['R']).transpose((2,1,0)) 

    U[np.isnan(U)] = 0
    V[np.isnan(V)] = 0
    W[np.isnan(W)] = 0
        
    if POSTPROCESS == True:

        L2 = np.array(data['L2']).transpose((2,1,0))
        Q  = np.array(data['Q']).transpose((2,1,0))
        #Mach    =  np.array(data['Ma']).transpose((2,1,0))
        #Vortk   =  np.array(data['Vortk']).transpose((2,1,0))
        #Dilatation = np.array(data['Dilatation']).transpose((2,1,0))

        Q[np.isnan(Q)] = 0
        L2[np.isnan(L2)] = 0
        #Mach[np.isnan(Mach)] = 0
        #Vortk[np.isnan(Vortk)] = 0
        #Dilatation[np.isnan(Dilatation)] = 0


    ## PART 1

    # Create object vtkRectiliniarGrid
    my_vtk_dataset = vtk.vtkRectilinearGrid()
    my_vtk_dataset_postporcessing = vtk.vtkRectilinearGrid()

    ## CONFIGURE POINTS

    # Set grid dimensions (IMPORTANT)
    my_vtk_dataset.SetDimensions(len(x),len(y),len(z))
    my_vtk_dataset_postporcessing.SetDimensions(len(x),len(y),len(z))

    # Transform coordinate vectors into VTK format
    x_vtk = dsa.numpyTovtkDataArray(x)
    y_vtk = dsa.numpyTovtkDataArray(y)
    z_vtk = dsa.numpyTovtkDataArray(z)

    # Set grid coordinates using VTK vectors
    my_vtk_dataset.SetXCoordinates(x_vtk)
    my_vtk_dataset.SetYCoordinates(y_vtk)
    my_vtk_dataset.SetZCoordinates(z_vtk)
    my_vtk_dataset_postporcessing.SetXCoordinates(x_vtk)
    my_vtk_dataset_postporcessing.SetYCoordinates(y_vtk)
    my_vtk_dataset_postporcessing.SetZCoordinates(z_vtk)
    # IT IS NOT NECESSARY TO CREATE A POINTS OBJECT FOR RECTILINIARGRID OBJECT
    #pts = vtk.vtkPoints()
    #points = algs.make_vector(X.flatten('F'),
    #                          Y.flatten('F'),
    #                          Z.flatten('F'))
    #pts.SetData(dsa.numpyTovtkDataArray(points, "Points"))
    #my_vtk_dataset.GetPoints(pts)


    ## CONFIGURE DATA
    print('Configure data')
    # Create vector matrix on VTK format
    if POSTPROCESS== True:
        
        my_vtk_dataset_postporcessing.GetPointData().AddArray(dsa.numpyTovtkDataArray(wall.flatten('F'), "wall"))
        my_vtk_dataset_postporcessing.GetPointData().AddArray(dsa.numpyTovtkDataArray(L2.flatten('F'), "L2"))
        my_vtk_dataset_postporcessing.GetPointData().AddArray(dsa.numpyTovtkDataArray(Q.flatten('F'), "Q"))
        #my_vtk_dataset.GetPointData().AddArray(dsa.numpyTovtkDataArray(Dilatation.flatten('F'), "Dilatation"))
        #my_vtk_dataset_postporcessing.GetPointData().AddArray(dsa.numpyTovtkDataArray(Vortk.flatten('F'), "Vortk"))
        #my_vtk_dataset_postporcessing.GetPointData().AddArray(dsa.numpyTovtkDataArray(Mach.flatten('F'), "Mach"))
    
    vectors = algs.make_vector(U.flatten('F'),
                            V.flatten('F'),
                            W.flatten('F'))

    # Use AddArray to add the data
    my_vtk_dataset.GetPointData().AddArray(dsa.numpyTovtkDataArray(wall.flatten('F'), "wall"))
    my_vtk_dataset.GetPointData().AddArray(dsa.numpyTovtkDataArray(vectors, "Velocity"))
    my_vtk_dataset.GetPointData().AddArray(dsa.numpyTovtkDataArray(U.flatten('F'), "U"))
    my_vtk_dataset.GetPointData().AddArray(dsa.numpyTovtkDataArray(V.flatten('F'), "V"))
    my_vtk_dataset.GetPointData().AddArray(dsa.numpyTovtkDataArray(W.flatten('F'), "W"))

    print('Write data on vtk')
    # Create Writer object to save the VTK file
    writer = vtk.vtkRectilinearGridWriter()
    writer.SetFileName(outputfilepath)
    writer.SetInputData(my_vtk_dataset)
    writer.Update()
    writer.Write()
    
    if POSTPROCESS == True:
    
        print('Write postprocessing data on vtk')
        writer = vtk.vtkRectilinearGridWriter()
        writer.SetFileName(outputfilepath_postporcessing)
        writer.SetInputData(my_vtk_dataset_postporcessing)
        writer.Update()
        writer.Write()


mesh = scipy.io.loadmat(r,files_path,"mesh.mat")


# file names
file_base_name = "flow_{:010d}.mat"
file_base_name_vtk = "flow_{:010d}.vtk"
file_base_name_vtk_postporcessing= "postprocessing_{:010d}.vtk"
files_number = range(InitialFlow,EndFlow,1)

for file_number in files_number:
    #clear()
    print('flow:'+str(file_number))
    file_name = file_base_name.format(file_number)
    file_name_vtk = file_base_name_vtk.format(file_number)
    file_name_postprocess_vtk = file_base_name_vtk_postporcessing.format(file_number)
    file_path = os.path.join(files_path, file_name)
    file_path_vtk = os.path.join(output_path, file_name_vtk)
    file_path_postporcessing_vtk = os.path.join(output_path, file_name_postprocess_vtk)
    if os.path.exists(file_path):
        print(f" File: {file_name}. Loading and converting.")
        data = h5py.File(file_path,'r')
        vtkconvert(data,mesh,file_path_vtk,file_path_postporcessing_vtk)

       
    else:
        print(f" Not found or skipped the looping ")
   
        continue