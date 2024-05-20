import os
import sys
sys.path.append('/home/felipe/pythonFunc/mainfuncsPy')
import vtk
from vtk.numpy_interface import algorithms as algs
from vtk.numpy_interface import dataset_adapter as dsa
import numpy as np
import scipy
import h5py
import multiprocessing 



global POSTPROCESS

InitialFlow = int(sys.argv[1])
EndFlow = int(sys.argv[2])
POSTPROCESS = int(sys.argv[3])
numprocs = int(sys.argv[4])

if POSTPROCESS == 1:
    POSTPROCESS = True
if POSTPROCESS == 0:
    POSTPROCESS = False



files_path = '/home/felipe/autolst/m05-FINAL/'    

clear = lambda: os.system('clear')  # On Windows System
clear()

   

def vtkconvert(data,mesh,outputfilepath,outputfilepath_postporcessing):
    #` Load mesh data
    print('load mesh')
    x = np.array(mesh['X'][0],dtype=np.float64)
    y = np.array(mesh['Y'][0],dtype=np.float64)
    z = np.array(mesh['Z'][0],dtype=np.float64)
    #x1 = np.abs(mesh['x1'][0])
    #x2 = np.abs(mesh['x2'][0])
    #y1 = np.abs(mesh['y1'][0])
    #y2 = np.abs(mesh['y2'][0])
    #z1 = np.abs(mesh['z1'][0])
    #z2 = np.abs(mesh['z2'][0])
    #L = np.abs(mesh['L'][0])
    #D = np.abs(mesh['D'][0])
    x1 = 250
    x2 = 500
    y1 = -6
    y2 = 20
    z1 = -20
    z2 = 20
    L = 12.22
    D = 6.11


    # Set display interval
    x_start = np.abs(x - (x1-2*L)).argmin()
    x_end = np.abs(x - (x2)).argmin()
    y_start = np.abs(y -y1).argmin()
    y_end = np.abs(y - y2).argmin()
    z_start = np.abs(z - z1).argmin()
    z_end = np.abs(z - z2).argmin()

    x_factor = 1
    y_factor = 1
    z_factor = 1

    x = x[x_start:x_end][::x_factor]
    y = y[y_start:y_end][::y_factor]
    z = z[z_start:z_end][::z_factor]

    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    wall = np.array(mesh['wall'])[x_start:x_end,y_start:y_end,z_start:z_end][::x_factor,::y_factor,::z_factor]

    print('load data')

    U = np.array(data['U']).transpose((2,1,0))[x_start:x_end,y_start:y_end,z_start:z_end][::x_factor,::y_factor,::z_factor]
    V = np.array(data['V']).transpose((2,1,0))[x_start:x_end,y_start:y_end,z_start:z_end][::x_factor,::y_factor,::z_factor]
    W = np.array(data['W']).transpose((2,1,0))[x_start:x_end,y_start:y_end,z_start:z_end][::x_factor,::y_factor,::z_factor]
    R = np.array(data['E']).transpose((2,1,0))[x_start:x_end,y_start:y_end,z_start:z_end][::x_factor,::y_factor,::z_factor]
    E = np.array(data['R']).transpose((2,1,0))[x_start:x_end,y_start:y_end,z_start:z_end][::x_factor,::y_factor,::z_factor]    

    U[np.isnan(U)] = 0
    V[np.isnan(V)] = 0
    W[np.isnan(W)] = 0
        
    if POSTPROCESS == True:

        L2 = np.array(data['L2']).transpose((2,1,0))[x_start:x_end,y_start:y_end,z_start:z_end][::x_factor,::y_factor,::z_factor]
        Q  = np.array(data['Q']).transpose((2,1,0))[x_start:x_end,y_start:y_end,z_start:z_end][::x_factor,::y_factor,::z_factor]
        #Mach    =      calculate_mach_number(U, E, R)
        #Vortk   =      calculate_vorticity_z(U, V, dx, dy)
        #Dilatation = np.array(data['Dilatation']).transpose((2,1,0))[x_start:x_end,y_start:y_end,z_start:z_end][::x_factor,::y_factor,::z_factor]

        Q[np.isnan(Q)] = 0
        L2[np.isnan(L2)] = 0
        #Mach[np.isnan(Mach)] = 0
        #Vortk[np.isnan(Vortk)] = 0
       # Dilatation[np.isnan(Dilatation)] = 0


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
    print('configure data')
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

    print('write data on vtk')
    # Create Writer object to save the VTK file
    writer = vtk.vtkRectilinearGridWriter()
    writer.SetFileName(outputfilepath)
    writer.SetInputData(my_vtk_dataset)
    writer.Update()
    writer.Write()
    
    if POSTPROCESS == True:
    
        print('write postprocessing data on vtk')
        writer = vtk.vtkRectilinearGridWriter()
        writer.SetFileName(outputfilepath_postporcessing)
        writer.SetInputData(my_vtk_dataset_postporcessing)
        writer.Update()
        writer.Write()


mesh = scipy.io.loadmat(r"/home/felipe/autolst/m05-FINAL/mesh.mat")


# Nome base dos arquivos
base_nome_arquivo = "flow_{:010d}.mat"
base_nome_arquivovtk = "flow_{:010d}.vtk"
base_nome_arquivovtk_postporcessing= "postprocessing_{:010d}.vtk"
numeros_arquivos = range(InitialFlow,EndFlow,1)

for numero_arquivo in numeros_arquivos:
    clear()
    print('flow:'+str(numero_arquivo))
    nome_arquivo = base_nome_arquivo.format(numero_arquivo)
    nome_arquivovtk = base_nome_arquivovtk.format(numero_arquivo)
    nome_arquivovtk_postporcessing = base_nome_arquivovtk_postporcessing.format(numero_arquivo)
    caminho_arquivo = os.path.join(r"/home/felipe/autolst/m05-FINAL/", nome_arquivo)
    caminho_arquivovtk = os.path.join(r"/home/felipe/OLIVER-POSTPROCESS/caseD/lambda2/", nome_arquivovtk)
    caminho_arquivo_postporcessing_vtk = os.path.join(r"/home/felipe/OLIVER-POSTPROCESS/caseD/lambda2/", nome_arquivovtk_postporcessing)
    if os.path.exists(caminho_arquivo):
        print(f"O arquivo {nome_arquivo} existe. Carregando e chamando a função convertvtk.")
        data = h5py.File(caminho_arquivo,'r')
        vtkconvert(data,mesh,caminho_arquivovtk,caminho_arquivo_postporcessing_vtk)

       
    else:
        print(f" Não entrou no looping anterior / Arquivo não encontrado ")
   
        continue