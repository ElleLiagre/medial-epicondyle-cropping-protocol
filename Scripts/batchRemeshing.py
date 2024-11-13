import slicer
import qt
import os
import vtk
import logging
from typing import List, Optional
import pymeshlab

# Configure logging
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

for handler in logger.handlers[:]:
    logger.removeHandler(handler)

# Define a formatter for consistency
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

# Create a console handler for real-time logging output
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.DEBUG)
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

# Constants for remeshing parameters
TARGET_LENGTH = 0.3
MAX_SURFACE_DISTANCE = 0.3

def ensure_pymeshlab_installed():
    try:
        import pymeshlab
        logger.info("pymeshlab is already installed.")
    except ImportError:
        logger.info("Installing pymeshlab...")
        slicer.util.pip_install("pymeshlab")
        logger.info("pymeshlab installation complete.")

def process_mesh_file(mesh_file_path: str, res: int) -> pymeshlab.MeshSet:
    ms = pymeshlab.MeshSet()
    ms.load_new_mesh(mesh_file_path)

    # Step 1: Apply isotropic remeshing
    ms.meshing_isotropic_explicit_remeshing(
        targetlen=pymeshlab.PureValue(TARGET_LENGTH), 
        maxsurfdist=pymeshlab.PureValue(MAX_SURFACE_DISTANCE)
    )

    # Step 2: Get geometric measures
    geometric_measures = ms.get_geometric_measures()
    surface_area = geometric_measures['surface_area']
    current_faces = ms.current_mesh().face_number()

    # Step 3: Calculate the target number of faces
    target_faces = int(surface_area * res)
    logger.info(f"Current number of faces: {current_faces}, Target number of faces: {target_faces}")

    # Step 4: Apply Quadric Edge Collapse Decimation
    ms.meshing_decimation_quadric_edge_collapse(
        targetfacenum=target_faces, 
        preservenormal=True, 
        preservetopology=True, 
        preserveboundary=True, 
        optimalplacement=True, 
        qualitythr=1.0
    )

    logger.info(f"Mesh processing complete for {os.path.basename(mesh_file_path)}")
    return ms

def load_model_from_file(inputFilePath: str) -> Optional[vtk.vtkPolyData]:
    modelNode = slicer.util.loadModel(inputFilePath)
    if not modelNode:
        logger.error(f"Failed to load model: {inputFilePath}")
    return modelNode

def save_processed_mesh(ms: pymeshlab.MeshSet, outputModelName: str):
    ms.save_current_mesh(outputModelName)
    logger.info(f"Saved remeshed model as: {outputModelName}")

def is_valid_folder(folderPath: str) -> bool:
    if not os.path.isdir(folderPath):
        logger.error(f"The folder path '{folderPath}' does not exist.")
        return False
    if not os.listdir(folderPath):
        logger.error(f"The folder '{folderPath}' is empty.")
        return False
    return True

def handle_filename_conflict(outputModelName: str) -> str:
    base, ext = os.path.splitext(outputModelName)
    counter = 1
    while os.path.exists(outputModelName):
        outputModelName = f"{base}_{counter}{ext}"
        counter += 1
    return outputModelName 

def batchRemeshing(folderPath: Optional[str] = None, 
                   res: int = 22, 
                   valid_extensions: List[str] = ['.obj', '.ply', '.stl'], 
                   output_format: str = 'ply') -> None:
    if not folderPath:
        folderDialog = qt.QFileDialog()
        folderDialog.setFileMode(qt.QFileDialog.Directory)
        folderDialog.setOption(qt.QFileDialog.ShowDirsOnly, True)
        if folderDialog.exec_() == qt.QFileDialog.Accepted:
            folderPath = folderDialog.selectedFiles()[0]
            logger.info(f"Selected folder: {folderPath}")
        else:
            logger.warning("Folder selection canceled.")
            return

    if not is_valid_folder(folderPath):
        return

    originalFolderName = os.path.basename(folderPath)
    remeshedFolder = os.path.join(folderPath, f'{originalFolderName}_Remeshed')
    os.makedirs(remeshedFolder, exist_ok=True)
    logger.info(f"Created remeshed folder: {remeshedFolder}")

    # Add file logging in the new folder
    log_file_path = os.path.join(remeshedFolder, "batch_remeshing.log")
    file_handler = logging.FileHandler(log_file_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # Loop over and process each valid file in the selected folder
    for inputFilename in os.listdir(folderPath):
        fn, ext = os.path.splitext(inputFilename)
        if ext.lower() not in valid_extensions:
            logger.warning(f"Skipping unsupported file: {inputFilename}")
            continue

        inputFilePath = os.path.join(folderPath, inputFilename)
        modelNode = slicer.util.loadModel(inputFilePath)
        
        if modelNode is None:
            continue

        try:
            ms = process_mesh_file(inputFilePath, res)
            outputModelName = os.path.join(remeshedFolder, f"{fn}_remeshed.{output_format}")
            outputModelName = handle_filename_conflict(outputModelName)
            save_processed_mesh(ms, outputModelName)

        except Exception as e:
            logger.error(f"Error processing file '{inputFilename}': {e}")
        finally:
            if modelNode is not None:
                slicer.mrmlScene.RemoveNode(modelNode)

    logger.info(f'All models remeshed and saved in the "{originalFolderName}_Remeshed" folder.')

