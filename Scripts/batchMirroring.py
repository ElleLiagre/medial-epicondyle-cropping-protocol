import slicer
import qt
import os
import vtk
import logging
from typing import List, Optional

# Configure logging
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

# Remove existing handlers if they exist to avoid duplicates
for handler in logger.handlers[:]:
    logger.removeHandler(handler)

# Define a console handler for real-time logging output
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


def load_model_from_file(inputFilePath: str) -> Optional[vtk.vtkPolyData]:
    modelNode = slicer.util.loadModel(inputFilePath)
    if not modelNode:
        logging.error(f"Failed to load model: {inputFilePath}")
    return modelNode

def apply_mirror_transform(modelNode: vtk.vtkPolyData) -> None:
    if modelNode:
        mirrorTransform = vtk.vtkTransform()
        mirrorTransform.Scale(-1, 1, 1)  # Mirror along the x-axis
        modelNode.ApplyTransform(mirrorTransform)
        logging.info(f'Applied mirror transformation to model: {modelNode.GetName()}.')

def recalculate_normals(modelNode: vtk.vtkPolyData) -> None:
    if modelNode:
        normals = vtk.vtkPolyDataNormals()
        normals.SetInputData(modelNode.GetPolyData())
        normals.SetFlipNormals(True)
        normals.Update()
        modelNode.SetAndObservePolyData(normals.GetOutput())
        logging.info(f'Recalculated normals for model: {modelNode.GetName()}.')

def save_model(modelNode: vtk.vtkPolyData, outputModelName: str) -> None:
    if modelNode:
        slicer.util.saveNode(modelNode, outputModelName)
        logging.info(f'Successfully exported mirrored model as: {outputModelName}')

def is_valid_folder(folderPath: str) -> bool:
    if not os.path.isdir(folderPath):
        logging.error(f"The folder path '{folderPath}' does not exist.")
        return False
    if not os.listdir(folderPath):
        logging.error(f"The folder '{folderPath}' is empty.")
        return False
    return True

def handle_filename_conflict(outputModelName: str) -> str:
    base, ext = os.path.splitext(outputModelName)
    counter = 1
    while os.path.exists(outputModelName):
        outputModelName = f"{base}_{counter}{ext}"
        counter += 1
    return outputModelName

def batchMirroring(folderPath: Optional[str] = None, 
                   valid_extensions: List[str] = ['.obj', '.ply', '.stl'], 
                   output_format: str = 'ply') -> None:
    """
    Batch process 3D models in a selected folder to apply a mirror transformation.
    The mirrored models will be saved in a new folder named 'OriginalFolderName_Mirrored'.

    Parameters:
        folderPath (str): Path to the folder containing models. If None, the user will be prompted.
        valid_extensions (list): List of valid file extensions to process.
        output_format (str): The format to save mirrored models (e.g., 'ply', 'stl').
    """

    # Prompt the user to select a folder if not provided
    if not folderPath:
        folderDialog = qt.QFileDialog()
        folderDialog.setFileMode(qt.QFileDialog.Directory)
        folderDialog.setOption(qt.QFileDialog.ShowDirsOnly, True)
        if folderDialog.exec_() == qt.QFileDialog.Accepted:
            folderPath = folderDialog.selectedFiles()[0]
            logging.info(f"Selected folder: {folderPath}")
        else:
            logging.warning("Folder selection canceled.")
            return

    # Validate folder and files
    if not is_valid_folder(folderPath):
        return

    # Create a new folder named "OriginalFolderName_Mirrored" inside the selected folder
    originalFolderName = os.path.basename(folderPath)
    mirroredFolder = os.path.join(folderPath, f'{originalFolderName}_Mirrored')
    os.makedirs(mirroredFolder, exist_ok=True)
    logging.info(f"Created mirrored folder: {mirroredFolder}")

    # Add file logging in the new folder
    log_file_path = os.path.join(mirroredFolder, "batch_mirroring.log")
    file_handler = logging.FileHandler(log_file_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # Loop over and process each valid file in the selected folder
    for inputFilename in os.listdir(folderPath):
        fn, ext = os.path.splitext(inputFilename)
        if ext.lower() not in valid_extensions:
            logging.warning(f"Skipping unsupported file: {inputFilename}")
            continue

        # Load model with exception handling
        inputFilePath = os.path.join(folderPath, inputFilename)
        modelNode = load_model_from_file(inputFilePath)
        
        if modelNode is None:
            continue

        try:
            apply_mirror_transform(modelNode)
            recalculate_normals(modelNode)

            # Save mirrored model with conflict handling
            outputModelName = os.path.join(mirroredFolder, f"{fn}_mirrored.{output_format}")
            outputModelName = handle_filename_conflict(outputModelName)
            save_model(modelNode, outputModelName)

        except vtk.vtkError as e:
            logging.error(f"VTK error processing file '{inputFilename}': {e}")
        except Exception as e:
            logging.error(f"Error processing file '{inputFilename}': {e}")
        finally:
            if modelNode is not None:
                slicer.mrmlScene.RemoveNode(modelNode)

    logging.info(f'All models mirrored and saved in the "{originalFolderName}_Mirrored" folder.')

