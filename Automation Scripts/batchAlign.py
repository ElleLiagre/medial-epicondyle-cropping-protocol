# ===========================================================================
# Title: A Standardized, Three-Dimensional Cropping Protocol for Analyzing 
#        the Medial Epicondyle of the Humerus - Script for Batch Alignment
# 
# Author: Elle B. K. Liagre
# Email: elle.liagre@u-bordeaux.fr
# ORCID: https://orcid.org/0000-0002-8993-3266
# Version: 1.0
# Date: 2025-03-24
# 
# Description: Batch function to align all 3D models in a folder to a reference 
# model. Note: Requires the SlicerMorph extension in 3D Slicer. Load the 
# reference model in 3D Slicer and name it "target" before running the function.
# You can use the reference model provided, or your own model. The aligned 
# models will be saved in a new folder within the original folder. 
# More information can be found on: 
# https://github.com/ElleLiagre/medial-epicondyle-cropping-protocol
#
# Software Requirements: 
# - 3D Slicer (https://www.slicer.org/)
# - SlicerMorph extension (https://github.com/SlicerMorph)
#
# License: This script is licensed under the GNU General Public License v3.0.
# See https://www.gnu.org/licenses/gpl-3.0.html for more details.
# ===========================================================================  



import slicer
import qt
import os
import io
import re
import logging
from typing import List, Optional
import sys

# Parameters for alignment and thresholds
parameters_dict = {
    'pointDensity': 1.0,
    'normalSearchRadius': 2.0,
    'FPFHNeighbors': 100,
    'FPFHSearchRadius': 5.0,
    'distanceThreshold': 3.0,
    'maxRANSAC': 1000000,
    'ICPDistanceThreshold': 1.5
}

FITNESS_THRESHOLD = 0.95
RMSE_THRESHOLD = 3.4
MAX_RETRIES = 3

# Check for the SlicerMorph extension
extensionsModel = slicer.app.extensionsManagerModel()
managedExtensions = extensionsModel.managedExtensions

extensionName = 'SlicerMorph'
if extensionName in managedExtensions:
    print(f"The extension '{extensionName}' is installed.")
    import ALPACA
else:
    print(f"The extension '{extensionName}' is not installed.")
    sys.exit()

# Function definitions
def capture_registration_output(logic, sourceModel, targetModel, parameters_dict):
    old_stdout = sys.stdout  
    sys.stdout = io.StringIO()

    try:
        logic.ITKRegistration(sourceModel, targetModel, skipScalingOption=False, parameterDictionary=parameters_dict, usePoisson=False)
    except Exception as e:
        logging.error(f"Error during registration: {e}")
        return None, None
    finally:
        output = sys.stdout.getvalue()
        sys.stdout = old_stdout

    # Parse the fitness and RMSE from the captured output
    fitness_match = re.search(r'Fitness\s*=\s*([\d\.]+)', output)
    rmse_match = re.search(r'RMSE\s*=\s*([\d\.]+)', output)

    fitness = float(fitness_match.group(1)) if fitness_match else None
    rmse = float(rmse_match.group(1)) if rmse_match else None

    return fitness, rmse

def realign_model(sourceModel, targetModel, parameters_dict):
    logic = slicer.util.getModuleLogic('FastModelAlign')
    return capture_registration_output(logic, sourceModel, targetModel, parameters_dict)

def is_successful_alignment(fitness, rmse):
    return (fitness is not None and fitness >= FITNESS_THRESHOLD) and (rmse is not None and rmse <= RMSE_THRESHOLD)

def process_model(sourceModel, targetModel, parameters_dict):
    results = []
    fitness, rmse = realign_model(sourceModel, targetModel, parameters_dict)
    results.append({'attempt': 'Initial', 'fitness': fitness, 'rmse': rmse})

    if is_successful_alignment(fitness, rmse):
        logging.info("Alignment was successful on the initial attempt.")
    else:
        logging.info("Initial alignment did not meet quality thresholds, starting retry attempts.")
        for attempt in range(1, MAX_RETRIES + 1):
            parameters_dict['distanceThreshold'] *= 2
            fitness, rmse = realign_model(sourceModel, targetModel, parameters_dict)
            results.append({'attempt': f'Retry {attempt}', 'fitness': fitness, 'rmse': rmse})

            if is_successful_alignment(fitness, rmse):
                logging.info(f"Alignment was successful on attempt {attempt}.")
                break
            else:
                logging.info(f"Retry {attempt} did not meet quality thresholds: Fitness = {fitness}, RMSE = {rmse}")

    parameters_dict['distanceThreshold'] /= (2 ** min(MAX_RETRIES, len(results) - 1))
    return results

def save_aligned_model(sourceModel, orientedFolder, inputFilename):
    output_name = os.path.splitext(inputFilename)[0] + '_oriented.ply'
    sourceOutput = os.path.join(orientedFolder, output_name)

    try:
        slicer.util.saveNode(sourceModel, sourceOutput)
        logging.info(f'Successfully saved aligned model: {sourceOutput}')
    except Exception as e:
        logging.error(f"Failed to save model: {e}")

def is_valid_folder(folderPath: str) -> bool:
    if not os.path.isdir(folderPath):
        logging.error(f"The folder path '{folderPath}' does not exist.")
        return False
    if not os.listdir(folderPath):
        logging.error(f"The folder '{folderPath}' is empty.")
        return False
    return True

def batchAlign(folderPath: Optional[str] = None) -> None:
    """
    Batch aligns 3D models in a specified folder to a target model.
    Aligned models will be saved in a new folder named 'OriginalFolderName_Oriented'.
    
    Parameters:
        folderPath (str): Path to the folder containing models. If None, the user will be prompted.
    """

    # Check if target model exists in the scene
    targetModel = slicer.util.getNode("target")
    if not targetModel:
        logging.error("Target model named 'target' does not exist in the scene. Please load or set the target model as 'target'.")
        print("Target model named 'target' does not exist in the scene. Please load or set the target model as 'target'.")
        return  

    # Prompt the user to select a folder if not provided
    if not folderPath:
        folderDialog = qt.QFileDialog()
        folderDialog.setFileMode(qt.QFileDialog.Directory)
        folderDialog.setOption(qt.QFileDialog.ShowDirsOnly, True)
        if folderDialog.exec_() != qt.QFileDialog.Accepted:
            logging.info("Folder selection cancelled.")
            return
        folderPath = folderDialog.selectedFiles()[0]

    # Validate folder and files
    if not is_valid_folder(folderPath):
        return

    # Create the "Oriented" output folder and the log file path within it
    originalFolderName = os.path.basename(folderPath)
    orientedFolder = os.path.join(folderPath, f'{originalFolderName}_Oriented')
    try:
        os.makedirs(orientedFolder, exist_ok=True)
    except Exception as e:
        logging.error(f"Could not create output folder: {e}")
        return
    
    # Configure logging to save in the orientedFolder
    log_file_path = os.path.join(orientedFolder, "batch_align.log")
    logger = logging.getLogger()
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    # File handler for the orientedFolder
    file_handler = logging.FileHandler(log_file_path)
    file_handler.setLevel(logging.DEBUG)
    
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)
    
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)
    
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    logging.info(f"Logging to file: {log_file_path}")

    # Loop over and process each valid file in the selected folder
    valid_extensions = ['.obj', '.ply', '.stl']
    for inputFilename in os.listdir(folderPath):
        fn, ext = os.path.splitext(inputFilename)
        if ext.lower() not in valid_extensions:
            logging.info(f"Skipping {inputFilename} due to unsupported extension.")
            continue

        inputFilePath = os.path.join(folderPath, inputFilename)
        sourceModel = slicer.util.loadModel(inputFilePath)
        
        if not sourceModel:
            logging.error(f"Failed to load model: {inputFilePath}")
            continue

        sourceModel.SetName('source')
        logging.info(f'Successfully loaded: {inputFilePath} as "source"')

        # Process the model for alignment
        results = process_model(sourceModel, targetModel, parameters_dict)
        for result in results:
            logging.info(f"Attempt: {result['attempt']}, Fitness: {result['fitness']}, RMSE: {result['rmse']}")

        # Save the aligned model
        save_aligned_model(sourceModel, orientedFolder, inputFilename)

        # Cleanup
        try:
            transformNodes = slicer.util.getNodesByClass('vtkMRMLTransformNode')
            for transformNode in transformNodes:
                slicer.mrmlScene.RemoveNode(transformNode)
            slicer.mrmlScene.RemoveNode(sourceModel)
        except Exception as e:
            logging.error(f"Error during cleanup: {e}")

    logging.info(f'All models aligned and saved in the "{originalFolderName}_Oriented" folder.')
