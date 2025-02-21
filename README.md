# Medial Epicondyle Cropping Protocol

A semi-automated, standardized 3D cropping protocol for analyzing the medial epicondyle of the human humerus, developed to facilitate the study of entheseal changes and skeletal adaptations in biological anthropology. The method ensures high repeatability and reproducibility for defining the region of interest.
<br> 
<br>

**Author and contact details**  
Author &nbsp;&nbsp;&nbsp; *Elle B. K. Liagre*  
Email &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  *elle.liagre@u-bordeaux.fr*  
ORCID <span>&nbsp;&nbsp;&nbsp;</span>   *[https://orcid.org/0000-0002-8993-3266](https://orcid.org/0000-0002-8993-3266)*

<br>

---

## Table of Contents

- [Description](#description)
- [Dependencies and libraries](#dependencies-and-libraries)
- [How to Use](#how-to-use)
- [Data Availability](#data-availability)
- [Contributing](#contributing)
- [License and Citation](#license-and-citation)
<br>

---

## Description

This project introduces a semi-automated, standardized 3D cropping protocol for analyzing the medial epicondyle of the human humerus, with a focus on its entheses. The protocol, implemented in **3DSlicer** using the **PyMeshLab** and **SlicerMorph** extensions, standardizes input models —such as mesh resolution and anatomical orientation—ensuring high repeatability and reproducibility in the analysis. Anatomical landmarks define the boundaries of the medial epicondyle, accounting for variations in humeral size, and ensuring consistent capture of the entire entheseal surface, even amidst anatomical variability.

To improve efficiency, batch processing is incorporated, enabling the protocol to be applied to multiple humeral specimens in a streamlined manner. This feature is particularly beneficial for handling larger datasets, ensuring consistent results across specimens while minimizing manual intervention. The methodology enhances the accuracy and comparability of entheseal surface analysis, offering a valuable tool for studies on skeletal adaptations, physical activity, and human evolution, with strong potential for large-scale applications.

This protocol is presented in the study by **Liagre, E.B.K.; Remy, F.; Knüsel, C.J.; Villotte, S. (Submitted), "A Standardized, Three-Dimensional Cropping Protocol for Analyzing the Medial Epicondyle of the Humerus."** Methodological choices and further information can be found there.
<br>
<br>
<br>

## Dependencies and libraries
This project was built with the following software and library versions. Please ensure these dependencies are installed to maintain compatibility and ensure proper functionality of the protocol. Refer to the [How to Use](#how-to-use) section for installation steps.

### Software and Extensions

- **3D Slicer**: >= 5.6.2
- **SlicerMorph**: >= 5151748
- **Meshlab**: >= 2022.02

### Python Libraries
- **PyMeshLab**: >= 2023.12.post2
- **VTK**: >= 9.2.20230607

### Citations

- *3D Slicer*. (2024, September 9). 3D Slicer image computing platform | 3D Slicer. Retrieved September 16, 2024, from 3D Slicer is a free, open-source software for visualization, processing, segmentation, registration, and analysis of medical, biomedical, and other 3D images and meshes; and planning and navigating image-guided procedures. website: [https://www.slicer.org/](https://www.slicer.org/)
- Fedorov, A., Beichel, R., Kalpathy-Cramer, J., Finet, J., Fillion-Robin, J.-C., Pujol, S., … Kikinis, R. (2012). 3D Slicer as an image computing platform for the Quantitative Imaging Network. *Magnetic Resonance Imaging, 30*(9), 1323–1341. [doi: 10.1016/j.mri.2012.05.001](https://doi.org/10.1016/j.mri.2012.05.001)
- Rolfe, S., Pieper, S., Porto, A., Diamond, K., Winchester, J., Shan, S., … Maga, A. M. (2021). SlicerMorph: An open and extensible platform to retrieve, visualize and analyse 3D morphology. *Methods in Ecology and Evolution, 12*(10), 1816–1825. [doi: 10.1111/2041-210X.13669](https://doi.org/10.1111/2041-210X.13669)
- Porto, A., Rolfe, S., & Maga, A. M. (2021). ALPACA: A fast and accurate computer vision approach for automated landmarking of three‐dimensional biological structures. *Methods in Ecology and Evolution, 12*(11), 2129–2144. [doi: 10.1111/2041-210X.13689](https://doi.org/10.1111/2041-210X.13689)
- Cignoni, P., Callieri, M., Corsini, M., Dellepiane, M., Ganovelli, F., & Ranzuglia, G. (2008). MeshLab: An Open-Source Mesh Processing Tool. *Europgraphics Italian Chapter Conference*, 8.
- Muntoni, A., & Cignoni, P. (2024). *PyMeshLab: PyMeshLab v2023.12.post2*. Zenodo. [doi: 10.5281/zenodo.13768931](https://doi.org/10.5281/zenodo.13768931)
<br>


## How to Use

### Step 1: Install Required Software
1. **3D Slicer**: Download from the [3D Slicer website](https://www.slicer.org).
2. **Meshlab**: Download from the [Meshlab website](https://www.meshlab.net/).  
   *Required for remeshing functions.*
3. **SlicerMorph Extension**: Install via the **Extension Manager** in 3D Slicer. See the [SlicerMorph installation guide](https://slicermorph.github.io/) for details.  
   *Required for orientation functions.*

> **Note**: It’s recommended to restart 3D Slicer after installing these tools.
<br>


### Step 2: Set Up Your 3D Slicer Environment 
1. **Open 3D Slicer.**  
2. By default, 3D Slicer uses the **"Conventional" layout**, which includes a **3D view** and three slice views (**red, green, and yellow**).

   Since the following functions require **only the 3D view**, you may prefer switching to the **"3D only" layout** for better visibility. You can change the layout by:  
   - Clicking on the **"Layout selection" icon** in the toolbar (it mimics the current layout).  
   - Or, navigating to **View** → **Layout** → **3D only** in the menu bar.  

3. **Open the Python Console**, as it will be needed for the next steps. You can do this by:  
   - Clicking the **Python logo** <img src="https://upload.wikimedia.org/wikipedia/commons/c/c3/Python-logo-notext.svg" alt="Python Logo" width="20">  in the toolbar.  
   - Or, going to **View** → **Python Console** in the menu bar.  
   - Or, using the keyboard shortcut **CTRL + 3**.  
<br>


### Step 3: Add Functions to 3D Slicer
1. Open 3D Slicer and go to the python console.
2. Copy your Python script directly into the console, or load a `.py` file by running (replace `'path_to_your_script/your_script.py'` with the actual path to your script):
   ```python
   exec(open(r'path_to_your_script/your_script.py').read())
3. Press **Enter** to execute the script.
<br>


### Step 4: Run functions
To execute each function, enter the function call in the Python console and press **Enter**.

**Note**:  
- Valid input file extensions are: `.obj`, `.ply`, and `.stl`.  
- Output files will be saved in `.ply` format.
- When providing a file or folder path, be sure to use the **raw string format** (i.e., prefix the path with `r''`) to avoid issues with special characters in file paths. For example, use `r'C:\path\to\your\folder'` instead of `'C:\path\to\your\folder'`, otherwise the function may not work correctly.



| Function                                | Description                                                                                                                                                                                                                                                                                   |
|-----------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `batchRemeshing(r'/path/to/your/folder', res = 22)` | Batch function to remesh all 3D models in a folder. Specify a folder path, or if omitted, a window will prompt you to select one. Mesh resolution (faces per mm²) can be set, with a default of 22. A new folder with the remeshed models will be created within the original folder.        |
| `batchMirroring(r'/path/to/your/folder')` | Batch function to mirror all 3D models in a specified folder. If no folder path is provided, a window will prompt you to select one. Mirrored models will be saved in a new folder within the original folder.                                                                                 |
| `batchAlign(r'/path/to/your/folder')` | Batch function to align all 3D models in a folder to a reference model. **Note**: Requires the SlicerMorph extension (see Step 1). Load the reference model in 3D Slicer and name it "target" before running the function. A reference model is provided (see [Data Availability](#data-availability)), or use your own. The aligned models will be saved in a new folder within the original folder. |
| `crop(r'/path/to/your/file')`            | Crops the humeral medial epicondyle as described in the associated article. Specify the model’s file path, or if omitted, a window will prompt you to select one. Instructions are provided when the function runs. **Note**: This function is for right humeri only due to orientation requirements.                                   |
| `batchCrop(r'/path/to/your/folder')`     | Batch function to crop all 3D models in a specified folder. If no folder path is provided, a window will prompt you to select one. Cropped models will be saved in a new folder within the original folder.                                            |
<br>

> **Note**:  
> - **Step 1** (installing software) needs to be done only once.  
> - **Step 2** (adding functions to 3D Slicer) must be repeated each time you open 3D Slicer.
<br>

### Manual Execution of Functions 
   The `Remeshing`, `Mirroring`, and `FastModelAlign` functions can also be run manually for a single model:
   - **Remeshing**: Uses "Remeshing: Isotropic Explicit Remeshing" and "Simplification: Quadric Edge Collapse Decimation" filters in Meshlab.
   - **Mirroring**: Available in 3D Slicer’s **Surface Toolbox** module.
   - **FastModelAlign**: Available through the **SlicerMorph** extension in 3D Slicer (requires download and installation as noted in Step 1). Detailed instructions on how to use this function can be found in the [FastModelAlign tutorial](https://github.com/SlicerMorph/Tutorials/tree/main/FastModelAlign).
<br>


## Data Availability

The reference 3D model used in the alignment process is publicly available for download from the [MorphoSource 3D data repository](https://www.morphosource.org/concern/media/000691022). This model can be used to align new humeral specimens in the same orientation. For detailed instructions on how to use the model, refer to the [How to Use](#how-to-use) section and the associated article.

<br>


## Contributing

We welcome contributions to this project! If you encounter any issues or have ideas for enhancements, please reach out to the corresponding author. Any feedback, questions, or suggestions are welcomed and can be directed to **Elle B. K. Liagre** at *elle.liagre@u-bordeaux.fr*.

<br>


## License and citation

This project is licensed under the GPL 3.0-License - see the [LICENSE](LICENSE) file for details.

Please cite this repository as: Liagre, E.B.K. (2024). Medial Epicondyle Cropping Protocol. https://github.com/ElleLiagre/medial-epicondyle-cropping-protocol 

