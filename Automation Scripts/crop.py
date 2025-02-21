import slicer
import vtk
import os
import qt
import numpy as np
import logging
import traceback

# Configure logging
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

# Remove existing handlers if they exist to avoid duplicates
for handler in logger.handlers[:]:
    logger.removeHandler(handler)

# Create console handler
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


def crop(file_path: str = None):
    try:
        
        global originalName
        global originalDirectory
        
        # Define constants
        MODEL_NODE_NAME = "model"
        FILE_FILTER = "Model files (*.obj *.ply *.stl)"
        CURVE_NODE_NAME = "TrochlearBorder"
        MODEL_NODE_NAME = "model"
        DISTANCE_FACTOR = 8
        COLOR_YELLOW = [1.0, 1.0, 0.0]
        COLOR_GREY = [0.5, 0.5, 0.5]
        CROP_SUFFIX = "_crop.ply"
        
        if not file_path:
            file_dialog = qt.QFileDialog()
            file_dialog.setNameFilter(FILE_FILTER)
            
            if file_dialog.exec_():
                file_path = file_dialog.selectedUrls()[0].toString()
                file_path = qt.QUrl(file_path).toLocalFile()
            else:
                logging.warning('No file selected.')
                return
            
        originalDirectory = os.path.dirname(file_path)

        # Add file handler globally to log all messages from the start
        log_file_path = os.path.join(originalDirectory, "crop.log")
        file_handler = logging.FileHandler(log_file_path)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
            
            # Load the file
        try:
            modelNode = slicer.util.loadModel(file_path)
            originalName = modelNode.GetName()             
            
            modelNode.SetName(MODEL_NODE_NAME)
                
            logging.info(f'Successfully loaded: {file_path} as "{MODEL_NODE_NAME}" (original name: {originalName}, directory: {originalDirectory})')
                
        except Exception as e:
            logging.error(f'Failed to load model: {e}')
        
        ############# Set the scene #############
        
        modelNode = slicer.util.getNode(MODEL_NODE_NAME)
        
        if modelNode is not None:
            pointListNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode", "LM")
            logging.info(f'Created new point list node: {pointListNode.GetName()}')
        else:
            logging.warning('Model node is None. Cannot create point list.')
        
        
        ############# Defining the anatomical landmarks ################
        
        def set_view_and_confirm(landmark_number, axis, landmark_description):
            layoutManager = slicer.app.layoutManager()
            threeDView = layoutManager.threeDWidget(0).threeDView()
            threeDView.resetFocalPoint()
            threeDView.resetCamera()
            threeDView.rotateToViewAxis(axis)
        
            # Wait for user confirmation
            user_input = input(f"You are about to place landmark {landmark_number}, defined as '{landmark_description}'. "
                               "Please locate this landmark on the model. Click on the 'Place a control point'-icon "
                               "(red ball with blue arrow pointing up) on the toolbar to activate point placement. "
                               "You can move the placed point by clicking and dragging until it is placed in the right location. "
                               "Press spacebar and enter to continue the script: ")
        
            if user_input.strip() == '':
                logger.info("User confirmed to continue. Continuing the script...")
                return True
            else:
                logger.warning("User did not confirm. Stopping the script.")
                return False
        
        ### Placement of LM1
        if set_view_and_confirm(1, 2, 'Medial-most point of the olecranon fossa'):
            pass
        
        ### Placement of LM2
        if set_view_and_confirm(2, 3, 'Viewed anteriorly: proximomedial extreme of the trochlea'):
            pass
        
        ### Placement of LM3
        if set_view_and_confirm(3, 2, 'Viewed posteriorly: proximomedial extreme of the trochlea'):
            pass 
        
        logging.info("All manual landmarks have been placed. Continuing with automatic placement of landmarks.")
        

        ### Automatic Placement of Landmarks (LM4 to LM6)

        def automaticLandmarkPlacement():

            ## Getting the maximal values of the model bounding box
            polydata = modelNode.GetPolyData()
            points = polydata.GetPoints()
        
            # Initialize max and min values to negative and positive infinity, respectively
            max_x = max_y = max_z = float('-inf')
            min_x = min_y = min_z = float('inf')
        
            # Iterate through each point (vertex) in the PolyData
            for i in range(points.GetNumberOfPoints()):
                vertex = points.GetPoint(i)
        
                max_x = max(max_x, vertex[0])
                max_y = max(max_y, vertex[1])
                max_z = max(max_z, vertex[2])
        
                min_x = min(min_x, vertex[0])
                min_y = min(min_y, vertex[1])
                min_z = min(min_z, vertex[2])
        
            # Log bounding box limits
            logging.info("Max X: %s", max_x)
            logging.info("Max Y: %s", max_y)
            logging.info("Max Z: %s", max_z)
            logging.info("Min X: %s", min_x)
            logging.info("Min Y: %s", min_y)
            logging.info("Min Z: %s", min_z)
        
            
            ## Find most distal point (based on Min Z)
            # Calculate z-range
            tolerance = 1e-6
            z_coordinate_range = (min_z - tolerance, min_z + tolerance)
        
            # Find all points with z-coordinate within the specified range
            min_z_matching_points = []
            for i in range(points.GetNumberOfPoints()):
                vertex = points.GetPoint(i)
                if z_coordinate_range[0] <= vertex[2] <= z_coordinate_range[1]: 
                    min_z_matching_points.append(vertex)
        
            # Choose a point
            if min_z_matching_points:
                chosen_z = max(min_z_matching_points, key=lambda p: p[1])
                logging.info("Defined LM4.")
            else:
                logging.warning("No points within the specified Z-coordinate range found.")


            ## Find most anterior and posterior points
            # Retrieve LM2 from the point list by index
            coords_P2 = np.array(pointListNode.GetNthControlPointPosition(1))
            P2_x = coords_P2[0]
            
            # Retrieve LM3 from the point list by index
            coords_P3 = np.array(pointListNode.GetNthControlPointPosition(2))
            P3_x = coords_P3[0]
        
            # Calculate ranges
            x_range_LM5 = (min(chosen_z[0], P2_x) - 5.0, max(chosen_z[0], P2_x) + 1.0)
            x_range_LM6 = (min(chosen_z[0], P3_x), max(chosen_z[0], P3_x) + 1.0)
            z_range = (chosen_z[2], coords_P3[2])
        
            # Subset points that fell within the x and z ranges
            polyData = modelNode.GetPolyData()
            filteredPoints_LM5 = vtk.vtkPoints()
            filteredPoints_LM6 = vtk.vtkPoints()
            
            # Initialize variables to track the max and min y-values
            max_y = float('-inf')
            min_y = float('inf')
        
            # Iterate through the points in the model
            total_points_checked = 0
            for i in range(polyData.GetNumberOfPoints()):
                point = polyData.GetPoint(i)
                x, y, z = point
        
                # Count the total points checked
                total_points_checked += 1
        
                # Check if the point is within the specified ranges
                if x_range_LM5[0] <= x <= x_range_LM5[1] and z_range[0] <= z <= z_range[1]:
                    filteredPoints_LM5.InsertNextPoint(point)
            
                    if y > max_y:
                        max_y = y
            
            # Iterate through the points in the model
            total_points_checked = 0
            for i in range(polyData.GetNumberOfPoints()):
                point = polyData.GetPoint(i)
                x, y, z = point
        
                # Count the total points checked
                total_points_checked += 1
                
                # Check if the point is within the specified ranges
                if x_range_LM6[0] <= x <= x_range_LM6[1] and z_range[0] <= z <= z_range[1]:
                    filteredPoints_LM6.InsertNextPoint(point)

                    if y < min_y:
                        min_y = y
        
            # Create a new polydata to hold the filtered points
            filtered_LM5PolyData = vtk.vtkPolyData()
            filtered_LM5PolyData.SetPoints(filteredPoints_LM5)
            filtered_LM6PolyData = vtk.vtkPolyData()
            filtered_LM6PolyData.SetPoints(filteredPoints_LM6)
            
            # Check the number of points in the filtered polydata
            num_filtered_points_LM5 = filtered_LM5PolyData.GetNumberOfPoints()
            logging.info('Number of points in the anteriorly filtered polydata: %d', num_filtered_points_LM5)

            num_filtered_points_LM6 = filtered_LM6PolyData.GetNumberOfPoints()
            logging.info('Number of points in the posteriorly filtered polydata: %d', num_filtered_points_LM6)
        
            logging.info('Total points checked: %d', total_points_checked)
        
            # Log the max and min y-values found
            if num_filtered_points_LM5 > 0:
                logging.info('Max y-value among anteriorly filtered points: %f', max_y)
                logging.info('Min y-value among anteriorly filtered points: %f', min_y)
            else:
                logging.info('No points were found in the specified ranges.')
            
            # Log the max and min y-values found
            if num_filtered_points_LM6 > 0:
                logging.info('Max y-value among posteriorly filtered points: %f', max_y)
                logging.info('Min y-value among posteriorly filtered points: %f', min_y)
            else:
                logging.info('No points were found in the specified ranges.')
              

            ## Place LM4: Distomedial extreme of the trochlea
            pointListNode.AddControlPoint([chosen_z[0], chosen_z[1], chosen_z[2]])
            logging.info("Added LM4 to point list node.")
        

            ## Place LM5: Viewed medially: Anterior extreme of the medial edge/ridge of the trochlea
            # Set a range for y-coordinate comparison, but using max_y this time
            y_coordinate_range_max = (max_y - tolerance, max_y + tolerance)

            # Find all points with y-coordinate within the specified range within the subset
            max_y_matching_points = []
            for i in range(filtered_LM5PolyData.GetNumberOfPoints()):
                vertex = filtered_LM5PolyData.GetPoint(i)
                if y_coordinate_range_max[0] <= vertex[1] <= y_coordinate_range_max[1]:
                    max_y_matching_points.append(vertex)
        
            # Choose a point
            if max_y_matching_points:
                chosen_y_max = max(max_y_matching_points, key=lambda p: p[1])
                pointListNode.AddControlPoint([chosen_y_max[0], chosen_y_max[1], chosen_y_max[2]])
                logging.info("Added LM5 to point list node.")
            else:
                logging.warning("No points within the specified Y-coordinate range found. Not adding LM5.")
        

            ## Place LM6: Viewed medially: Posterior extreme of the medial edge/ridge of the trochlea
            # Set a range for y-coordinate comparison, but using min_y this time
            y_coordinate_range_min = (min_y - tolerance, min_y + tolerance)
        
            # Find all points with y-coordinate within the specified range within the subset
            min_y_matching_points = []
            for i in range(filtered_LM6PolyData.GetNumberOfPoints()):
                vertex = filtered_LM6PolyData.GetPoint(i)
                if y_coordinate_range_min[0] <= vertex[1] <= y_coordinate_range_min[1]:
                    min_y_matching_points.append(vertex)
        
            # Choose a point
            if min_y_matching_points:
                chosen_y_min = min(min_y_matching_points, key=lambda p: p[1])
                pointListNode.AddControlPoint([chosen_y_min[0], chosen_y_min[1], chosen_y_min[2]])
                logging.info("Added LM6 to point list node.")
            else:
                logging.warning("No points within the specified Y-coordinate range found. Not adding LM6.")
        
            logging.info("Automatic landmark placement completed.")
        
        # Call the automatic placement function
        automaticLandmarkPlacement()
        
        ### Overview all landmarks
        # Set 3D view to the specified axis
        layoutManager = slicer.app.layoutManager()
        threeDView = layoutManager.threeDWidget(0).threeDView()
        threeDView.resetFocalPoint()
        threeDView.resetCamera()
        threeDView.rotateToViewAxis(0)
        
        # Wait for user confirmation in the console
        user_input = input("All landmarks have been placed. Please verify if they have been placed accurately. Press spacebar and enter to continue the script:")
        
        # Check if the user confirmed
        if user_input.lower() == ' ':
            logging.info("User confirmed to continue. Continuing the script...")
        else:
            logging.info("User did not confirm. Stopping the script.")
        

        ############# Preparation for cropping in markups module #############
        
        ### Calculate buffer zone

        def tricircumcenter3d(a, b, c):
            """
            Find the circumcenter of a triangle in 3D space.
        
            Parameters:
            a, b, c: 3D points (list or numpy array) defining the triangle.
        
            Returns:
            circumcenter: 3D coordinates of the circumcenter.
            xi: Local xi-coordinate relative to triangle point 'a'.
            eta: Local eta-coordinate relative to triangle point 'a'.
            distances: Distances from circumcenter to points a, b, and c.
            """
            # Convert inputs to numpy arrays for easier calculations
            a, b, c = np.array(a), np.array(b), np.array(c)
        
            # Vectors relative to point `a`
            vector_ab = b - a
            vector_ac = c - a
        
            # Squared lengths of the edges incident to `a`
            ab_length_squared = np.dot(vector_ab, vector_ab)
            ac_length_squared = np.dot(vector_ac, vector_ac)
        
            # Cross product of the edges
            cross_product = np.cross(vector_ab, vector_ac)
        
            # Calculate the denominator for the circumcenter formula
            denominator = 0.5 / np.dot(cross_product, cross_product)
        
            # Calculate the circumcenter's offset from point `a`
            x_offset = ((ab_length_squared * vector_ac[1] - ac_length_squared * vector_ab[1]) * cross_product[2] -
                (ab_length_squared * vector_ac[2] - ac_length_squared * vector_ab[2]) * cross_product[1]) * denominator
            y_offset = ((ab_length_squared * vector_ac[2] - ac_length_squared * vector_ab[2]) * cross_product[0] -
                (ab_length_squared * vector_ac[0] - ac_length_squared * vector_ab[0]) * cross_product[2]) * denominator
            z_offset = ((ab_length_squared * vector_ac[0] - ac_length_squared * vector_ab[0]) * cross_product[1] -
                (ab_length_squared * vector_ac[1] - ac_length_squared * vector_ab[1]) * cross_product[0]) * denominator
 
            circumcenter = a + np.array([x_offset, y_offset, z_offset])
 
            # Calculate distances from circumcenter to each vertex
            distances = np.array([
                np.linalg.norm(circumcenter - a),
                np.linalg.norm(circumcenter - b),
                np.linalg.norm(circumcenter - c)
            ])
        
            return circumcenter, x_offset, y_offset, distances
        
        def get_coords(pointListNode, landmark_id):
            if not pointListNode:
                logger.error("Point list node not found.")
                return None
        
            # Find the index for the specified landmark
            index = pointListNode.GetControlPointIndexByLabel(landmark_id)
            if index == -1:
                logger.error(f"Landmark '{landmark_id}' not found in point list node.")
                return None
        
            # Get the fiducial position
            coords = np.array(pointListNode.GetNthControlPointPosition(index))
            return coords
        
        
        ## Calculating diameter circumcircle

        if not pointListNode:
            logger.error("Point list node not found.")
        else:
            coords_LM4 = get_coords(pointListNode, "LM-4")
            coords_LM5 = get_coords(pointListNode, "LM-5")
            coords_LM6 = get_coords(pointListNode, "LM-6")

        
            # Calculate circumcenter and distances
            circumcenter, xi, eta, distances = tricircumcenter3d(coords_LM4, coords_LM5, coords_LM6)
            logger.info(f"Circumcenter: {circumcenter}")
            logger.info(f"Distance from circumcenter to points LM4, LM5, LM6: {distances[1]}")
        
            # Diameter of the circle defined by the points LM4, LM5, and LM6
            diam = distances[1] * 2
            logger.info(f"Diameter of the circumcircle: {diam}")
        
        
        ### Defining cutting plane
        # Constants for landmark IDs
        LANDMARK_IDS = ['LM-1', 'LM-2', 'LM-3']
        
        def create_plane_from_landmarks(pointListNode):
            if not pointListNode:
                logger.error("Point list node not found.")
                return
        
            landmark_coords = []
            
            # Get coordinates of necessary landmarks
            for landmark_id in LANDMARK_IDS:
                index = pointListNode.GetControlPointIndexByLabel(landmark_id)
                if index == -1:
                    logger.error(f"Landmark '{landmark_id}' not found in point list node.")
                    return
                
                coords = np.array(pointListNode.GetNthControlPointPosition(index))
                landmark_coords.append(coords)
                logger.info(f"Retrieved coordinates for Landmark '{landmark_id}': {coords}")
        
            # Create a plane node and add points
            planeNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsPlaneNode')
            planeNode.SetPlaneType(0)
        
            for coords in landmark_coords:
                planeNode.AddControlPoint(coords)
        
            logger.info("Plane created successfully with the provided landmarks.")
            slicer.app.processEvents()
        
        # Main execution
        if __name__ == "__main__":
            create_plane_from_landmarks(pointListNode)
        
        
        ### Defining trochlear border
        ## Creating curve
        
        # Get landmark coords
        if not pointListNode:
            logger.error("Point list node not found.")
        else:
            coords_LM1 = get_coords(pointListNode, 'LM-1')
            coords_LM2 = get_coords(pointListNode, 'LM-2')
            coords_LM3 = get_coords(pointListNode, 'LM-3')
        
        
        # Main logic to create and configure the trochlear border curve
        try:
            if pointListNode is None:
                logger.error("Point list node not found. Aborting script.")
                raise ValueError("Point list node not found. Aborting script.")
        
            # Create a new markups curve node
            curveNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsCurveNode")
            curveNode.SetName(CURVE_NODE_NAME)
            logger.info("Curve node created: %s", CURVE_NODE_NAME)
        
            # Add control points to the curve node
            coords_LM = [coords_LM2, coords_LM5, coords_LM4, coords_LM6, coords_LM3]
            for point in coords_LM:
                curveNode.AddControlPointWorld(point)
            logger.info("Added control points to curve node: %s", CURVE_NODE_NAME)
        
            # Get the model node
            model = slicer.util.getNode(MODEL_NODE_NAME)
            if model is None:
                logger.error("Model node '%s' not found. Aborting script.", MODEL_NODE_NAME)
                slicer.mrmlScene.RemoveNode(curveNode)
                raise ValueError("Model node not found. Aborting script.")
        
            # Set the model as a surface constraint for the curve node
            curveNode.SetAndObserveSurfaceConstraintNode(model)
            logger.info("Surface constraint set for curve node: %s", CURVE_NODE_NAME)
        
            # Update the display to show the curve
            slicer.app.processEvents()
            logger.info("Trochlear border curve created successfully.")
        
        except Exception as e:
            logger.exception("An error occurred during the process: %s", e)
        
        
        ## Resample according to diameter circumscribed circle
        """ 
        The distance between sample points should be smaller than the smallest radius(r1) used for the 
        cropping of the trochlear border. It ensures that there is sufficient overlap between buffers later on. 
        This value (diam) is calculated as the diameter of the circumscribed distance / 6. 
        Therefore the distance between the sample points will be defined arbitrarily as diam / 8. 
        """
 
        distanceBetweenControlPoints = diam / DISTANCE_FACTOR 
        
        try:
            logging.info("Starting curve resampling...")
            curveNode.ResampleCurveWorld(distanceBetweenControlPoints)
            logging.info("Curve resampling completed successfully.")
        except Exception as e:
            slicer.mrmlScene.RemoveNode(curveNode)
            logging.error(f"Error during curve resampling: {str(e)}")
            raise RuntimeError(f"Error during curve resampling: {str(e)}")
 
 
        ## Creating points lists for the Trochlear border
        # Trochlear border A
        if curveNode is None:
            logging.error("Curve node not found.")
        else:
            logging.info("Creating Trochlear Border A fiducial node...")
            trochBorderA = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
            trochBorderA.SetName("TrochlearBorderA")
        
            for i in range(curveNode.GetNumberOfControlPoints()):
                point = [0, 0, 0]
                try:
                    curveNode.GetNthControlPointPositionWorld(i, point)
                    trochBorderA.AddControlPointWorld(vtk.vtkVector3d(*point))
                    logging.debug(f"Control point {i + 1} added to Trochlear Border A.")
                except Exception as e:
                    logging.error(f"Error adding control point {i + 1} to Trochlear Border A: {str(e)}")
        
            trochBorderA.SetAttribute('CurveSourceNodeID', curveNode.GetID())
            logging.info("Trochlear Border A fiducial node created successfully.")
        
        # Trochlear border B
        if curveNode is None:
            logging.error("Curve node not found.")
        else:
            logging.info("Creating Trochlear Border B fiducial node...")
            trochBorderB = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
            trochBorderB.SetName("TrochlearBorderB")
        
            # Define the name range for the subset
            index_T2 = curveNode.GetControlPointIndexByLabel("TrochlearBorder-2")
            index_T3 = curveNode.GetControlPointIndexByLabel("TrochlearBorder-3")
        
            for i in range(index_T2, min(index_T3 + 1, curveNode.GetNumberOfControlPoints())):
                point = [0, 0, 0]
                try:
                    curveNode.GetNthControlPointPositionWorld(i, point)
                    trochBorderB.AddControlPointWorld(vtk.vtkVector3d(*point))
                    logging.debug(f"Control point {i + 1} added to Trochlear Border B.")
                except Exception as e:
                    logging.error(f"Error adding control point {i + 1} to Trochlear Border B: {str(e)}")
        
            trochBorderB.SetAttribute('CurveSourceNodeID', curveNode.GetID())
            logging.info("Trochlear Border B fiducial node created successfully.")
        
        
        ############# Cropping of model #############
        
        ### Execute first crop of the model
        modelNode = slicer.util.getNode("model")
        planeNode = slicer.util.getNodesByClass("vtkMRMLMarkupsPlaneNode")
        planeNode = planeNode[0]
        
        def execute_plane_cut(model, plane_node, operation_type="Union"):
            
            # Check and retrieve nodes
            if isinstance(model, str):
                model = slicer.util.getNode(model)
            if isinstance(plane_node, str):
                plane_node = slicer.util.getNodesByClass("vtkMRMLMarkupsPlaneNode")
                plane_node = plane_node[0]
        
            if not model or not plane_node:
                raise ValueError("Model or Plane Node not found. Please check the names provided.")
        
            # Create new models for output
            positive_model = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode")
            negative_model = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode")
        
            # Set up the plane cut tool
            plane_cut = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLDynamicModelerNode")
            plane_cut.SetToolName("Plane cut")
            plane_cut.SetAttribute("CapSurface", "0")
            plane_cut.SetNodeReferenceID("PlaneCut.InputModel", model.GetID())
            plane_cut.SetNodeReferenceID("PlaneCut.InputPlane", plane_node.GetID())
            plane_cut.SetNodeReferenceID("PlaneCut.OutputPositiveModel", positive_model.GetID())
            plane_cut.SetNodeReferenceID("PlaneCut.OutputNegativeModel", negative_model.GetID())
            plane_cut.SetAttribute("OperationType", operation_type) 
        
            # Run the dynamic modeler tool
            model.SetDisplayVisibility(1)
            slicer.modules.dynamicmodeler.logic().RunDynamicModelerTool(plane_cut)
        
            positive_model.SetName("firstCrop")
            slicer.mrmlScene.RemoveNode(negative_model)

            # Log success message
            logging.info("Plane cut executed successfully. Positive model: %s", positive_model.GetName())

            return positive_model, negative_model
        
        positive_model, negative_model = execute_plane_cut(modelNode, planeNode, operation_type="Difference")
        
        
        ### Execute second crop of the model
        class TrochleaProcessor:
            def __init__(self):
                pass
        
            def process_crop(self, node_name, troch_border_name, distance, output_name):
                try:
                    crop = self.get_node(node_name)
                    troch_border = self.get_node(troch_border_name)
        
                    logging.info(f"Processing crop for nodes: {node_name}, {troch_border_name}")
        
                    output_scalars = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode")
                    self.setup_select_points_tool(crop, troch_border, output_scalars, distance)
        
                    logging.info(f"Dynamic modeler tool run successfully for {output_name}")
        
                    output_scalars.SetName(output_name)
                    model_node = self.get_node(output_name)
        
                    thresholded_model = self.create_thresholded_model(model_node, 0.5)
                    thresholded_model.SetName(f"{output_name}Scalars")
        
                    self.reset_scalars(thresholded_model)
                    self.copy_normals(model_node.GetPolyData(), thresholded_model.GetPolyData())
        
                    logging.info(f"Processing completed successfully for {output_name}")
        
                except Exception as e:
                    logging.error(f"Error: {str(e)}")
                    logging.error(traceback.format_exc())
        
            def setup_select_points_tool(self, crop, troch_border, output_scalars, distance):
                select_points = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLDynamicModelerNode")
                select_points.SetToolName("Select by points")
                select_points.SetNodeReferenceID("SelectByPoints.InputModel", crop.GetID())
                select_points.SetNodeReferenceID("SelectByPoints.InputFiducial", troch_border.GetID())
                select_points.SetNodeReferenceID("SelectByPoints.SelectionScalarsModel", output_scalars.GetID())
                select_points.SetAttribute("SelectionDistance", str(distance))
                crop.SetDisplayVisibility(1)
                slicer.modules.dynamicmodeler.logic().RunDynamicModelerTool(select_points)
        
            def get_node(self, node_name):
                node = slicer.util.getNode(node_name)
                if not node:
                    raise ValueError(f"Node '{node_name}' not found.")
                return node
        
            def create_thresholded_model(self, input_model, lower_threshold):
                threshold = vtk.vtkThreshold()
                threshold.SetInputConnection(input_model.GetPolyDataConnection())
                threshold.SetThresholdFunction(threshold.THRESHOLD_LOWER)
                threshold.SetLowerThreshold(lower_threshold)
                threshold.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, "Selection")
                threshold.Update()
        
                surface = vtk.vtkGeometryFilter()
                surface.SetInputConnection(threshold.GetOutputPort())
                surface.Update()
        
                return slicer.modules.models.logic().AddModel(surface.GetOutput())
        
            def reset_scalars(self, model_node):
                poly_data = model_node.GetPolyData()
                if poly_data:
                    point_data = poly_data.GetPointData()
                    if point_data:
                        num_arrays = point_data.GetNumberOfArrays()
                        for i in range(num_arrays - 1, -1, -1):
                            point_data.RemoveArray(i)
                        poly_data.Modified()
        
                display_node = model_node.GetDisplayNode()
                if display_node:
                    display_node.SetScalarVisibility(0)
        
            def copy_normals(self, source_poly_data, target_poly_data):
                source_normals = source_poly_data.GetPointData().GetNormals()
                if source_normals:
                    target_poly_data.GetPointData().SetNormals(source_normals)
        
        # Execution: removal of trochlea
        if __name__ == "__main__":
            processor = TrochleaProcessor()
            
            processor.process_crop("firstCrop", "TrochlearBorderA", diam / 6, "secondCrop")

            processor.process_crop("secondCropScalars", "TrochlearBorderB", diam / 4, "thirdCrop")
        
        
        ## Keep largest component
        # Get model
        modelNode = slicer.util.getNode("thirdCropScalars")
        if modelNode:
            logger.info("Model node 'thirdCropScalars' retrieved successfully.")
        else:
            logger.error("Model node 'thirdCropScalars' not found.")
            raise RuntimeError("Model node 'thirdCropScalars' not found.")
        
        # Create new empty model for cropped medial epicondyle
        finalCrop = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelNode', 'finalCrop')
        finalCrop.CreateDefaultDisplayNodes()
        finalCrop.GetDisplayNode().SetVisibility2D(True)
        logger.info("New model node 'finalCrop' created.")
        
        # Perform "keep largest connected component" operation on the model
        logic = slicer.util.getModuleLogic('SurfaceToolbox')
        if logic:
            logic.extractLargestConnectedComponent(modelNode, finalCrop)
            logger.info("Extracted the largest connected component.")
        else:
            logger.error("SurfaceToolbox logic could not be retrieved.")
            raise RuntimeError("SurfaceToolbox logic not found.")
        
        # Clean up
        slicer.app.processEvents()
        cleanFilter = vtk.vtkCleanPolyData()
        cleanFilter.SetInputData(finalCrop.GetPolyData())
        cleanFilter.Update()
        logger.info("Unused vertices cleaned up.")
        
        finalCrop.SetAndObservePolyData(cleanFilter.GetOutput())
        
        
        ##  Keep the second largest component
        # Get the model node by name
        modelNode = slicer.util.getNode("thirdCropScalars")
        
        # Create function to get the second largest component
        def keep_second_largest_component(modelNode):    
            inputPolyData = modelNode.GetPolyData()
        
            # Apply the connectivity filter to extract connected components
            connectivityFilter = vtk.vtkPolyDataConnectivityFilter()
            connectivityFilter.SetInputData(inputPolyData)
            connectivityFilter.SetExtractionModeToAllRegions()
            connectivityFilter.ColorRegionsOn()
            connectivityFilter.Update()
            
            # Get the number of connected regions
            numberOfRegions = connectivityFilter.GetNumberOfExtractedRegions()
            logger.info(f"Number of connected regions: {numberOfRegions}")
            
            regionSizes = [] 

            if numberOfRegions == 1:
                logger.info("Only one connected component found. Using the only component available as the second largest.")
                secondLargestRegionIndex = 0

            elif numberOfRegions >= 2:
                for regionIndex in range(numberOfRegions):
                    # Set the filter to extract one region at a time
                    connectivityFilter.SetExtractionModeToSpecifiedRegions()
                    connectivityFilter.InitializeSpecifiedRegionList()
                    connectivityFilter.AddSpecifiedRegion(regionIndex)
                    connectivityFilter.Update()
                    
                    # Get the output data for the current region
                    regionPolyData = vtk.vtkPolyData()
                    regionPolyData.DeepCopy(connectivityFilter.GetOutput())
                    
                    # Calculate the size of the region (number of points)
                    regionSize = regionPolyData.GetNumberOfPoints()
                    regionSizes.append((regionSize, regionIndex))    
                # Sort the regions by size in descending order
                regionSizes.sort(reverse=True, key=lambda x: x[0])

                # Extract the second largest region
                if len(regionSizes) > 1:
                    secondLargestRegionIndex = regionSizes[1][1]
                else:
                    logger.warning("Only one connected region found, using the only region as the second largest.")

            else:
                raise ValueError("Unexpected number of regions")
        
            # Set the filter to extract only the second-largest region
            connectivityFilter.SetExtractionModeToSpecifiedRegions()
            connectivityFilter.InitializeSpecifiedRegionList()
            connectivityFilter.AddSpecifiedRegion(secondLargestRegionIndex)
            connectivityFilter.Update()
            
            # Create a new polydata to store the second-largest component
            secondLargestPolyData = vtk.vtkPolyData()
            secondLargestPolyData.DeepCopy(connectivityFilter.GetOutput())
            
            # Create a new model node to store the second-largest component
            secondLargestModelNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelNode', 'finalCrop_2')
            secondLargestModelNode.SetAndObservePolyData(secondLargestPolyData)
            logger.info("Second-largest component succesfully extracted.")
        
            # Update the display properties (optional)
            secondLargestModelNode.CreateDefaultDisplayNodes()
        
            # Ensure the display node does not expect any scalar visibility
            display_node = secondLargestModelNode.GetDisplayNode()
            if display_node:
                display_node.SetScalarVisibility(0)
        
        # Extract second largest component
        try:
            keep_second_largest_component(modelNode)
        except Exception as e:
            logger.error(f"An error occurred while extracting the second largest component: {e}")
        
        # Get the model node by name
        modelNode = slicer.util.getNode("finalCrop_2")
        
        # Clean up unused vertices
        cleanFilter = vtk.vtkCleanPolyData()
        cleanFilter.SetInputData(modelNode.GetPolyData())
        cleanFilter.Update()
        
        modelNode.SetAndObservePolyData(cleanFilter.GetOutput())
        
        
        ############# Exporting the cropped model #############
        ## Export cropped model: function
        
        # Function to export the model node as a .ply file and the point list node as .fcsv
        def exportModelAsPLY(modelNode, pointListNode):
            global originalName
            global originalDirectory
        
            if modelNode:

                exportedFilePath = os.path.join(originalDirectory, originalName + CROP_SUFFIX)

                try:
                    slicer.util.saveNode(modelNode, exportedFilePath)
                    logger.info('Successfully exported model as: %s (original name: %s)', exportedFilePath, originalName)
                except Exception as e:
                    logger.error('Failed to export model as PLY: %s', str(e))
                
                if pointListNode:
                    exportedFilePath_fcsv = os.path.join(originalDirectory, originalName + ".fcsv")
                    try:
                        slicer.util.saveNode(pointListNode, exportedFilePath_fcsv)
                        logger.info('Successfully exported points as: %s', exportedFilePath_fcsv)
                    except Exception as e:
                        logger.error('Failed to export points as FCSV: %s', str(e))
                else:
                    logger.warning("Point list node 'pointListNode' not found.")
            else:
                logger.warning("Model node 'model' not found.")
        
        # Get the full model and final crop model nodes
        fullmodelNode = slicer.util.getNode("model")
        modelNode = slicer.util.getNode("finalCrop")
        
        # Function to set visibility and color for nodes
        def setNodeVisibilityAndColor(node, color):
            if node and node.IsA('vtkMRMLModelNode') and node.GetDisplayNode():
                node.GetDisplayNode().SetVisibility(True)
                node.GetDisplayNode().SetColor(color)
            else:
                logger.warning("Node '%s' not found or does not have a display node.", node.GetID() if node else 'unknown')
        
        # Iterate over all nodes in the scene and hide displayable nodes
        for node in slicer.mrmlScene.GetNodes():
            if node.IsA('vtkMRMLDisplayableNode') and node.GetDisplayNode():
                node.GetDisplayNode().SetVisibility(False)
        
        # Set visibility and colors for the crop model and full model
        setNodeVisibilityAndColor(modelNode, COLOR_YELLOW)
        setNodeVisibilityAndColor(fullmodelNode, COLOR_GREY)
        
        ## Model selection and export
        # Wait for user confirmation in the console
        user_input = input("Please verify if the correct model is selected for export. Press 'yes' to continue to export the visible model. Press 'no' if you want to examine an alternative model.")
        
        if user_input.lower() == 'yes':
            logging.info("The selected crop will be exported.")
            modelNode = slicer.util.getNode("finalCrop")
            exportModelAsPLY(modelNode, pointListNode)
            slicer.mrmlScene.Clear()
        else:

            fullmodelNode = slicer.util.getNode("model")
            modelNode = slicer.util.getNode("finalCrop_2")
        
            # Hide all displayable nodes
            for node in slicer.mrmlScene.GetNodes():
                if node.IsA('vtkMRMLDisplayableNode') and node.GetDisplayNode():
                    node.GetDisplayNode().SetVisibility(False)
        
            # Set visibility and color for the alternative crop model
            setNodeVisibilityAndColor(modelNode, COLOR_YELLOW)
        
            if modelNode:
                modelNode.GetDisplayNode().SetColor(COLOR_YELLOW)
            else:
                logging.warning("Final crop model node not found or does not have a display node.")
        
            # Set visibility for the full model node
            setNodeVisibilityAndColor(fullmodelNode, COLOR_GREY)
        
            if not fullmodelNode:
                logging.warning("Full model node 'model' not found or does not have a display node.")
        
            # Ask for confirmation again
            user_input = input("Please verify if the correct model is selected for export. Press 'yes' to continue to export the visible model. Press 'no' to abort the script.")
            
            if user_input.lower() == 'yes':
                logging.info("The selected crop will be exported.")
                modelNode = slicer.util.getNode("finalCrop_2")
                exportModelAsPLY(modelNode, pointListNode)
                slicer.mrmlScene.Clear()
            else:
                logging.info("User did not confirm. Stopping the script.")
 
    except Exception as e:
        logging.error(f"Error occurred: {str(e)}")
        logging.error(traceback.format_exc())
 

