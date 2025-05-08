# Protocol for integrating orthogonal information with Xenium data at single-cell resolution
This includes all the code used in the manuscript "Protocol for integrating orthogonal information with Xenium data at single-cell resolution".

The dataset used in the protocol can be found here:

1. Image Registration 
    - **IHC Channel Separation.groovy**: Extracts IHC channels as individual .tiff files 
    - **Binary Mask Creation.groovy**: Creates a binary mask from annotations    
    

2. Integration
     - **Mean_Signal_Workflow.ipynb**: Calculates the mean intensity of the IHC channel/mask of each cell
     - **Image_Gen_Protocol.ipynb**: Visualizes the integrated data 

