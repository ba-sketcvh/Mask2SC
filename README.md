# Protocol for integrating orthogonal information with Xenium data at single-cell resolution
![STAR_AI_Template-01](https://github.com/user-attachments/assets/c0a01923-6ba4-4e3a-a2f1-6cd7ceec8d35)

This includes all the code used in the manuscript "Protocol for integrating orthogonal information with Xenium data at single-cell resolution".

The dataset used in the protocol can be found here:

1. Image Registration 
    - **IHC Channel Separation.groovy**: Extracts IHC channels as individual .tiff files 
    - **Binary Mask Creation.groovy**: Creates a binary mask from annotations    
    

2. Integration
     - **Mean_Signal_Workflow.ipynb**: Calculates the mean intensity of the IHC channel/mask of each cell
     - **Image_Gen_Protocol.ipynb**: Visualizes the integrated data 

