import qupath.lib.images.servers.LabeledImageServer
import qupath.lib.common.GeneralTools

// Get current image data
def imageData = getCurrentImageData()
def server = imageData.getServer()

// Get original image dimensions
int width = server.getWidth()
int height = server.getHeight()

// Get the original pixel size (in microns)
double originalPixelSize = server.getPixelCalibration().getAveragedPixelSize()

// Define output folder
def outputFolder = buildFilePath(PROJECT_BASE_DIR, 'export')
mkdirs(outputFolder)  // Ensure the folder exists

// Define output TIFF file path explicitly
def name = GeneralTools.getNameWithoutExtension(server.getMetadata().getName())
def pathOutput = buildFilePath(outputFolder, name + "_binary_mask.tif")

// Set the downsample factor to maintain full resolution
double downsample = 1

// Create an ImageServer for binary mask extraction
def labelServer = new LabeledImageServer.Builder(imageData)
    .backgroundLabel(0, 0)  // Background is black (0)
    .downsample(downsample)  // Full resolution
    .addLabel('Tumor', 255)  // Tumor is white (255)
    .multichannelOutput(false)  // Single-channel binary image
    .build()

// Save the binary mask as a TIFF file
writeImage(labelServer, pathOutput)