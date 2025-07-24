import qupath.lib.images.servers.ImageServer
import qupath.lib.regions.RegionRequest
import qupath.lib.images.servers.PixelType
import java.awt.image.BufferedImage
import javax.imageio.ImageIO
import java.io.File

// Define the channels to extract (zero-based indexing: 5 = 6th channel, 6 = 7th channel)
def channelsToExtract = [5, 6]

def imageData = getCurrentImageData()
def server = imageData.getServer()

def outputDir = buildFilePath(PROJECT_BASE_DIR, "Exported_Channels")
mkdirs(outputDir)

def width = server.getWidth()
def height = server.getHeight()

for (channelIndex in channelsToExtract) {
    println "Extracting channel ${channelIndex + 1}..."

    def region = RegionRequest.createInstance(server.getPath(), 1.0, 0, 0, width, height)

    def imgChannel = server.readRegion(region).getRaster().getSamples(0, 0, width, height, channelIndex, (double[]) null)

    def img = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY)
    def raster = img.getRaster()

    raster.setPixels(0, 0, width, height, imgChannel)

    def outputFile = new File(outputDir, "Channel_${channelIndex + 1}.tiff")

    ImageIO.write(img, "TIFF", outputFile)

    println "Saved Channel ${channelIndex + 1} to ${outputFile.getAbsolutePath()}"
}
