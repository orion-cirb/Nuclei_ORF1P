package Julia_Nucleus_Analysis;


import Tools.Jnucleus_Tools3D;

import ij.IJ;
import ij.ImagePlus;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Object3D;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.measurements.Measure2Colocalisation;
import mcib3d.geom2.measurements.MeasureIntensity;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.image3d.ImageHandler;
import org.apache.commons.io.FilenameUtils;


public class YH2AX implements PlugIn {
    
    Jnucleus_Tools3D tools = new Jnucleus_Tools3D();
    
   private final boolean canceled = false;
    private String imageDir = "";
    public  String outDirResults = "";
    public  String rootName = "";
    public BufferedWriter nucleus_Analyze;
    
    @Override
    public void run(String arg) {
        try {
            if (canceled) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            imageDir = IJ.getDirectory("Choose Directory Containing Image Files...");
            if (imageDir == null) {
                return;
            }
            
            // Find images with nd extension
            ArrayList<String> imageFile = tools.findImages(imageDir, "nd");
            if (imageFile == null) {
                IJ.showMessage("Error", "No images found with nd extension");
                return;
            }
            // create output folder
            outDirResults = imageDir + File.separator+ "Results"+ File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            
            // create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            reader.setId(imageFile.get(0));
            
            // Find channel names
            List<String> channels = tools.findChannels(imageFile.get(0));
            
            // Find image calibration
            Calibration cal = tools.findImageCalib(meta);


            // write headers
            String header= "Image Name\t# Nucleus\tNucleus Volume\tnucleus Intensity\t# Dots\tDots Nucleus intensity\tDots volume\n";
            FileWriter fwNucleus = new FileWriter(outDirResults + "nucleus_Results.xls", false);
            nucleus_Analyze = new BufferedWriter(fwNucleus);
            nucleus_Analyze.write(header);
            nucleus_Analyze.flush();
            
            
            // Channels dialog
            List<String> chs = new ArrayList();
            List<String> channelsName = new ArrayList();
            channelsName.add("Nucleus");
            channelsName.add("YH2AX");
            if (channels.size() > 1) {
                chs = tools.dialog(channels, channelsName);
                if ( chs == null) {
                    IJ.showStatus("Plugin cancelled");
                    return;
                }
            }
            
            for (String f : imageFile) {
                rootName = FilenameUtils.getBaseName(f);
                reader.setId(f);
                ImporterOptions options = new ImporterOptions();
                    
                        
                /** 
                 * read nd
                 * Detect nucleus in 405 channel Tomato cells channel3, measure intensity in GFP channel1 and PV channel2
                */

                options.setId(f);
                options.setSplitChannels(true);
                options.setQuiet(true);
                options.setCrop(true);
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);

                       
                // Find nucleus
                System.out.println("Opening nucleus channel " + channels.get(0) +" ...");
                int channel = channels.indexOf(chs.get(0));
                ImagePlus imgNucleus = BF.openImagePlus(options)[channel];
                
                //section volume in mm^3
                double sectionVol = (imgNucleus.getWidth() * cal.pixelWidth * imgNucleus.getHeight() * cal.pixelHeight * imgNucleus.getNSlices() * cal.pixelDepth)/1e9;
                
                // Find nucleus
                Objects3DIntPopulation nucPop =  tools.stardistNucleiPop(imgNucleus);
                
                // open YH2AX Channel
                System.out.println("Opening YH2AX channel " + channels.get(1)+ " ...");
                channel = channels.indexOf(chs.get(1));
                ImagePlus imgYH2AX = BF.openImagePlus(options)[channel];
                ImageHandler imhYH2AX = ImageHandler.wrap(imgYH2AX);
                
                // read nucleus intensity in OFR1P channel 
                ArrayList<Double> nucOFR1P_intensity = tools.readIntensity(imgYH2AX, nucPop);
                
                // Find YH2AX dots
                Objects3DIntPopulation dotsPop = tools.find_dots(imgYH2AX, 2, 1, "Triangle");
                System.out.println(dotsPop.getNbObjects()+" dots found");
                
                // Save image objects
                tools.saveImageObjects(null, dotsPop, nucPop, imgYH2AX, outDirResults+rootName+"_objects.tif", 3, 40);
                
                // write parameters
                for (Object3DInt nucObj : nucPop.getObjects3DInt()) {
                    double nucObjVol = new MeasureVolume(nucObj).getVolumeUnit();
                    int dotsNuc = 0;
                    double dotsNucInt = 0;
                    double dotsNucVol = 0;
                    double nucInt = new MeasureIntensity(nucObj, imhYH2AX).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
                    for (Object3DInt dotObj : dotsPop.getObjects3DInt()) {
                        // find dots inside nucleus
                        Measure2Colocalisation coloc = new Measure2Colocalisation(dotObj, nucObj);
                        if (coloc.getValue(Measure2Colocalisation.COLOC_PC) >= 80) {
                            dotsNuc++;
                            dotsNucVol += new MeasureVolume(dotObj).getVolumeUnit();
                           dotsNucInt += new MeasureIntensity(dotObj, imhYH2AX).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
                        }
                    }
                    nucleus_Analyze.write(rootName+"\t"+(int)nucObj.getLabel()+"\t"+nucObjVol+"\t"+nucInt+"\t"+dotsNuc+"\t"+dotsNucInt+"\t"+dotsNucVol+"\n");
                }
                
                tools.closeImages(imgYH2AX);
                tools.closeImages(imgNucleus);
                
            }
            nucleus_Analyze.close();
        } catch (DependencyException | ServiceException | FormatException | IOException ex) {
           Logger.getLogger(YH2AX.class.getName()).log(Level.SEVERE, null, ex);
       }
       IJ.showStatus("Process done ..."); 
    }    
}
