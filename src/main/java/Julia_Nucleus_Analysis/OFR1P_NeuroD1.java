package Julia_Nucleus_Analysis;



import Tools.Jnucleus_Tools3D;
import Tools.Nucleus;
import ij.IJ;
import ij.ImagePlus;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.util.ImageProcessorReader;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.plugin.PlugIn;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;
import loci.plugins.in.ImporterOptions;
import mcib3d.geom.Objects3DPopulation;
import org.apache.commons.io.FilenameUtils;


/*
 * Find Nucleus intensity in OFR1P channel (488)
 * and cytoplasmic intensity in NeuroD1 channel (561)
 * 
 */

/**
 *
 * @author phm
 */
public class OFR1P_NeuroD1 implements PlugIn {
    
    Jnucleus_Tools3D tools = new Jnucleus_Tools3D();
    
    private final boolean canceled = false;
    private String imageDir = "";
    public  String outDirResults = "";
    public  String rootName = "";
    public BufferedWriter nucleus_Analyze;
    public BufferedWriter nucleusGlobal_Analyze;
    
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
            tools.cal = tools.findImageCalib(meta);


            // write headers
            String header= "Image Name\t# Nucleus\tNucleus Volume (µm3)\tNucleus intensity\tNucleus corrected intensity\tNucleus dots number\tNucleus dots volume (µm3)\tNucleus dots intensity\t"
                    + "Nucleus inner volume (µm3)" + "\t" + "Nucleus inner intensity\tNucleus inner corrected intensity\tNucleus inner dots number\tNucleus inner dots volume (µm3)\tNucleus inner dots intensity\t"
                    + "Nucleus inner ring volume (µm3)" + "\t" + "Nucleus inner ring intensity\tNucleus inner ring corrected intensity\tNucleus inner ring dots number\tNucleus inner ring dots volume (µm3)\tNucleus inner ring dots intensity\t"
                    + "Nucleus outer ring volume (µm3)" + "\t" + "Nucleus outer ring intensity\tNucleus outer ring corrected intensity\tNucleus outer ring dots number\tNucleus outer ring dots volume (µm3)\tNucleus outer ring dots intensity\t"
                    + "Cytoplasm Vol (µm3)\tCytoplasm intensity\tCytoplasm corrected Cell intensity\tCytoplasm dots number\tCytoplasm dots volume (µm3)\tCytoplasm dots intensity\n";
            FileWriter fwNucleus = new FileWriter(outDirResults + "nucleus_Results.xls", false);
            nucleus_Analyze = new BufferedWriter(fwNucleus);
            nucleus_Analyze.write(header);
            nucleus_Analyze.flush();
            header= "Image Name\tNucleus number\tSection volume (µm3)\tBackground Intensity\tStd BG\tCorrected sum intensity of nucleus masked\tCorrected of sum intensity of cells processes\n";
            FileWriter fwNucleusGlobal = new FileWriter(outDirResults + "nucleus_GlobalResults.xls", false);
            nucleusGlobal_Analyze = new BufferedWriter(fwNucleusGlobal);
            nucleusGlobal_Analyze.write(header);
            nucleusGlobal_Analyze.flush();
            
            // Channels dialog
            List<String> chs = new ArrayList();
            List<String> channelsName = new ArrayList();
            channelsName.add("Nucleus");
            channelsName.add("OFR1P");
            channelsName.add("NeuroD1");
            if (channels.size() > 1) {
                chs = tools.dialog(channels, channelsName, true);
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
               
                // Find Z where nucleus stack intensity max
                ImagePlus imgZCropNuc = (tools.zCrop) ? tools.cropZmax(imgNucleus) : new Duplicator().run(imgNucleus);
                tools.closeImages(imgNucleus);
                
                //section volume in µm^3
                double volPix = tools.cal.pixelWidth*tools.cal.pixelHeight*tools.cal.pixelDepth;
                double sectionVol = (imgZCropNuc.getWidth() * imgZCropNuc.getHeight() * imgZCropNuc.getNSlices() * volPix);
                
                // Find nucleus
                Objects3DPopulation nucPop = new Objects3DPopulation();
                if (tools.nucleusDetector.equals("StarDist"))
                    nucPop = tools.stardistNucleiPop(imgZCropNuc);
                else
                    nucPop = tools.findNucleus(imgZCropNuc, 18, 20, 2, "Triangle");
                System.out.println(nucPop.getNbObjects()+" nucleus founds");
                
                // Find inner/outer ring nucleus
                // outer
                System.out.println("Finding outer ring ....");
                Objects3DPopulation outerRingPop = tools.createDonutPop(nucPop, imgZCropNuc, tools.outerNucDil, true);
                // inner
                System.out.println("Finding inner ring ....");
                Objects3DPopulation innerRingPop = tools.createDonutPop(nucPop, imgZCropNuc, tools.innerNucDil, false);
                // inner nucleus
                System.out.println("Finding inner nucleus ....");
                Objects3DPopulation innerNucPop = tools.getInnerNucleus(nucPop, imgZCropNuc, tools.innerNucDil, false);
                
                tools.closeImages(imgZCropNuc);
                   
                // open OFRP1 Channel
                System.out.println("Opening OFR1P channel " + channels.get(1)+ " ...");
                channel = channels.indexOf(chs.get(1));
                ImagePlus imgOFR1P = BF.openImagePlus(options)[channel];
                
                // Take same stack as nucleus
                ImagePlus imgZCropOFR1P = new Duplicator().run(imgOFR1P, tools.zMax, imgOFR1P.getNSlices());
                tools.closeImages(imgOFR1P);
                
                // Find background
                double[] bgOFR1P = tools.find_background(imgZCropOFR1P);
                
                // Find cell cytoplasm
                Objects3DPopulation cellPop = tools.findCells(imgZCropOFR1P, nucPop);
                
                // mask OFR1P with nucleus object
                ImagePlus imgOFR1P_NucleusMask = tools.maskImage(imgZCropOFR1P, nucPop, outDirResults, rootName+"_nucleusMasked.tif");
                
                 // read all intensity in OFR1P nucleus masked channel
                double sumNucMaskedIntensity = tools.readSumIntensity(imgOFR1P_NucleusMask, 0, null);
                
                // mask OFR1P with nucleus and cell object
                ImagePlus imgOFR1P_CellMask = tools.maskImage(imgOFR1P_NucleusMask, cellPop, outDirResults, rootName+"_CellsMasked.tif");
                
                tools.closeImages(imgOFR1P_NucleusMask);
                // read all intensity in OFR1P nucleus cells masked channel
                double sumCellMaskedIntensity = tools.readSumIntensity(imgOFR1P_CellMask, bgOFR1P[0]+bgOFR1P[1], outDirResults+rootName+"_CellProcesses.tif");
                tools.closeImages(imgOFR1P_CellMask);
                
                // Dots detections
                // All dots population
                Objects3DPopulation allDotsPop = new Objects3DPopulation();
                if (tools.dotsDetect)
                    allDotsPop = tools.find_dots(imgZCropOFR1P, 2, 1, "Triangle");

                // Save image objects
                tools.saveImageObjects(nucPop, outerRingPop, null, imgZCropOFR1P, outDirResults+rootName+"_OuterRingObjects.tif", 40);
                tools.saveImageObjects(innerRingPop, innerNucPop, null, imgZCropOFR1P, outDirResults+rootName+"_innerRingObjects.tif", 40);
                if (tools.dotsDetect)
                    tools.saveImageObjects(nucPop, allDotsPop, cellPop, imgZCropOFR1P, outDirResults+rootName+"_dotsObjects.tif", 40);
                
                // tags nucleus with parameters
                ArrayList<Nucleus> nucleus = tools.tagsNuclei(imgZCropOFR1P, nucPop, innerNucPop, innerRingPop, outerRingPop, cellPop, allDotsPop);
                             
                // Write results
                for (Nucleus nuc : nucleus) {
                    nucleus_Analyze.write(rootName+"\t"+nuc.getIndex()+"\t"+nuc.getNucVol()+"\t"+nuc.getNucInt()+"\t"+(nuc.getNucInt() - bgOFR1P[0] * (nuc.getNucVol()/volPix)) +"\t"+nuc.getNucDots()+"\t"+nuc.getNucDotsVol()+"\t"+nuc.getNucDotsInt()+
                            "\t"+nuc.getInnerNucVol()+"\t"+nuc.getInnerNucInt()+"\t"+(nuc.getInnerNucInt() - bgOFR1P[0] * (nuc.getInnerNucVol()/volPix))+"\t"+nuc.getInnerNucDots()+"\t"+nuc.getInnerNucDotsVol()+"\t"+nuc.getInnerNucDotsInt()+
                            "\t"+nuc.getInnerRingVol()+"\t"+nuc.getInnerRingInt()+"\t"+(nuc.getInnerRingInt() - bgOFR1P[0] * (nuc.getInnerRingVol()/volPix))+"\t"+nuc.getInnerRingDots()+"\t"+nuc.getInnerRingDotsVol()+"\t"+nuc.getInnerRingDotsInt()+
                            "\t"+nuc.getOuterRingVol()+"\t"+nuc.getOuterRingInt()+"\t"+(nuc.getOuterRingInt() - bgOFR1P[0] * (nuc.getOuterRingVol()/volPix))+"\t"+nuc.getOuterRingDots()+"\t"+nuc.getOuterRingDotsVol()+"\t"+nuc.getOuterRingDotsInt()+
                            "\t"+nuc.getCytoVol()+"\t"+nuc.getCytoInt()+"\t"+(nuc.getCytoInt() - bgOFR1P[0] * (nuc.getCytoVol()/volPix))+"\t"+nuc.getCytoDots()+"\t"+nuc.getCytoDotsVol()+"\t"+nuc.getCytoDotsInt()+"\n");
                    nucleus_Analyze.flush();
                }
                // Global measurements
                
                nucleusGlobal_Analyze.write(rootName+"\t"+nucPop.getNbObjects()+"\t"+sectionVol+"\t"+bgOFR1P[0]+"\t"+bgOFR1P[1]+"\t"+sumNucMaskedIntensity+"\t"+sumCellMaskedIntensity+"\n");
                nucleusGlobal_Analyze.flush();               
                tools.closeImages(imgZCropOFR1P);
            }
            nucleus_Analyze.close();
            nucleusGlobal_Analyze.close();
        } catch (IOException | DependencyException | ServiceException | FormatException  ex) {
            Logger.getLogger(OFR1P_NeuroD1.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
        IJ.showStatus("Process done ...");
    }
}
