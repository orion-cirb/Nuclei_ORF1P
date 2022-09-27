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
import ij.plugin.PlugIn;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;
import loci.plugins.in.ImporterOptions;
import mcib3d.geom2.Objects3DIntPopulation;
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
            String header= "Image Name\t# Nucleus\tNucleus Volume (µm3)\tNucleus Sphericity\tNucleus Compactness\tNucleus elongation\tNucleus flatness\tNucleus intensity\tNucleus cor. intensity\tNucleus dots number\tNucleus dots volume (µm3)\tNucleus dots intensity\tNucleus dots cor. intensity\t"
                + "Nucleus inner volume (µm3)\tNucleus inner intensity\tNucleus inner cor. intensity\tNucleus inner dots number\tNucleus inner dots volume (µm3)\tNucleus inner dots intensity\tNucleus inner dots cor. intensity\t"
                + "Nucleus inner ring volume (µm3)\tNucleus inner ring intensity\tNucleus inner ring cor. intensity\tNucleus inner ring dots number\tNucleus inner ring dots volume (µm3)\tNucleus inner ring dots intensity\tNucleus inner ring dots cor. intensity\t"
                + "Nucleus outer ring volume (µm3)\tNucleus outer ring intensity\tNucleus outer ring cor. intensity\tNucleus outer ring dots number\tNucleus outer ring dots volume (µm3)\tNucleus outer ring dots intensity\tNucleus outer ring dots cor. intensity\t"
                + "Cytoplasm Vol (µm3)\tCytoplasm intensity\tCytoplasm cor. Cell intensity\tCytoplasm dots number\tCytoplasm dots volume (µm3)\tCytoplasm dots intensity\tCytoplasm dots cor. intensity\n";
            FileWriter fwNucleus = new FileWriter(outDirResults + "nucleus_Results.xls", false);
            nucleus_Analyze = new BufferedWriter(fwNucleus);
            nucleus_Analyze.write(header);
            nucleus_Analyze.flush();
            
            // Channels dialog
            List<String> chs = new ArrayList();
            List<String> channelsName = new ArrayList();
            channelsName.add("Nucleus");
            channelsName.add("ORF1P");
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
                
                //section volume in µm^3
                double volPix = tools.cal.pixelWidth*tools.cal.pixelHeight*tools.cal.pixelDepth;
                double sectionVol = (imgNucleus.getWidth() * imgNucleus.getHeight() * imgNucleus.getNSlices() * volPix);
                
                // Find nucleus
                Objects3DIntPopulation nucPop =  tools.stardistObjectsPop(imgNucleus, "nucleus");
                System.out.println(nucPop.getNbObjects()+" nucleus founds");

                // Find inner/outer ring nucleus
                // outer
                System.out.println("Finding outer ring ....");
                Objects3DIntPopulation outerRingPop = tools.createDonutPop(nucPop, imgNucleus, tools.outerNucDil, true);
                // inner
                System.out.println("Finding inner ring ....");
                Objects3DIntPopulation innerRingPop = tools.createDonutPop(nucPop, imgNucleus, tools.innerNucDil, false);
                // inner nucleus
                System.out.println("Finding inner nucleus ....");
                Objects3DIntPopulation innerNucPop = tools.getInnerNucleus(nucPop, imgNucleus);
                
                // Save image objects
                tools.saveImageObjects(null, outerRingPop, nucPop, imgNucleus, outDirResults+rootName+"_OuterRingObjects.tif", 0, 40);
                tools.closeImages(imgNucleus);
                
                // open OFRP1 Channel
                System.out.println("Opening OFR1P channel " + channels.get(1)+ " ...");
                channel = channels.indexOf(chs.get(1));
                ImagePlus imgOFR1P = BF.openImagePlus(options)[channel];
                
                // Find background
                double bgOFR1P = tools.find_background(imgOFR1P);
                
                // Find cells cytoplasm
                Objects3DIntPopulation cellPop = new Objects3DIntPopulation();
                cellPop = tools.cellPoseCellsPop(imgOFR1P, nucPop);
                System.out.println(cellPop.getNbObjects()+" cells found with nucleus");
                
                // Save cells cytoplasm
                tools.saveImageObjects(null, cellPop, nucPop, imgOFR1P, outDirResults+rootName+"_CellsCytoplasmObjects.tif", 2, 40);

                // Dots detections
                // All dots population
                Objects3DIntPopulation allDotsPop = new Objects3DIntPopulation();
                if (tools.dotsDetect)
                    allDotsPop = tools.stardistObjectsPop(imgOFR1P, "dots");
                System.out.println(allDotsPop.getNbObjects()+" dots found");

                // Save image objects
                tools.saveImageObjects(null, outerRingPop, nucPop, imgOFR1P, outDirResults+rootName+"_OuterRingObjects.tif", 0, 40);
                tools.saveImageObjects(innerNucPop, innerRingPop, null, imgOFR1P, outDirResults+rootName+"_innerRingObjects.tif", 0, 40);
                if (tools.dotsDetect)
                    tools.saveImageObjects(allDotsPop, cellPop, nucPop, imgOFR1P, outDirResults+rootName+"_dotsObjects.tif", 3, 40);
                
                // tags nucleus with parameters
                ArrayList<Nucleus> nucleus = tools.tagsNuclei(imgOFR1P, nucPop, innerNucPop, innerRingPop, outerRingPop, cellPop, allDotsPop);
                             
                // Write resultsSearch (Ctrl+I)
                for (Nucleus nuc : nucleus) {
                    nucleus_Analyze.write(rootName+"\t"+nuc.getIndex()+"\t"+nuc.getNucVol()+"\t"+nuc.getNucComp()+"\t"+nuc.getNucSph()+"\t"+nuc.getNucEllElong()+"\t"+
                        nuc.getNucEllFlat()+"\t"+nuc.getNucInt()+"\t"+(nuc.getNucInt() - bgOFR1P*nuc.getNucVol()) +"\t"+nuc.getNucDots()+"\t"+nuc.getNucDotsVol()+"\t"+
                        nuc.getNucDotsInt()+"\t"+(nuc.getNucDotsInt()- bgOFR1P*nuc.getNucDotsVol())+"\t"+nuc.getInnerNucVol()+"\t"+nuc.getInnerNucInt()+"\t"+(nuc.getInnerNucInt() - bgOFR1P * nuc.getInnerNucVol())+"\t"+
                        nuc.getInnerNucDots()+"\t"+nuc.getInnerNucDotsVol()+"\t"+nuc.getInnerNucDotsInt()+"\t"+(nuc.getInnerNucDotsInt() - bgOFR1P*nuc.getInnerNucDotsVol())+"\t"+nuc.getInnerRingVol()+"\t"+nuc.getInnerRingInt()+"\t"+
                        (nuc.getInnerRingInt() - bgOFR1P * nuc.getInnerRingVol())+"\t"+nuc.getInnerRingDots()+"\t"+nuc.getInnerRingDotsVol()+"\t"+nuc.getInnerRingDotsInt()+"\t"+(nuc.getInnerRingDotsInt()- bgOFR1P * nuc.getInnerRingDotsVol())+"\t"+
                        nuc.getOuterRingVol()+"\t"+nuc.getOuterRingInt()+"\t"+(nuc.getOuterRingInt() - bgOFR1P * nuc.getOuterRingVol())+"\t"+nuc.getOuterRingDots()+"\t"+
                        nuc.getOuterRingDotsVol()+"\t"+nuc.getOuterRingDotsInt()+"\t"+(nuc.getOuterRingDotsInt() - bgOFR1P*nuc.getOuterRingDotsVol())+"\t"+nuc.getCytoVol()+"\t"+nuc.getCytoInt()+"\t"+
                        (nuc.getCytoInt() - bgOFR1P * nuc.getCytoVol())+"\t"+nuc.getCytoDots()+"\t"+nuc.getCytoDotsVol()+"\t"+nuc.getCytoDotsInt()+"\t"+(nuc.getCytoDotsInt()-bgOFR1P*nuc.getCytoDotsVol())+"\n");
                    nucleus_Analyze.flush();
                }
                        
                tools.closeImages(imgOFR1P);
            }
            nucleus_Analyze.close();
        } catch (IOException | DependencyException | ServiceException | FormatException  ex) {
            Logger.getLogger(OFR1P_NeuroD1.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
        IJ.showStatus("Process done ...");
    }
}
