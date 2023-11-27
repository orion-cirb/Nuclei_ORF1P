import Nuclei_ORF1P_Tools.Tools;
import Nuclei_ORF1P_Tools.Cell;
import ij.IJ;
import ij.ImagePlus;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.HashMap;
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
import loci.plugins.in.ImporterOptions;
import mcib3d.geom2.Objects3DIntPopulation;
import org.apache.commons.io.FilenameUtils;
import org.scijava.util.ArrayUtils;


/*
 * Find ORF1P cells and their corresponding DAPI nuclei
 * Measure nucleus and cytoplasm intensity in ORF1P channel (488)
 */

/**
 * @author phm
 */
public class Nuclei_ORF1P implements PlugIn {
    
    Tools tools = new Tools();
    private final boolean canceled = false;
    private String imageDir = "";
    public  String outDirResults = "";
    public  String rootName = "";
    public BufferedWriter results;
    
    public void run(String arg) {
        try {
            if (canceled) {
                IJ.showMessage("Plugin canceled");
                return;
            }
            if ((!tools.checkInstalledModules()) || (!tools.checkStarDistModels())) {
                return;
            } 

                        
            imageDir = IJ.getDirectory("Choose directory containing image files...");
            if (imageDir == null) {
                return;
            }
            // Find images with nd extension
            ArrayList<String> imageFile = tools.findImages(imageDir, "nd");
            if (imageFile == null) {
                IJ.showMessage("Error", "No images found with nd extension");
                return;
            }
            
            // Create output folder
            outDirResults = imageDir + File.separator + "Results" + File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            // Write header in results file
            String header = "Image name\tImage background\tNucleus ID\tNucleus volume (µm3)\tNucleus compactness\tNucleus sphericity\tNucleus elongation\tNucleus flatness\tNucleus intensity\tNucleus cor. intensity\tNucleus dots number\tNucleus dots volume (µm3)\tNucleus dots intensity\tNucleus dots cor. intensity\t"
                + "Nucleus inner volume (µm3)\tNucleus inner intensity\tNucleus inner cor. intensity\tNucleus inner dots number\tNucleus inner dots volume (µm3)\tNucleus inner dots intensity\tNucleus inner dots cor. intensity\t"
                + "Nucleus inner ring volume (µm3)\tNucleus inner ring intensity\tNucleus inner ring cor. intensity\tNucleus inner ring dots number\tNucleus inner ring dots volume (µm3)\tNucleus inner ring dots intensity\tNucleus inner ring dots cor. intensity\t"
                + "Nucleus outer ring volume (µm3)\tNucleus outer ring intensity\tNucleus outer ring cor. intensity\tNucleus outer ring dots number\tNucleus outer ring dots volume (µm3)\tNucleus outer ring dots intensity\tNucleus outer ring dots cor. intensity\t"
                + "Cytoplasm volume (µm3)\tCytoplasm intensity\tCytoplasm cor. intensity\tCytoplasm dots number\tCytoplasm dots volume (µm3)\tCytoplasm dots intensity\tCytoplasm dots cor. intensity\n";
            FileWriter fwResults = new FileWriter(outDirResults + "results.xls", false);
            results = new BufferedWriter(fwResults);
            results.write(header);
            results.flush();
            
            // Create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            reader.setId(imageFile.get(0));
            
            // Find image calibration
            tools.cal = tools.findImageCalib(meta);
            
            // Find channel names
            String[] channels = tools.findChannels(imageFile.get(0), meta, reader);

            // Channels dialog
            String[] chs = tools.dialog(channels);
            if (chs == null) {
                IJ.showStatus("Plugin cancelled");
                return;
            }
            
            for (String f : imageFile) {
                rootName = FilenameUtils.getBaseName(f);
                tools.print("--- ANALYZING IMAGE " + rootName + " ------");
                reader.setId(f);
                
                ImporterOptions options = new ImporterOptions();
                options.setId(f);
                options.setSplitChannels(true);
                options.setQuiet(true);
                options.setCrop(true);
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);

                // Open DAPI channel
                tools.print("- Analyzing " + tools.channelsName[0] + " channel -");
                int indexCh = ArrayUtils.indexOf(channels, chs[0]);
                ImagePlus imgNucleus = BF.openImagePlus(options)[indexCh];
                
                // Find DAPI nuclei
                System.out.println("Finding " + tools.channelsName[0] + " nuclei....");
                Objects3DIntPopulation nucPop = tools.cellposeDetection(imgNucleus, true, "cyto2", 1, 100, 0.5, 
                        false, tools.minNucVol, tools.maxNucVol);
                System.out.println(nucPop.getNbObjects() + " " + tools.channelsName[0] + " nuclei found");
                //tools.drawPop(imgNucleus, nucPop, outDirResults, "nuclei");
                
                // Open ORF1P channel
                tools.print("- Analyzing " + tools.channelsName[1] + " channel -");
                indexCh = ArrayUtils.indexOf(channels, chs[1]);
                ImagePlus imgORF1P = BF.openImagePlus(options)[indexCh];
                
                // Find background
                double bgORF1P = tools.findBackground(imgORF1P);
                
                // Find cells cytoplasm
                System.out.println("Finding " + tools.channelsName[1] + " cells....");
                Objects3DIntPopulation cellPop = tools.cellposeDetection(imgORF1P, true, "cyto2", 1, 100, 0.5, 
                        true, tools.minCellVol, tools.maxCellVol);
                System.out.println(cellPop.getNbObjects() + " " + tools.channelsName[1] + " cells found");
                //tools.drawPop(imgORF1P, cellPop, outDirResults, "cells");
                 
                // Colocalization
                System.out.println("Finding " + tools.channelsName[1] + " cells colocalizing with a " + tools.channelsName[1] + " nuclei...");
                ArrayList<Cell> colocPop = tools.colocalization(cellPop, nucPop);
                System.out.println(colocPop.size() + " " + tools.channelsName[1] + " cells colocalized with " + tools.channelsName[0] + " nuclei");
                tools.resetLabels(colocPop);
                
                // Find nuclei outer ring
                System.out.println("Finding outer rings...");
                tools.setNucleiRing(colocPop, tools.outerNucDil, true);
                // Find nuclei inner ring and inner nucleus
                System.out.println("Finding inner rings and inner nuclei...");
                tools.setNucleiRing(colocPop, tools.innerNucDil, false);
                
                // Find dots
                Objects3DIntPopulation allDotsPop = new Objects3DIntPopulation();
                if (tools.dotsDetect) {
                    System.out.println("Finding " + tools.channelsName[1] + " dots...");
                    allDotsPop = tools.stardistDetection(imgORF1P, false, tools.stardistModelDots, 
                        tools.stardistProbThreshDots, tools.stardistOverlayThreshDots, tools.minDotVol, tools.maxDotVol);
                    System.out.println(allDotsPop.getNbObjects() + " " + tools.channelsName[1] + " dots found");
                } else {
                    System.out.println("No dots detection performed");
                }
                
                // Tag nuclei with parameters
                tools.print("- Measuring cells parameters -");
                tools.tagCells(imgORF1P, colocPop, allDotsPop);
                
                // Save image objects
                tools.print("- Saving results -");
                tools.drawResults(colocPop, allDotsPop, imgNucleus, imgORF1P, rootName, outDirResults, 40);
                             
                // Write results
                for (Cell cell: colocPop) {
                    HashMap<String, Double> params = cell.params;
                    results.write(rootName+"\t"+bgORF1P+"\t"+(int)((double)params.get("index"))+"\t"+params.get("nucVol")+"\t"+params.get("nucComp")+"\t"+params.get("nucSph")+"\t"+params.get("nucEllElong")+"\t"+params.get("nucEllFlat")+"\t"+params.get("nucInt")+"\t"+(params.get("nucInt")-bgORF1P*params.get("nucVol")) +"\t"+(int)((double)params.get("nucDots"))+"\t"+params.get("nucDotsVol")+"\t"+params.get("nucDotsInt")+"\t"+(params.get("nucDotsInt")-bgORF1P*params.get("nucDotsVol"))+
                            "\t"+params.get("innerNucVol")+"\t"+params.get("innerNucInt")+"\t"+(params.get("innerNucInt")-bgORF1P*params.get("innerNucVol"))+"\t"+(int)((double)params.get("innerNucDots"))+"\t"+params.get("innerNucDotsVol")+"\t"+params.get("innerNucDotsInt")+"\t"+(params.get("innerNucDotsInt")-bgORF1P*params.get("innerNucDotsVol"))+
                            "\t"+params.get("innerRingVol")+"\t"+params.get("innerRingInt")+"\t"+(params.get("innerRingInt")-bgORF1P*params.get("innerRingVol"))+"\t"+(int)((double)params.get("innerRingDots"))+"\t"+params.get("innerRingDotsVol")+"\t"+params.get("innerRingDotsInt")+"\t"+(params.get("innerRingDotsInt")-bgORF1P*params.get("innerRingDotsVol"))+
                            "\t"+params.get("outerRingVol")+"\t"+params.get("outerRingInt")+"\t"+(params.get("outerRingInt")-bgORF1P*params.get("outerRingVol"))+"\t"+(int)((double)params.get("outerRingDots"))+"\t"+params.get("outerRingDotsVol")+"\t"+params.get("outerRingDotsInt")+"\t"+(params.get("outerRingDotsInt")-bgORF1P*params.get("outerRingDotsVol"))+"\t"+
                            params.get("cytoVol")+"\t"+params.get("cytoInt")+"\t"+(params.get("cytoInt")-bgORF1P*params.get("cytoVol"))+"\t"+(int)((double)params.get("cytoDots"))+"\t"+params.get("cytoDotsVol")+"\t"+params.get("cytoDotsInt")+"\t"+(params.get("cytoDotsInt")-bgORF1P*params.get("cytoDotsVol"))+"\n");
                    results.flush();
                }
                  
                tools.flush_close(imgNucleus);
                tools.flush_close(imgORF1P);
            }
            
            results.close();
        } catch (IOException | DependencyException | ServiceException | FormatException  ex) {
            Logger.getLogger(Nuclei_ORF1P.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        tools.print("--- All done! ---");
    }
}
