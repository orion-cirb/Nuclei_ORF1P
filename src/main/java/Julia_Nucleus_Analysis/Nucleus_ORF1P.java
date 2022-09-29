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
public class Nucleus_ORF1P implements PlugIn {
    
    Jnucleus_Tools3D tools = new Jnucleus_Tools3D();
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
            String header = "Image Name\t# Nucleus\tNucleus Volume (µm3)\tNucleus Sphericity\tNucleus Compactness\tNucleus elongation\tNucleus flatness\tNucleus intensity\tNucleus cor. intensity\tNucleus dots number\tNucleus dots volume (µm3)\tNucleus dots intensity\tNucleus dots cor. intensity\t"
                + "Nucleus inner volume (µm3)\tNucleus inner intensity\tNucleus inner cor. intensity\tNucleus inner dots number\tNucleus inner dots volume (µm3)\tNucleus inner dots intensity\tNucleus inner dots cor. intensity\t"
                + "Nucleus inner ring volume (µm3)\tNucleus inner ring intensity\tNucleus inner ring cor. intensity\tNucleus inner ring dots number\tNucleus inner ring dots volume (µm3)\tNucleus inner ring dots intensity\tNucleus inner ring dots cor. intensity\t"
                + "Nucleus outer ring volume (µm3)\tNucleus outer ring intensity\tNucleus outer ring cor. intensity\tNucleus outer ring dots number\tNucleus outer ring dots volume (µm3)\tNucleus outer ring dots intensity\tNucleus outer ring dots cor. intensity\t"
                + "Cytoplasm Vol (µm3)\tCytoplasm intensity\tCytoplasm cor. Cell intensity\tCytoplasm dots number\tCytoplasm dots volume (µm3)\tCytoplasm dots intensity\tCytoplasm dots cor. intensity\n";
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
                tools.print("- Analyzing " + channels[0] + " channel -");
                int indexCh = ArrayUtils.indexOf(channels, chs[0]);
                ImagePlus imgNucleus = BF.openImagePlus(options)[indexCh];
                
                // Find DAPI nuclei
                Objects3DIntPopulation nucPop = tools.cellposeDetection(imgNucleus, true, "cyto2", 1, 100, 0.5, 
                        false, tools.minNucVol, tools.maxNucVol);
                System.out.println(nucPop.getNbObjects() + " " + channels[0] + " nuclei found");
                tools.drawPop(imgNucleus, nucPop, outDirResults, "NUCLEI");

                // Find nucleus outer ring
                System.out.println("Finding outer ring ....");
                Objects3DIntPopulation outerRingPop = tools.createDonutPop(nucPop, imgNucleus, tools.outerNucDil, true);
                // Find nucleus inner ring
                System.out.println("Finding inner ring ....");
                Objects3DIntPopulation innerRingPop = tools.createDonutPop(nucPop, imgNucleus, tools.innerNucDil, false);
                // Find inner nucleus
                System.out.println("Finding inner nucleus ....");
                Objects3DIntPopulation innerNucPop = tools.getInnerNucleus(nucPop, imgNucleus);
                
                // Open ORF1P channel
                tools.print("- Analyzing " + channels[1] + " channel -");
                indexCh = ArrayUtils.indexOf(channels, chs[1]);
                ImagePlus imgORF1P = BF.openImagePlus(options)[indexCh];
                
                // Find background
                double bgORF1P = tools.findBackground(imgORF1P);
                
                // Find cells cytoplasm
                Objects3DIntPopulation cellPop = tools.cellposeDetection(imgORF1P, true, "cyto2", 1, 100, 0.5, 
                        true, tools.minCellVol, tools.maxCellVol);
                cellPop = tools.colocalization(cellPop, nucPop);
                System.out.println(cellPop.getNbObjects() + " " + channels[1] + " cells found with a nucleus");
                tools.drawPop(imgORF1P, cellPop, outDirResults, "CELLS");
                
                // Find dots
                Objects3DIntPopulation allDotsPop = new Objects3DIntPopulation();
                if (tools.dotsDetect) {
                    allDotsPop = tools.stardistDetection(imgORF1P, false, tools.stardistModelDots, 
                        tools.stardistProbThreshDots, tools.stardistOverlayThreshDots, tools.minDotVol, tools.maxDotVol);
                    System.out.println(allDotsPop.getNbObjects() + " " + channels[1] + " dots found");
                } else {
                    System.out.println("No dots detection performed");
                }
                
                // Tag nuclei with parameters
                tools.print("- Measuring cells parameters -");
                ArrayList<Nucleus> nuclei = tools.tagNuclei(imgORF1P, nucPop, innerNucPop, innerRingPop, outerRingPop, cellPop, allDotsPop);
                
                // Save image objects
                tools.print("- Saving results -");
                tools.saveImageObjects(null, cellPop, nucPop, imgORF1P, outDirResults+rootName+"_CellsCytoplasmObjects.tif", 2, 40);
                tools.saveImageObjects(null, outerRingPop, nucPop, imgORF1P, outDirResults+rootName+"_OuterRingObjects.tif", 0, 40);
                tools.saveImageObjects(innerNucPop, innerRingPop, null, imgORF1P, outDirResults+rootName+"_innerRingObjects.tif", 0, 40);
                if (tools.dotsDetect)
                    tools.saveImageObjects(allDotsPop, cellPop, nucPop, imgORF1P, outDirResults+rootName+"_dotsObjects.tif", 3, 40);
                             
                // Write results
                for (Nucleus nuc : nuclei) {
                    results.write(rootName+"\t"+nuc.getIndex()+"\t"+nuc.getNucVol()+"\t"+nuc.getNucComp()+"\t"+nuc.getNucSph()+"\t"+nuc.getNucEllElong()+"\t"+
                        nuc.getNucEllFlat()+"\t"+nuc.getNucInt()+"\t"+(nuc.getNucInt() - bgORF1P*nuc.getNucVol()) +"\t"+nuc.getNucDots()+"\t"+nuc.getNucDotsVol()+"\t"+
                        nuc.getNucDotsInt()+"\t"+(nuc.getNucDotsInt()- bgORF1P*nuc.getNucDotsVol())+"\t"+nuc.getInnerNucVol()+"\t"+nuc.getInnerNucInt()+"\t"+(nuc.getInnerNucInt() - bgORF1P * nuc.getInnerNucVol())+"\t"+
                        nuc.getInnerNucDots()+"\t"+nuc.getInnerNucDotsVol()+"\t"+nuc.getInnerNucDotsInt()+"\t"+(nuc.getInnerNucDotsInt() - bgORF1P*nuc.getInnerNucDotsVol())+"\t"+nuc.getInnerRingVol()+"\t"+nuc.getInnerRingInt()+"\t"+
                        (nuc.getInnerRingInt() - bgORF1P * nuc.getInnerRingVol())+"\t"+nuc.getInnerRingDots()+"\t"+nuc.getInnerRingDotsVol()+"\t"+nuc.getInnerRingDotsInt()+"\t"+(nuc.getInnerRingDotsInt()- bgORF1P * nuc.getInnerRingDotsVol())+"\t"+
                        nuc.getOuterRingVol()+"\t"+nuc.getOuterRingInt()+"\t"+(nuc.getOuterRingInt() - bgORF1P * nuc.getOuterRingVol())+"\t"+nuc.getOuterRingDots()+"\t"+
                        nuc.getOuterRingDotsVol()+"\t"+nuc.getOuterRingDotsInt()+"\t"+(nuc.getOuterRingDotsInt() - bgORF1P*nuc.getOuterRingDotsVol())+"\t"+nuc.getCytoVol()+"\t"+nuc.getCytoInt()+"\t"+
                        (nuc.getCytoInt() - bgORF1P * nuc.getCytoVol())+"\t"+nuc.getCytoDots()+"\t"+nuc.getCytoDotsVol()+"\t"+nuc.getCytoDotsInt()+"\t"+(nuc.getCytoDotsInt()-bgORF1P*nuc.getCytoDotsVol())+"\n");
                    results.flush();
                }
                
                tools.flush_close(imgNucleus);
                tools.flush_close(imgORF1P);
            }
            
            results.close();
        } catch (IOException | DependencyException | ServiceException | FormatException  ex) {
            Logger.getLogger(Nucleus_ORF1P.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        tools.print("--- All done! ---");
    }
}
