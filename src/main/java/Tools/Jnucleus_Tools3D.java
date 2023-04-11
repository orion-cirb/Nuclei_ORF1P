package Tools;

import Cellpose.CellposeSegmentImgPlusAdvanced;
import Cellpose.CellposeTaskSettings;
import StardistOrion.StarDist2D;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.plugin.RGBStackMerge;
import ij.plugin.ZProjector;
import ij.process.ImageProcessor;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.BoundingBox;
import mcib3d.geom2.Object3DComputation;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.Objects3DIntPopulationComputation;
import mcib3d.geom2.measurements.MeasureCompactness;
import mcib3d.geom2.measurements.MeasureEllipsoid;
import mcib3d.geom2.measurements.MeasureIntensity;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.geom2.measurementsPopulation.MeasurePopulationColocalisation;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import org.apache.commons.io.FilenameUtils;
import org.scijava.util.ArrayUtils;



/**
 * @author phm
 */
public class Jnucleus_Tools3D {
    public final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
    
    public String[] channelsName = {"DAPI", "ORF1P"};  
    public Calibration cal = new Calibration();
    public float pixVol = 0;
    public double minNucVol= 50;
    public double maxNucVol = 1000;
    public float innerNucDil = 1;
    public float outerNucDil = 1;
    public double minCellVol= 100;
    public double maxCellVol = 2000;
    public boolean dotsDetect = true;
    public double minDotVol= 0.1;
    public double maxDotVol = 500;
    
    // StarDist
    public File modelsPath = new File(IJ.getDirectory("imagej")+File.separator+"models");
    public String stardistModelDots = "pmls2.zip";
    public Object syncObject = new Object();
    public final double stardistPercentileBottom = 0.2;
    public final double stardistPercentileTop = 99.8;
    public final double stardistProbThreshNuc = 0.5;
    public final double stardistOverlayThreshNuc = 0.0;
    public double stardistProbThreshDots = 0.2;
    public final double stardistOverlayThreshDots = 0.25;
    public String stardistOutput = "Label Image"; 
    
    // Cellpose
    public int cellPoseDiameter = 100;
    public String cellPoseModel = "cyto2";
    public String cellPoseEnvDirPath = (IJ.isWindows()) ? System.getProperty("user.home")+"\\miniconda3\\envs\\CellPose" : "/opt/miniconda3/envs/cellpose";
      
 
    /**
     * Display a message in the ImageJ console and status bar
     */
    public void print(String log) {
        System.out.println(log);
        IJ.showStatus(log);
    }
    
    
    /**
     * Check that needed modules are installed
     */
    public boolean checkInstalledModules() {
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("mcib3d.geom.Object3D");
        } catch (ClassNotFoundException e) {
            IJ.showMessage("Error", "3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        return true;
    }
    
    
    /**
     * Check that required StarDist models are present in Fiji models folder
     */
    public boolean checkStarDistModels() {
        FilenameFilter filter = (dir, name) -> name.endsWith(".zip");
        File[] modelList = modelsPath.listFiles(filter);
        int index = ArrayUtils.indexOf(modelList, new File(modelsPath+File.separator+stardistModelDots));
        if (index == -1) {
            IJ.showMessage("Error", stardistModelDots + " StarDist model not found, please add it in Fiji models folder");
            return false;
        }
        return true;
    }

    
    /**
     * Find images in folder
     */
    public ArrayList findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No image found in " + imagesFolder);
            return null;
        }
        ArrayList<String> images = new ArrayList();
        for (String f : files) {
            // Find images with extension
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals(imageExt) && !f.startsWith("."))
                images.add(imagesFolder + File.separator + f);
        }
        Collections.sort(images);
        return(images);
    }
    
   
    /**
     * Find image calibration
     */
    public Calibration findImageCalib(IMetadata meta) {
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        System.out.println("XY calibration = " + cal.pixelWidth + ", Z calibration = " + cal.pixelDepth);
        return(cal);
    }
    
    
    /**
     * Find channels name
     * @throws loci.common.services.DependencyException
     * @throws loci.common.services.ServiceException
     * @throws loci.formats.FormatException
     * @throws java.io.IOException
     */
    public String[] findChannels (String imageName, IMetadata meta, ImageProcessorReader reader) throws DependencyException, ServiceException, FormatException, IOException {
        int chs = reader.getSizeC();
        String[] channels = new String[chs];
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                for (int n = 0; n < chs; n++) 
                {
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n).toString();
                }
                break;
            case "lif" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null || meta.getChannelName(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n).toString();
                break;
            case "czi" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelFluor(0, n).toString();
                break;
            case "ics" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelExcitationWavelength(0, n).value().toString();
                break;    
            default :
                for (int n = 0; n < chs; n++)
                    channels[n] = Integer.toString(n);
        }
        return(channels);         
    }
    
    
    /**
     * Generate dialog box
     */
    public String[] dialog(String[] channels) {     
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 80, 0);
        gd.addImage(icon);
          
        gd.addMessage("Channels", new Font(Font.MONOSPACED , Font.BOLD, 12), Color.blue);
        int index = 0;
        for (String chName: channelsName) {
            gd.addChoice(chName + ": ", channels, channels[index]);
            index++;
        }
        
        gd.addMessage("Nuclei detection", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min nucleus volume (µm3):", minNucVol);
        gd.addNumericField("Max nucleus volume (µm3):", maxNucVol);   
        
        gd.addMessage("Nuclei donut", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Nucleus outer ring (µm):", outerNucDil);
        gd.addNumericField("Nucleus inner ring (µm):", innerNucDil);
        
        gd.addMessage("Cells detection", Font.getFont("Monospace"), Color.blue);
        gd.addDirectoryField("Cellpose environment directory: ", cellPoseEnvDirPath);
        gd.addNumericField("Min cell volume (µm3): ", minCellVol);
        gd.addNumericField("Max cell volume (µm3): ", maxCellVol);
        
        gd.addMessage("Dots detection", Font.getFont("Monospace"), Color.blue);
        gd.addCheckbox("Detect dots", dotsDetect);
        gd.addNumericField("Stardist probability threshold (decrease = more dots): ", stardistProbThreshDots);
        gd.addNumericField("Min dot volume (µm3):", minDotVol);
        gd.addNumericField("Max dot volume (µm3):", maxDotVol);  
        
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("XY calibration (µm):", cal.pixelWidth);
        gd.addNumericField("Z calibration (µm):", cal.pixelDepth);
        gd.showDialog();
        
        String[] ch = new String[channelsName.length];
        for (int i = 0; i < channelsName.length; i++)
            ch[i] = gd.getNextChoice();
        if(gd.wasCanceled())
            ch = null;
       
        minNucVol = (float) gd.getNextNumber();
        maxNucVol = (float) gd.getNextNumber();
        outerNucDil = (float) gd.getNextNumber();
        innerNucDil = (float) gd.getNextNumber();
        cellPoseEnvDirPath = gd.getNextString();
        minCellVol= (float) gd.getNextNumber();
        maxCellVol = (float) gd.getNextNumber();
        dotsDetect = gd.getNextBoolean();
        stardistProbThreshDots = gd.getNextNumber();
        minDotVol= (float) gd.getNextNumber();
        maxDotVol = (float) gd.getNextNumber();
        cal.pixelWidth = gd.getNextNumber();
        cal.pixelDepth = gd.getNextNumber();
        cal.pixelHeight = cal.pixelWidth;
        pixVol = (float) (cal.pixelWidth*cal.pixelHeight*cal.pixelDepth);

        return(ch);
    }
    
    
    /**
     * Flush and close an image
     */
    public void flush_close(ImagePlus img) {
        img.flush();
        img.close();
    }
    
    
    /**
     * Look for all 3D cells in a Z-stack: 
     * - apply CellPose in 2D slice by slice 
     * - let CellPose reconstruct cells in 3D using the stitch threshold parameters
     */
    public Objects3DIntPopulation cellposeDetection(ImagePlus img, boolean resize, String cellposeModel, int channel, int diameter, double stitchThreshold, boolean zFilter, double volMin, double volMax) throws IOException{
        float resizeFactor;
        ImagePlus imgResized;
        if (resize && (img.getWidth() > 1024)) {
            resizeFactor = 0.5f;
            imgResized = img.resize((int)(img.getWidth()*resizeFactor), (int)(img.getHeight()*resizeFactor), 1, "none");
        } else {
            resizeFactor = 1;
            imgResized = new Duplicator().run(img);
            resize = false;
        }

        // Define CellPose settings
        CellposeTaskSettings settings = new CellposeTaskSettings(cellposeModel, channel, (int)(diameter*resizeFactor), cellPoseEnvDirPath);
        settings.setStitchThreshold(stitchThreshold);
        settings.useGpu(true);
       
        // Run CellPose
        CellposeSegmentImgPlusAdvanced cellpose = new CellposeSegmentImgPlusAdvanced(settings, imgResized);
        ImagePlus imgOut = cellpose.run();
        if(resize) imgOut = imgOut.resize(img.getWidth(), img.getHeight(), "none");
        imgOut.setCalibration(cal);
        //imgOut.show();
        //new WaitForUserDialog("test").show();
       
        // Get cells as a population of objects
        ImageHandler imgH = ImageHandler.wrap(imgOut);
        Objects3DIntPopulation pop = new Objects3DIntPopulation(imgH);
        if (zFilter)
            pop = zFilterPop(pop);
        System.out.println(pop.getNbObjects() + " CellPose detections");
       
        // Filter cells by size
        Objects3DIntPopulationComputation popComputation = new Objects3DIntPopulationComputation​(pop);
        Objects3DIntPopulation popFilter = popComputation.getFilterSize​(volMin/pixVol, volMax/pixVol);
        popFilter.resetLabels();
        System.out.println(popFilter.getNbObjects() + " detections remaining after size filtering (" + (pop.getNbObjects()-popFilter.getNbObjects()) + " filtered out)");
       
        flush_close(imgOut);
        imgH.closeImagePlus();
        return(popFilter);
    } 
    
    
    /*
     * Remove objects present in only one z slice from population 
     */
    public Objects3DIntPopulation zFilterPop (Objects3DIntPopulation pop) {
        Objects3DIntPopulation popZ = new Objects3DIntPopulation();
        for (Object3DInt obj : pop.getObjects3DInt()) {
            int zmin = obj.getBoundingBox().zmin;
            int zmax = obj.getBoundingBox().zmax;
            if (zmax != zmin)
                popZ.addObject(obj);
        }
        return popZ;
    }
   
    
    /**
     * Find cells colocalizing with a nucleus
     */
    public ArrayList<Cell> colocalization(Objects3DIntPopulation cellsPop, Objects3DIntPopulation nucleiPop) {
        ArrayList<Cell> colocPop = new ArrayList<Cell>();
        if (cellsPop.getNbObjects() > 0 && nucleiPop.getNbObjects() > 0) {
            MeasurePopulationColocalisation coloc = new MeasurePopulationColocalisation(nucleiPop, cellsPop);
            for (Object3DInt cell: cellsPop.getObjects3DInt()) {
                for (Object3DInt nucleus: nucleiPop.getObjects3DInt()) {
                    double colocVal = coloc.getValueObjectsPair(nucleus, cell);
                    if (colocVal > 0.5*nucleus.size()) {
                        Object3DComputation objComputation = new Object3DComputation​(cell);
                        Object3DInt cytoplasm = objComputation.getObjectSubtracted(nucleus);
                        cytoplasm.setLabel(nucleus.getLabel());
                        colocPop.add(new Cell(cell, nucleus, cytoplasm));
                        break;
                    }
                }
            }
        }
        return(colocPop);
    }
    
    
    /*
     * Reset labels of cells in population
     */
    public void resetLabels(ArrayList<Cell> cellPop) {
        float label = 1;
        for (Cell cell: cellPop) {
            cell.cell.setLabel(label);
            cell.nucleus.setLabel(label);
            cell.cytoplasm.setLabel(label);
            label++;
        }
    }
    
        
    /**
     * Apply StarDist 2D slice by slice
     * Label detections in 3D
     * @return objects population
     */
    public Objects3DIntPopulation stardistDetection(ImagePlus img, boolean resize, String modelName, double probThresh, double overlayThresh, 
            double minVol, double maxVol) throws IOException{
        
        ImagePlus imgDup = null;
        if (resize && (img.getWidth() > 1024)) {
            float factor = 0.5f;
            imgDup = img.resize((int)(img.getWidth()*factor), (int)(img.getHeight()*factor), 1, "none");
        } else {
            imgDup = new Duplicator().run(img);
            resize = false;
        }
        // StarDist
        File starDistModelFile = new File(modelsPath+File.separator+modelName);
        StarDist2D star = new StarDist2D(syncObject, starDistModelFile);
        star.loadInput(imgDup);
        star.setParams(stardistPercentileBottom, stardistPercentileTop, probThresh, overlayThresh, stardistOutput);
        star.run();
        flush_close(imgDup);
        
        // Label detections in 3D
        ImagePlus imgLabels = (resize) ? star.associateLabels().resize(img.getWidth(), img.getHeight(), 1, "none") : star.associateLabels();
        ImageInt label3D = ImageInt.wrap(imgLabels);
        label3D.setCalibration(cal);
        flush_close(imgLabels);        
        Objects3DIntPopulation pop = new Objects3DIntPopulation(label3D);
        System.out.println(pop.getNbObjects() + " StarDist detections");
        Objects3DIntPopulation popFilter = new Objects3DIntPopulationComputation(pop).getFilterSize(minVol/pixVol, maxVol/pixVol);
        System.out.println(popFilter.getNbObjects() + " detections remaining after size filtering (" + (pop.getNbObjects()-popFilter.getNbObjects()) + " filtered out)");
        popFilter.resetLabels();
        return(popFilter);
    }
    
    
    /**
     * Set the inner/outer ring of nuclei in cell population
     * @param pop
     * @param ringXY
     * @param dil
     */
    public void setNucleiRing(ArrayList<Cell> cellsPop, float dilCoef, boolean dil) {
        dilCoef = (float) (dilCoef / cal.pixelWidth);
        for (Cell cell: cellsPop) {
            if (dil) {
                Object3DInt nucDil = new Object3DComputation(cell.nucleus).getObjectDilated(dilCoef, dilCoef, 0);
                Object3DComputation objComputation = new Object3DComputation​(nucDil);
                Object3DInt donut = objComputation.getObjectSubtracted(cell.nucleus);
                donut.setLabel(cell.nucleus.getLabel());
                cell.setOuterRing(donut);
            } else {
                Object3DInt nucErod = new Object3DComputation(cell.nucleus).getObjectEroded(dilCoef, dilCoef, 0);
                nucErod.setLabel(cell.nucleus.getLabel());
                cell.setInnerNucleus(nucErod);
                
                Object3DComputation objComputation = new Object3DComputation​(cell.nucleus);
                Object3DInt donut = objComputation.getObjectSubtracted(nucErod);
                donut.setLabel(cell.nucleus.getLabel());
                cell.setInnerRing(donut);
            }
        }
    }
       
    
    /**
     * Do Z projection
     * @param img
     * @param projection parameter
     */
    public ImagePlus doZProjection(ImagePlus img, int param) {
        ZProjector zproject = new ZProjector();
        zproject.setMethod(param);
        zproject.setStartSlice(1);
        zproject.setStopSlice(img.getNSlices());
        zproject.setImage(img);
        zproject.doProjection();
       return(zproject.getProjection());
    }
    
    
    /**
     * Find background image intensity:
     * Z projection over min intensity + read mean intensity
     * @param img
     */
    public double findBackground(ImagePlus img) {
      ImagePlus imgProj = doZProjection(img, ZProjector.MIN_METHOD);
      ImageProcessor imp = imgProj.getProcessor();
      double bg = imp.getStatistics().median;
      System.out.println("Background (median intensity of the min projection) = " + bg);
      flush_close(imgProj);
      return(bg);
    }
   
    
    /**
     * Find dots colocalizing with the object
     * @param allDots
     * @param obj
     */
     public Objects3DIntPopulation findDotsPop(Objects3DIntPopulation allDots, Object3DInt obj) {
        Objects3DIntPopulation objPop = new Objects3DIntPopulation();
        objPop.addObject(obj);
        Objects3DIntPopulation dotsIn = new Objects3DIntPopulation();
        MeasurePopulationColocalisation coloc = new MeasurePopulationColocalisation(allDots, objPop);
        for (Object3DInt dot: allDots.getObjects3DInt()) {
            if (coloc.getValueObjectsPair(dot, objPop.getFirstObject()) > 0.5*dot.size()) {
                dotsIn.addObject(dot);
            }
        }
        return(dotsIn);
    }
    
     
    /**
     * Find objects sum of volume
     */
    private double findObjectsVolume(Objects3DIntPopulation pop) {
        if (pop.getNbObjects() > 0) {
            double objsVol = 0;
            for (Object3DInt obj: pop.getObjects3DInt()) {
                objsVol += new MeasureVolume(obj).getVolumeUnit();
            }
            return(objsVol);
        }
        return 0;
    }
    
    
    /**
     * Find objects sum of intensity
     */
    private double findObjectsIntensity(Objects3DIntPopulation pop, ImageHandler imh) {
        if (pop.getNbObjects() > 0) {
            double objsInt = 0;
            for (Object3DInt obj : pop.getObjects3DInt()) {
                objsInt += new MeasureIntensity(obj, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
            }
            return(objsInt);
        }
        return 0;
    }
    
    
    /**
     * Tag cell with parameters....
     * @param img
     * @param cellsPop
     * @param dotsPop
     */
    public void tagCells(ImagePlus img, ArrayList<Cell> cellsPop, Objects3DIntPopulation dotsPop) {
        ImageHandler imh = ImageHandler.wrap(img);
        for (Cell cell: cellsPop) {
            // Get nucleus parameters
            double nucVol = new MeasureVolume(cell.nucleus).getVolumeUnit();
            double nucComp = new MeasureCompactness(cell.nucleus).getValueMeasurement(MeasureCompactness.COMP_CORRECTED);
            double nucSph = new MeasureCompactness(cell.nucleus).getValueMeasurement(MeasureCompactness.SPHER_CORRECTED);
            double nucElongation = new MeasureEllipsoid(cell.nucleus).getValueMeasurement(MeasureEllipsoid.ELL_ELONGATION);
            double nucFlatness = new MeasureEllipsoid(cell.nucleus).getValueMeasurement(MeasureEllipsoid.ELL_FLATNESS);
            double nucInt = new MeasureIntensity(cell.nucleus, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
            // Find dots in nucleus
            Objects3DIntPopulation nucDotsPop = findDotsPop(dotsPop, cell.nucleus);
            double nucDotsVol = findObjectsVolume(nucDotsPop);
            double nucDotsInt = findObjectsIntensity(nucDotsPop, imh);
            
            // Get inner nucleus parameters
            double innerNucVol = new MeasureVolume(cell.innerNucleus).getVolumeUnit();
            double innerNucInt = new MeasureIntensity(cell.innerNucleus, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
            // Find dots in inner nucleus
            Objects3DIntPopulation innerNucDotsPop = findDotsPop(dotsPop, cell.innerNucleus);
            double innerNucDotsVol = findObjectsVolume(innerNucDotsPop);
            double innerNucDotsInt = findObjectsIntensity(innerNucDotsPop, imh);
            
            // Get inner ring parameters
            double innerRingNucVol = new MeasureVolume(cell.innerRing).getVolumeUnit();
            double innerRingNucInt = new MeasureIntensity(cell.innerRing, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
            // Find dots in inner ring
            Objects3DIntPopulation innerRingNucDotsPop = findDotsPop(dotsPop, cell.innerRing);
            double innerRingNucDotsVol = findObjectsVolume(innerRingNucDotsPop);
            double innerRingNucDotsInt = findObjectsIntensity(innerRingNucDotsPop, imh);
            
            // Get outer ring parameters
            double outerRingVol = new MeasureVolume(cell.outerRing).getVolumeUnit();
            double outerRingInt = new MeasureIntensity(cell.outerRing, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
            // Find dots in outer ring
            Objects3DIntPopulation outerRingDotsPop = findDotsPop(dotsPop, cell.outerRing);
            double outerRingDotsVol = findObjectsVolume(outerRingDotsPop);
            double outerRingDotsInt = findObjectsIntensity(outerRingDotsPop, imh);
            
            // Get cell cytoplasm parameters
            double cytoVol = new MeasureVolume(cell.cytoplasm).getVolumeUnit();
            double cytoInt = new MeasureIntensity(cell.cytoplasm, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
            // Find dots in cytoplasm
            Objects3DIntPopulation cytoDotsPop = findDotsPop(dotsPop, cell.cytoplasm);
            double cytoDotsVol = findObjectsVolume(cytoDotsPop);
            double cytoDotsInt = findObjectsIntensity(cytoDotsPop, imh);
            
            // Add all parameters to cell
            cell.setParams(cell.nucleus.getLabel(), nucVol, nucComp, nucSph, nucElongation, nucFlatness, nucInt, 
                    nucDotsPop.getNbObjects(), nucDotsVol, nucDotsInt, innerNucVol, innerNucInt, innerNucDotsPop.getNbObjects(),
                    innerNucDotsVol, innerNucDotsInt, innerRingNucVol, innerRingNucInt, innerRingNucDotsPop.getNbObjects(), 
                    innerRingNucDotsVol, innerRingNucDotsInt, outerRingVol, outerRingInt, outerRingDotsPop.getNbObjects(), 
                    outerRingDotsVol, outerRingDotsInt, cytoVol, cytoInt, cytoDotsPop.getNbObjects(), cytoDotsVol, cytoDotsInt);
        }
        imh.closeImagePlus();
    }
    
    
    public void drawPop(ImagePlus imgIn, Objects3DIntPopulation pop, String outDir, String imgName) {
        ImagePlus img = new Duplicator().run(imgIn);
        ImageHandler imgObj1 = ImageHandler.wrap(img).createSameDimensions();
        pop.drawInImage(imgObj1);
        ImagePlus[] imgColors1 = {null, null, imgObj1.getImagePlus(), img};
        ImagePlus imgObjects1 = new RGBStackMerge().mergeHyperstacks(imgColors1, false);
        imgObjects1.setCalibration(img.getCalibration());
        FileSaver ImgObjectsFile1 = new FileSaver(imgObjects1);
        ImgObjectsFile1.saveAsTiff(outDir + imgName + ".tif");
    }
    
    
    /**
     * Label object
     * @param popObj
     * @param img 
     * @param fontSize 
     */
    public void labelObject(Object3DInt obj, ImagePlus img, int fontSize) {
        if (IJ.isMacOSX())
            fontSize *= 3;
        
        BoundingBox bb = obj.getBoundingBox();
        int z = (int)(bb.zmin + 0.5*(bb.zmax - bb.zmin)); //bb.zmin;
        int x = (int)(bb.xmin + 0.5*(bb.xmax - bb.xmin)); //bb.xmin - 1;
        int y = (int)(bb.ymin + 0.5*(bb.ymax - bb.ymin)); //bb.ymin - 1;
        img.setSlice(z); //z+1
        ImageProcessor ip = img.getProcessor();
        ip.setFont(new Font("SansSerif", Font.PLAIN, fontSize));
        ip.setColor(255);
        ip.drawString(Integer.toString((int)obj.getLabel()), x, y);
        img.updateAndDraw();
    }
    
    
    /**
     * Save image objects
     * @param cellsPop
     * @param dotsPop
     * @param imgDAPI
     * @param imgORF1P
     * @param imgName 
     * @param outDir
     * @param fontSize 
     */
    public void drawResults(ArrayList<Cell> cellsPop, Objects3DIntPopulation dotsPop, ImagePlus imgDAPI, ImagePlus imgORF1P, String imgName, String outDir, int fontSize) {
        ImageHandler imgObj1 = ImageHandler.wrap(new Duplicator().run(imgDAPI)).createSameDimensions();
        ImageHandler imgObj2 = imgObj1.createSameDimensions();
        ImageHandler imgObj3 = imgObj1.createSameDimensions();
        ImageHandler imgObj4 = imgObj1.createSameDimensions();
        ImageHandler imgObj5 = imgObj1.createSameDimensions();
        ImageHandler imgObj6 = imgObj1.createSameDimensions();
        ImageHandler imgObj7 = imgObj1.createSameDimensions();
        ImageHandler imgObj8 = imgObj1.createSameDimensions();
        if (cellsPop.size() > 0) {
            for (Cell cell: cellsPop) {
                cell.nucleus.drawObject(imgObj1, 255);
                cell.cytoplasm.drawObject(imgObj2, 255);
                
                cell.nucleus.drawObject(imgObj3);
                cell.cytoplasm.drawObject(imgObj4);
                labelObject(cell.nucleus, imgObj8.getImagePlus(), fontSize);
                
                cell.innerRing.drawObject(imgObj5, 255);
                cell.outerRing.drawObject(imgObj6, 255);
            }
        }
        if (dotsPop.getNbObjects() > 0) {
            for (Object3DInt dot: dotsPop.getObjects3DInt()) {
                dot.drawObject(imgObj7, 255);
            }
        }
       
        ImagePlus[] imgColors1 = {imgObj7.getImagePlus(), imgObj2.getImagePlus(), imgObj1.getImagePlus(), imgORF1P};
        ImagePlus imgObjects1 = new RGBStackMerge().mergeHyperstacks(imgColors1, true);
        imgObjects1.setCalibration(imgORF1P.getCalibration());
        FileSaver ImgObjectsFile1 = new FileSaver(imgObjects1);
        ImgObjectsFile1.saveAsTiff(outDir + imgName + "_cells.tif"); 
        flush_close(imgObjects1);
        
        ImagePlus[] imgColors2 = {imgObj7.getImagePlus(), null, null, imgDAPI, imgObj5.getImagePlus(), null, imgObj6.getImagePlus()};
        ImagePlus imgObjects2 = new RGBStackMerge().mergeHyperstacks(imgColors2, true);
        imgObjects2.setCalibration(imgDAPI.getCalibration());
        FileSaver ImgObjectsFile2 = new FileSaver(imgObjects2);
        ImgObjectsFile2.saveAsTiff(outDir + imgName + "_rings.tif");
        flush_close(imgObjects2);
        
        ImagePlus[] imgColors3 = {null, imgObj4.getImagePlus(), imgObj3.getImagePlus(), imgObj8.getImagePlus()};
        ImagePlus imgObjects3 = new RGBStackMerge().mergeHyperstacks(imgColors3, true);
        imgObjects3.setCalibration(imgORF1P.getCalibration());
        FileSaver ImgObjectsFile3 = new FileSaver(imgObjects3);
        ImgObjectsFile3.saveAsTiff(outDir + imgName + "_labels.tif"); 
        flush_close(imgObjects3);
        
        imgObj1.closeImagePlus();
        imgObj2.closeImagePlus();
        imgObj3.closeImagePlus();
        imgObj4.closeImagePlus();
        imgObj5.closeImagePlus();
        imgObj6.closeImagePlus();
        imgObj7.closeImagePlus();
    }
}
