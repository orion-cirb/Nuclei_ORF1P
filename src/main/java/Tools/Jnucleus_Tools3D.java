package Tools;

import Cellpose.CellposeSegmentImgPlusAdvanced;
import Cellpose.CellposeTaskSettings;
import StardistOrion.StarDist2D;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.PointRoi;
import ij.gui.WaitForUserDialog;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.plugin.RGBStackMerge;
import ij.plugin.ZProjector;
import ij.plugin.filter.RankFilters;
import ij.process.ImageProcessor;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Point3D;
import mcib3d.geom2.BoundingBox;
import mcib3d.geom2.Object3DComputation;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.Objects3DIntPopulationComputation;
import mcib3d.geom2.measurements.MeasureCentroid;
import mcib3d.geom2.measurements.MeasureCompactness;
import mcib3d.geom2.measurements.MeasureEllipsoid;
import mcib3d.geom2.measurements.MeasureIntensity;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.geom2.measurementsPopulation.MeasurePopulationColocalisation;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import org.apache.commons.io.FilenameUtils;
import org.scijava.util.ArrayUtils;



/**
 * @author phm
 */
public class Jnucleus_Tools3D {
    public final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
    
    public Calibration cal = new Calibration();
    public float pixVol = 0;
    public double minNucVol= 50;
    public double maxNucVol = 1000;
    public float innerNucDil = 1;
    public float outerNucDil = 2;
    public double minCellVol= 100;
    public double maxCellVol = 2000;
    public boolean dotsDetect = true;
    public double minDotVol= 0.5;
    public double maxDotVol = 800;
    
    // StarDist
    public File modelsPath = new File(IJ.getDirectory("imagej")+File.separator+"models");
    public String stardistModelNuc = "dsb2018_heavy_augment.zip";
    public String stardistModelDots = "pmls2.zip";
    public Object syncObject = new Object();
    public final double stardistPercentileBottom = 0.2;
    public final double stardistPercentileTop = 99.8;
    public final double stardistProbThreshNuc = 0.5;
    public final double stardistOverlayThreshNuc = 0.0;
    public final double stardistProbThreshDots = 0.02;
    public final double stardistOverlayThreshDots = 0.25;
    public String stardistOutput = "Label Image"; 
    
    // Cellpose
    public int cellPoseDiameter = 100;
    public String cellPoseModel = "cyto2";
    private String cellPoseEnvDirPath = "/opt/miniconda3/envs/cellpose";
    private String[] channelsName = {"DAPI", "ORF1P"};    
    
 
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
        int index = ArrayUtils.indexOf(modelList, new File(modelsPath+File.separator+stardistModelNuc));
        if (index == -1) {
            IJ.showMessage("Error", stardistModelNuc + " StarDist model not found, please add it in Fiji models folder");
            return false;
        }
        index = ArrayUtils.indexOf(modelList, new File(modelsPath+File.separator+stardistModelDots));
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
        System.out.println("XY calibration = " + cal.pixelWidth + ", Z calibration =" + cal.pixelDepth);
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
        if (IJ.isWindows())
            cellPoseEnvDirPath = System.getProperty("user.home")+"\\miniconda3\\envs\\CellPose";
        gd.addDirectoryField("Cellpose environment directory: ", cellPoseEnvDirPath);
        gd.addNumericField("Min cell volume (µm3): ", minCellVol);
        gd.addNumericField("Max cell volume (µm3): ", maxCellVol);
        
        gd.addMessage("Dots detection", Font.getFont("Monospace"), Color.blue);
        gd.addCheckbox("Detect dots", dotsDetect);
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
        
        imgOut.show();
        new WaitForUserDialog("test").show();
       
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
     * @return cells cytoplasm
     */
    public Objects3DIntPopulation colocalization(Objects3DIntPopulation cellsPop, Objects3DIntPopulation nucleiPop) {
        Objects3DIntPopulation colocPop = new Objects3DIntPopulation();
        if (cellsPop.getNbObjects() > 0 && nucleiPop.getNbObjects() > 0) {
            MeasurePopulationColocalisation coloc = new MeasurePopulationColocalisation(nucleiPop, cellsPop);
            for (Object3DInt cell: cellsPop.getObjects3DInt()) {
                for (Object3DInt nucleus: nucleiPop.getObjects3DInt()) {
                    double colocVal = coloc.getValueObjectsPair(nucleus, cell);
                    if (colocVal > 0.5*nucleus.size()) {
                        Object3DComputation objComputation = new Object3DComputation​(cell);
                        Object3DInt cytoplasm = objComputation.getObjectSubtracted(nucleus);
                        cytoplasm.setLabel(nucleus.getLabel());
                        colocPop.addObject(cytoplasm);
                        break;
                    }
                }
            }
        }
        colocPop.setVoxelSizeXY(cal.pixelWidth);
        colocPop.setVoxelSizeZ(cal.pixelDepth);
        return(colocPop);
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
        
        label3D.show();
        new WaitForUserDialog("test").show();
        
        Objects3DIntPopulation objPop = new Objects3DIntPopulationComputation(new Objects3DIntPopulation(label3D)).getFilterSize(minVol/pixVol, maxVol/pixVol);
        objPop.resetLabels();
        flush_close(imgLabels);
        return(objPop);
    }
    
    
    /**
     * Return objects population in a binary image
     */
    public Objects3DIntPopulation getPopFromImage(ImagePlus img) {
        ImageLabeller labeller = new ImageLabeller();
        ImageInt labels = labeller.getLabels(ImageHandler.wrap(img));
        Objects3DIntPopulation pop = new Objects3DIntPopulation(labels);
        return pop;
    }
    
    
    /**
     * Create donut object population
     * @param pop
     * @param img
     * @param ringXY
     * @param dil
     */
    public Objects3DIntPopulation createDonutPop(Objects3DIntPopulation pop, ImagePlus img, float dilCoef, boolean dil) {
        dilCoef = (float) (dilCoef / cal.pixelWidth);
        ImagePlus imgCopy = new Duplicator().run(img);
        ImageInt imgBin = ImageInt.wrap(imgCopy);
        Objects3DIntPopulation donutPop = new Objects3DIntPopulation();
        for (Object3DInt obj: pop.getObjects3DInt()) {
            imgBin.fill(0);
            if (dil) {
                Object3DInt objDil = new Object3DComputation(obj).getObjectDilated(dilCoef, dilCoef, 0);
                objDil.drawObject(imgBin, 255);
                obj.drawObject(imgBin, 0);
            } else {
                Object3DInt objErod = new Object3DComputation(obj).getObjectEroded(dilCoef, dilCoef, 0);
                obj.drawObject(imgBin, 255);
                objErod.drawObject(imgBin, 0);
            }
            Objects3DIntPopulation tmpPop = getPopFromImage(imgBin.getImagePlus());
            Object3DInt objD = tmpPop.getFirstObject();
            objD.setLabel(obj.getLabel());
            donutPop.addObject(objD);
        }
        flush_close(imgCopy);
        imgBin.closeImagePlus();
        return(donutPop);
    }
    
    
    /**
     * Get inner nucleus
     */
    public Objects3DIntPopulation getInnerNucleus(Objects3DIntPopulation pop, ImagePlus img) {
        Objects3DIntPopulation innerNucleusPop = new Objects3DIntPopulation();
        float erode = (float) (innerNucDil/cal.pixelWidth);
        for (Object3DInt obj: pop.getObjects3DInt()) {
            Object3DInt objEr = new Object3DComputation(obj).getObjectEroded(erode, erode, 0);
            objEr.setLabel(obj.getLabel());
            innerNucleusPop.addObject(objEr);
        }
        innerNucleusPop.setVoxelSizeXY(cal.pixelWidth);
        innerNucleusPop.setVoxelSizeZ(cal.pixelDepth);
        return(innerNucleusPop);
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
      double[] bg = new double[2];
      ImagePlus imgProj = doZProjection(img, ZProjector.MIN_METHOD);
      ImageProcessor imp = imgProj.getProcessor();
      bg[0] = imp.getStatistics().mean;
      bg[1] = imp.getStatistics().stdDev;
      System.out.println("Background = " + bg[0] + " +- " + bg[1]);
      flush_close(imgProj);
      return(bg[0]+bg[1]);
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
     * Tags cell with  parameters....
     * @param img
     * @param nucPop
     * @param innerNucPop
     * @param innerRingPop
     * @param outerRingPop
     * @param cellsPop
     * @param allDots
     */
    public ArrayList<Nucleus> tagNuclei(ImagePlus img, Objects3DIntPopulation nucPop, Objects3DIntPopulation innerNucPop, Objects3DIntPopulation innerRingPop,
            Objects3DIntPopulation outerRingPop, Objects3DIntPopulation cellsPop, Objects3DIntPopulation allDots) {
        
        if (allDots.getNbObjects() > 0)
            allDots.drawInImage(ImageHandler.wrap(img).createSameDimensions()); //???
        
        ArrayList<Nucleus> nuclei = new ArrayList<>();
        ImageHandler imh = ImageHandler.wrap(img);
        
        for (Object3DInt nucObj: nucPop.getObjects3DInt()) {
            // Get nucleus parameters
            double nucVol = new MeasureVolume(nucObj).getVolumeUnit();
            double nucComp = new MeasureCompactness(nucObj).getValueMeasurement(MeasureCompactness.COMP_CORRECTED);
            double nucSph = new MeasureCompactness(nucObj).getValueMeasurement(MeasureCompactness.SPHER_CORRECTED);
            double nucElongation = new MeasureEllipsoid(nucObj).getValueMeasurement(MeasureEllipsoid.ELL_ELONGATION);
            double nucFlatness = new MeasureEllipsoid(nucObj).getValueMeasurement(MeasureEllipsoid.ELL_FLATNESS);
            double nucInt = new MeasureIntensity(nucObj, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
            // Find dots in nucleus
            Objects3DIntPopulation nucDotsPop = findDotsPop(allDots, nucObj);
            double nucDotsVol = findObjectsVolume(nucDotsPop);
            double nucDotsInt = findObjectsIntensity(nucDotsPop, imh);
            
            // Get inner nucleus parameters
            Object3DInt innerNucObj = innerNucPop.getObjectByLabel(nucObj.getLabel());
            double innerNucVol = new MeasureVolume(innerNucObj).getVolumeUnit();
            double innerNucInt = new MeasureIntensity(innerNucObj, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
            // Find dots in inner nucleus
            Objects3DIntPopulation innerNucDotsPop = findDotsPop(allDots, innerNucObj);
            double innerNucDotsVol = findObjectsVolume(innerNucDotsPop);
            double innerNucDotsInt = findObjectsIntensity(innerNucDotsPop, imh);
            
            // Get inner ring parameters
            Object3DInt innerRingNucObj = innerRingPop.getObjectByLabel(nucObj.getLabel());
            double innerRingNucVol = new MeasureVolume(innerRingNucObj).getVolumeUnit();
            double innerRingNucInt = new MeasureIntensity(innerRingNucObj, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
            // Find dots in inner ring
            Objects3DIntPopulation innerRingNucDotsPop = findDotsPop(allDots, innerRingNucObj);
            double innerRingNucDotsVol = findObjectsVolume(innerRingNucDotsPop);
            double innerRingNucDotsInt = findObjectsIntensity(innerRingNucDotsPop, imh);
            
            // Get outer ring parameters
            Object3DInt outerRingObj = outerRingPop.getObjectByLabel(nucObj.getLabel());
            double outerRingVol = new MeasureVolume(outerRingObj).getVolumeUnit();
            double outerRingInt = new MeasureIntensity(outerRingObj, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
            // Find dots in outer ring
            Objects3DIntPopulation outerRingDotsPop = findDotsPop(allDots, outerRingObj);
            double outerRingDotsVol = findObjectsVolume(outerRingDotsPop);
            double outerRingDotsInt = findObjectsIntensity(outerRingDotsPop, imh);
            
            // Get cell cytoplasm parameters
            Object3DInt cellObj = cellsPop.getObjectByLabel(nucObj.getLabel());
            double cytoVol = 0;
            double cytoInt = 0;
            int cytoDotsNb = 0;
            double cytoDotsVol = 0;
            double cytoDotsInt = 0;
            if (cellObj != null) {
                cytoVol = new MeasureVolume(cellObj).getVolumeUnit();
                cytoInt = new MeasureIntensity(cellObj, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
                Objects3DIntPopulation cytoDotsPop = findDotsPop(allDots, cellObj);
                cytoDotsNb = cytoDotsPop.getNbObjects();
                cytoDotsVol = findObjectsVolume(cytoDotsPop);
                cytoDotsInt = findObjectsIntensity(cytoDotsPop, imh);
            }
            
            // Add all parameters to cell
            Nucleus nucleus = new Nucleus((int) nucObj.getLabel(), nucVol, nucComp, nucSph, nucElongation, nucFlatness, nucInt, 
                    nucDotsPop.getNbObjects(), nucDotsVol, nucDotsInt, innerNucVol, innerNucInt, innerNucDotsPop.getNbObjects(),
                    innerNucDotsVol, innerNucDotsInt, innerRingNucVol, innerRingNucInt, innerRingNucDotsPop.getNbObjects(), 
                    innerRingNucDotsVol, innerRingNucDotsInt, outerRingVol, outerRingInt, outerRingDotsPop.getNbObjects(), 
                    outerRingDotsVol, outerRingDotsInt, cytoVol, cytoInt, cytoDotsNb, cytoDotsVol, cytoDotsInt);            
            nuclei.add(nucleus);
        }
        
        imh.closeImagePlus();
        return(nuclei);
    }
    
    
    /**
     * Label object
     * @param popObj
     * @param img 
     */
    public void labelsObject(Objects3DIntPopulation popObj, ImagePlus img, int fontSize) {
        if (IJ.isMacOSX())
            fontSize *= 3;
        Font tagFont = new Font("SansSerif", Font.PLAIN, fontSize);
        for (Object3DInt obj : popObj.getObjects3DInt()) {
            BoundingBox bb = obj.getBoundingBox();
            int z = bb.zmax - bb.zmin;
            int x = bb.xmin - 2;
            int y = bb.ymin - 2;
            img.setSlice(z+1);
            ImageProcessor ip = img.getProcessor();
            ip.setFont(tagFont);
            ip.setColor(255);
            ip.drawString(Float.toString(obj.getLabel()), x, y);
            img.updateAndDraw();
        }
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
     * Save image objects
     * @param pop1
     * @param pop2
     * @param pop3
     * @param img
     * @param name 
     * @param fontSize 
     */
    public void saveImageObjects(Objects3DIntPopulation pop1, Objects3DIntPopulation pop2, Objects3DIntPopulation pop3, ImagePlus img, String name, int labelsCh, int fontSize) {
        ImagePlus imgObj = new Duplicator().run(img);
        ImageHandler imgObj1 = ImageHandler.wrap(imgObj).createSameDimensions();
        ImageHandler imgObj2 = imgObj1.duplicate();
        ImageHandler imgObj3 = imgObj1.duplicate();
        // draw obj population
        if (pop1 != null) {
            pop1.drawInImage(imgObj1);
            if (labelsCh == 1)
                 labelsObject(pop1, imgObj1.getImagePlus(), fontSize);
        }
        if (pop2 != null) {
            pop2.drawInImage(imgObj2);
            if (labelsCh == 2)
                labelsObject(pop2, imgObj2.getImagePlus(), fontSize);
        }
        if (pop3 != null) {
            pop3.drawInImage(imgObj3);
            if (labelsCh == 3)
                labelsObject(pop3, imgObj3.getImagePlus(), fontSize);
        }
       
        ImagePlus[] imgColors = {imgObj1.getImagePlus(), imgObj2.getImagePlus(), imgObj3.getImagePlus(), imgObj};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, true);
        FileSaver ImgObjectsFile = new FileSaver(imgObjects);
        ImgObjectsFile.saveAsTiff(name); 
        imgObj1.closeImagePlus(); 
        imgObj2.closeImagePlus();
        imgObj3.closeImagePlus();
        flush_close(imgObjects);
    }
}