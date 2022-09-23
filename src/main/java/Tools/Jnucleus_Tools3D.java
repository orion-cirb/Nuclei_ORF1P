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
import mcib3d.geom.Voxel3D;
import mcib3d.geom2.BoundingBox;
import mcib3d.geom2.Object3DComputation;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Object3DPlane;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.Objects3DIntPopulationComputation;
import mcib3d.geom2.VoxelInt;
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



/**
 *
 * @author phm
 */

public class Jnucleus_Tools3D {
    

    public double minNucVol= 25;
    public double maxNucVol = 500;
    public double minDotVol= 0.05;
    public double maxDotVol = 50;
    public float innerNucDil = 1;
    public float outerNucDil = 2;
    public boolean dotsDetect = true;
    public Calibration cal = new Calibration();
    public float pixVol = 0;
    
    // StarDist
    public Object syncObject = new Object();
    public final double stardistPercentileBottom = 0.2;
    public final double stardistPercentileTop = 99.8;
    public final double stardistProbThreshNuc = 0.65;
    public final double stardistOverlayThreshNuc = 0.01;
    public File modelsPath = new File(IJ.getDirectory("imagej")+File.separator+"models");
    public String stardistModel = "";
    public String stardistOutput = "Label Image"; 
    
    // Cellpose
    public int cellPoseDiameter = 80;
    private boolean useGpu = true;
    private String[] cellposeModels = {"cyto","nuclei","tissuenet","livecell", "cyto2", "general","CP", "CPx", "TN1", "TN2", "TN3", "LC1",
        "LC2", "LC3", "LC4"};
    public String cellModel = "";
    private String cellPoseEnvDirPath = "/home/phm/.conda/envs/cellpose";
    public String cellDetector = "";
    
    
    
    public final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
 
    /*
    Find starDist models in Fiji models folder
    */
    public String[] findStardistModels() {
        FilenameFilter filter = (dir, name) -> name.endsWith(".zip");
        File[] modelList = modelsPath.listFiles(filter);
        String[] models = new String[modelList.length];
        for (int i = 0; i < modelList.length; i++) {
            models[i] = modelList[i].getName();
        }
        Arrays.sort(models);
        return(models);
    } 
    
    
    /**
     *
     * @param img
     */
    public void closeImages(ImagePlus img) {
        img.flush();
        img.close();
    }
    
    
  /**
     * return objects population in an binary image
     * @param img
     * @return pop objects population
     */

    public  Objects3DIntPopulation getPopFromImage(ImagePlus img) {
        // label binary images first
        ImageLabeller labeller = new ImageLabeller();
        ImageInt labels = labeller.getLabels(ImageHandler.wrap(img));
        Objects3DIntPopulation pop = new Objects3DIntPopulation(labels);
        return pop;
    }
    
    
    /*Median filter 
     * 
     * @param img
     * @param size
     */ 
    public void median_filter(ImagePlus img, double size) {
        RankFilters median = new RankFilters();
        for (int s = 1; s <= img.getNSlices(); s++) {
            img.setZ(s);
            median.rank(img.getProcessor(), size, RankFilters.MEDIAN);
            img.updateAndDraw();
        }
    }
            
  
    /**
     * Ask for parameters
     * @param channels
     * @param channelsName
     * @return 
     */
    
    public ArrayList dialog(List<String> channels, List<String> channelsName) {
        ArrayList ch = new ArrayList();
        String[] models = findStardistModels();
        String[] cellsDetectors = {"CellPose", "CellOutiner"};
        if (IJ.isWindows())
            cellPoseEnvDirPath = System.getProperty("user.home")+"\\miniconda3\\envs\\CellPose";
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 80, 0);
        gd.addImage(icon);
        gd.addMessage("Choose channels", new Font(Font.MONOSPACED , Font.BOLD, 12), Color.blue);
        int index = 0;
        for (String chName : channelsName) {
            gd.addChoice(chName, channels.toArray(new String[0]), channels.get(index));
            index++;
        }
        gd.addMessage("Nucleus parameters", new Font(Font.MONOSPACED , Font.BOLD, 12), Color.blue);
        gd.addMessage("StarDist model", Font.getFont("Monospace"), Color.blue);
        if (models.length > 0) {
            gd.addChoice("StarDist model :",models, models[0]);
        }
        else {
            gd.addMessage("No StarDist model found in Fiji !!", Font.getFont("Monospace"), Color.red);
            gd.addFileField("StarDist model :", stardistModel);
        }
        gd.addNumericField("Min nucleus vol. :", minNucVol);
        gd.addNumericField("Max nucleus vol. :", maxNucVol);   
        gd.addMessage("Nucleus dounuts", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Nucleus cyto ring (µm) :", outerNucDil);
        gd.addNumericField("Nucleus inner ring (µm) :", innerNucDil);
        gd.addMessage("Cells parameters", new Font(Font.MONOSPACED , Font.BOLD, 12), Color.blue);
        gd.addChoice("Cells segmentation method :",cellsDetectors, cellsDetectors[0]);
        gd.addDirectoryField("Cellpose environment path : ", cellPoseEnvDirPath);
        gd.addChoice("Cellpose model : ", cellposeModels, cellposeModels[4]);
        gd.addNumericField("Cell size (µm3) : ", cellPoseDiameter, 2);
        gd.addCheckbox("  Do dots detection", dotsDetect);
        gd.addMessage("Image calibration", new Font(Font.MONOSPACED , Font.BOLD, 12), Color.blue);
        gd.addNumericField("XY cal. :", cal.pixelWidth);
        gd.addNumericField("Z cal.  :", cal.pixelDepth);
        gd.showDialog();
        for (int i = 0; i < index; i++)
            ch.add(i, gd.getNextChoice());
        if (models.length > 0) {
            stardistModel = modelsPath+File.separator+gd.getNextChoice();
        }
        else {
            stardistModel = gd.getNextString();
        }
        if (stardistModel.isEmpty()) {
            IJ.error("No model specify !!");
            return(null);
        }
        minNucVol = (float)gd.getNextNumber();
        maxNucVol = (float)gd.getNextNumber();
        outerNucDil = (float)gd.getNextNumber();
        innerNucDil = (float)gd.getNextNumber();
        cellDetector = gd.getNextChoice();
        cellPoseEnvDirPath = gd.getNextString();
        cellModel = gd.getNextChoice();
        cellPoseDiameter = (int)gd.getNextNumber();
        dotsDetect = gd.getNextBoolean();
        cal.pixelWidth = gd.getNextNumber();
        cal.pixelDepth = gd.getNextNumber();
        cal.pixelHeight =  cal.pixelWidth;
        pixVol = (float)(cal.pixelWidth*cal.pixelHeight*cal.pixelDepth);
        if(gd.wasCanceled())
            ch = null;
        return(ch);
    }
    
   
    /**
     * Find images in folder
     */
    public ArrayList findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No Image found in "+imagesFolder);
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
     * Find channels name
     * @param imageName
     * @param imageExt
     */
    public List<String> findChannels (String imageName) throws DependencyException, ServiceException, FormatException, IOException {
        List<String> channels = new ArrayList<>();
        // create OME-XML metadata store of the latest schema version
        ServiceFactory factory;
        factory = new ServiceFactory();
        OMEXMLService service = factory.getInstance(OMEXMLService.class);
        IMetadata meta = service.createOMEXMLMetadata();
        ImageProcessorReader reader = new ImageProcessorReader();
        reader.setMetadataStore(meta);
        reader.setId(imageName);
        int chs = reader.getSizeC();
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                String channelsID = meta.getImageName(0);
                channels = Arrays.asList(channelsID.replace("_", "-").split("/"));
                break;
            case "lif" :
                String[] ch = new String[chs];
                if (chs > 1) {
                    for (int n = 0; n < chs; n++) 
                        if (meta.getChannelExcitationWavelength(0, n) == null)
                            channels.add(Integer.toString(n));
                        else 
                            channels.add(meta.getChannelExcitationWavelength(0, n).value().toString());
                }
                break;
            default :
                chs = reader.getSizeC();
                for (int n = 0; n < chs; n++)
                    channels.add(Integer.toString(n));
        }
        return(channels);         
    }
    
    /**
     * Find image calibration
     */
    public Calibration findImageCalib(IMetadata meta) {
        // read image calibration
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        System.out.println("x cal = " +cal.pixelWidth+", z cal=" + cal.pixelDepth);
        return(cal);
    }
        
    /**
 * Find cells with cellpose
 * return cell cytoplasm
 * @param img
 * @param type
 * @return 
 */
    public Objects3DIntPopulation cellPoseCellsPop(ImagePlus img, Objects3DIntPopulation nucPop){
        CellposeTaskSettings settings = new CellposeTaskSettings(cellModel, 1, cellPoseDiameter, cellPoseEnvDirPath);
        settings.setCellProbTh(0);
        settings.setStitchThreshold(0.25); 
        settings.setFlowTh(0.6);
        settings.useGpu(true);
        ImagePlus imgIn = new Duplicator().run(img);
        CellposeSegmentImgPlusAdvanced cellpose = new CellposeSegmentImgPlusAdvanced(settings, imgIn);
        ImagePlus cellpose_img = cellpose.run(); 
        closeImages(imgIn);
        cellpose_img.setCalibration(cal);
        ImageHandler imh = ImageHandler.wrap(cellpose_img);
        Objects3DIntPopulation pop = new Objects3DIntPopulation(imh);
        imh.closeImagePlus();
        closeImages(cellpose_img);
        // take cell with nucleus
        Objects3DIntPopulation cellPop = new Objects3DIntPopulation();
        MeasurePopulationColocalisation coloc = new MeasurePopulationColocalisation(nucPop, pop);
        for (Object3DInt cell : pop.getObjects3DInt()) {
            for (Object3DInt nuc : nucPop.getObjects3DInt()) {
                if (coloc.getValueObjectsPair(nuc.getLabel(), cell.getLabel()) > 0.5*nuc.size()) {
                    Object3DComputation objComp = new Object3DComputation(cell);
                    Object3DInt cytoObj = objComp.getObjectSubtracted(nuc);
                    cytoObj.setLabel(nuc.getLabel());
                    cellPop.addObject(cytoObj);
                    break;
                }
            }
        }
        cellPop.setVoxelSizeXY(cal.pixelWidth);
        cellPop.setVoxelSizeZ(cal.pixelDepth);
        return(cellPop);
    }
    
    
    
    /*
    Find cell cytoplasm
    */
    public Objects3DIntPopulation findCells (ImagePlus img, Objects3DIntPopulation nucPop) {
        Objects3DIntPopulation cellPop = new Objects3DIntPopulation();
        ImagePlus imgCell = new Duplicator().run(img);
        CellOutliner cell = new CellOutliner();
        cell.cellRadius = 80;
        cell.darkEdge = false;
        cell.dilate = 8;
        cell.ellipticalFit = false;
        cell.iterations = 4;
        cell.kernelSmoothing = 1.2;
        cell.kernelWidth = 7;
        cell.polygonSmoothing = 1;
        cell.tolerance = 5;
        cell.weightingGamma = 2.5;
        cell.processAllSlices = true;
        cell.buildMaskOutput = true;
        for (Object3DInt nucObj : nucPop.getObjects3DInt()) {
            Point3D pt = new MeasureCentroid(nucObj).getCentroidAsPoint();
            imgCell.setSlice(pt.getRoundZ());
            PointRoi roi = new PointRoi(pt.getRoundX(), pt.getRoundY());
            imgCell.setRoi(roi, true);
            cell.setup("", imgCell);
            cell.run(imgCell.getProcessor());
            ImagePlus imgMask = cell.maskImp;
            imgMask.setCalibration(cal);
            nucObj.drawObject(ImageHandler.wrap(imgMask), 0);
            // get the cell cytoplasm
            Object3DInt cellCyto = getPopFromImage(imgMask).getFirstObject();
            cellCyto.setLabel(nucObj.getLabel());
            closeImages(imgMask);
            cellPop.addObject(cellCyto);
        }
        closeImages(imgCell);
        System.out.println(cellPop.getNbObjects()+" cells found");
        return(cellPop);
    }
    
    
    public Objects3DIntPopulation findNucleus(ImagePlus imgNuc, int blur1, int blur2, int radOut, String thMethod) {
        Objects3DIntPopulation nucPopOrg = find_nucleus2(imgNuc, blur1, blur2, radOut, thMethod);
        System.out.println("-- Total nucleus Population :"+nucPopOrg.getNbObjects());
        // size filter
        Objects3DIntPopulation nucPop = new Objects3DIntPopulationComputation(nucPopOrg).getFilterSize(minNucVol/pixVol, maxNucVol/pixVol);
        nucPop.resetLabels();
        int nbNucPop = nucPop.getNbObjects();
        System.out.println("-- Total nucleus Population after size filter: "+ nbNucPop);
        return(nucPop);
    }
        
    /**
     * Nucleus segmentation 2
     * @param imgNuc
     * @return cellPop
     */
    public Objects3DIntPopulation find_nucleus2(ImagePlus imgNuc, int blur1, int blur2, int radOut, String thMethod) {
        ImagePlus img = new Duplicator().run(imgNuc);
        ImageStack stack = new ImageStack(img.getWidth(), imgNuc.getHeight());
        for (int i = 1; i <= img.getStackSize(); i++) {
            IJ.showStatus("Finding nucleus section "+i+" / "+img.getStackSize());
            img.setZ(i);
            img.updateAndDraw();
            IJ.run(img, "Nuclei Outline", "blur="+blur1+" blur2="+blur2+" threshold_method="+thMethod+" outlier_radius="+radOut+" outlier_threshold=1 max_nucleus_size=500 "
                    + "min_nucleus_size=30 erosion=5 expansion_inner=5 expansion=5 results_overlay");
            img.setZ(1);
            img.updateAndDraw();
            ImagePlus mask = new ImagePlus("mask", img.createRoiMask().getBufferedImage());
            ImageProcessor ip =  mask.getProcessor();
            ip.invertLut();
            for (int n = 0; n < 3; n++) 
                ip.erode();
            stack.addSlice(ip);
        }
        ImagePlus imgStack = new ImagePlus("Nucleus", stack); 
        imgStack.setCalibration(cal);
        IJ.run(imgStack, "Watershed", "stack");
        Objects3DIntPopulation nucPop = getPopFromImage(imgStack);
        closeImages(img);
        closeImages(imgStack);
        return(nucPop);
    }
    
    
    /** Look for all nuclei
         Do z slice by slice stardist 
         * return nuclei population
         */
        public Objects3DIntPopulation stardistNucleiPop(ImagePlus imgNuc) throws IOException{
            // resize to be in a stardist-friendly scale
            ImagePlus img = null;
            int width = imgNuc.getWidth();
            int height = imgNuc.getHeight();
            float factor = 0.25f;
            boolean resized = false;
            if (imgNuc.getWidth() > 1024) {
                img = imgNuc.resize((int)(width*factor), (int)(height*factor), 1, "none");
                resized = true;
            }
            else
                img = new Duplicator().run(imgNuc);
            
            IJ.run(img, "Remove Outliers", "block_radius_x=20 block_radius_y=20 standard_deviations=1 which=Dark stack");
            // Clear unfocus Z plan
            Find_focused_slices focus = new Find_focused_slices();
            focus.run(img);
            // Go StarDist
            File starDistModelFile = new File(stardistModel);
            StarDist2D star = new StarDist2D(syncObject, starDistModelFile);
            star.loadInput(img);
            star.setParams(stardistPercentileBottom, stardistPercentileTop, stardistProbThreshNuc, stardistOverlayThreshNuc, stardistOutput);
            star.run();
            closeImages(img);
            // label in 3D
            ImagePlus nuclei = (resized) ? star.associateLabels().resize(width, height, 1, "none") : star.associateLabels();
            ImageInt label3D = ImageInt.wrap(nuclei);
            label3D.setCalibration(cal);
            Objects3DIntPopulation nucPop = new Objects3DIntPopulationComputation(new Objects3DIntPopulation(label3D)).
                    getFilterSize(minNucVol/pixVol, maxNucVol/pixVol);
            nucPop.resetLabels();
            closeImages(nuclei);
            return(nucPop);
        }
    
   

    public ArrayList<Double> readIntensity(ImagePlus img, Objects3DIntPopulation objPop) {
        ArrayList<Double> intensity = new ArrayList();
        ImageHandler imh = ImageHandler.wrap(img);
        for (Object3DInt obj : objPop.getObjects3DInt()) {
            intensity.add(new MeasureIntensity(obj, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM));
        }
        imh.closeImagePlus();
        return(intensity);
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
        closeImages(imgObjects);
    }

    
    /**
     * Do Z projection
     * @param img
     * @param projection parameter
     * @return 
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
    * Find background image intensity
    * Z project min intensity
    * read mean intensity
    * @param img 
    */
    public double find_background(ImagePlus img) {
      double[] bg = new double[2];
      ImagePlus imgProj = doZProjection(img, ZProjector.MIN_METHOD);
      ImageProcessor imp = imgProj.getProcessor();
      bg[0] = imp.getStatistics().mean;
      bg[1] = imp.getStatistics().stdDev;
      System.out.println("Background =  " + bg[0] + "+-" + bg[1]);
      closeImages(imgProj);
      return(bg[0]+bg[1]);
    }
    
    /**
     * Get inner nucleus
     * 
     */
    public Objects3DIntPopulation getInnerNucleus(Objects3DIntPopulation pop, ImagePlus img) {
        Objects3DIntPopulation innerNucleusPop = new Objects3DIntPopulation();
        float erode = (float)(innerNucDil/cal.pixelWidth);
        for (Object3DInt obj : pop.getObjects3DInt()) {
            Object3DInt objEr = new Object3DComputation(obj).getObjectEroded(erode, erode, 0);
            objEr.setLabel(obj.getLabel());
            innerNucleusPop.addObject(objEr);
        }
        innerNucleusPop.setVoxelSizeXY(cal.pixelWidth);
        innerNucleusPop.setVoxelSizeZ(cal.pixelDepth);
        return(innerNucleusPop);
    }
    
    
    /**
     * Create donut object population
     * 
     * @param pop
     * @param img
     * @param ringXY
     * @param dil
     * @return 
     */
    public Objects3DIntPopulation createDonutPop(Objects3DIntPopulation pop, ImagePlus img, float dilCoef, boolean dil) {
        dilCoef = (float)(dilCoef / cal.pixelWidth);
        ImagePlus imgCopy = new Duplicator().run(img);
        ImageInt imgBin = ImageInt.wrap(imgCopy);
        Objects3DIntPopulation donutPop = new Objects3DIntPopulation();
        for (Object3DInt obj : pop.getObjects3DInt()) {
            imgBin.fill(0);
            if (dil) {
                Object3DInt objDil = new Object3DComputation(obj).getObjectDilated(dilCoef, dilCoef, 0);
                objDil.drawObject(imgBin, 255);
                obj.drawObject(imgBin, 0);
            }
            else {
                Object3DInt objErod = new Object3DComputation(obj).getObjectEroded(dilCoef, dilCoef, 0);
                obj.drawObject(imgBin, 255);
                objErod.drawObject(imgBin, 0);
            }
            Objects3DIntPopulation tmpPop = getPopFromImage(imgBin.getImagePlus());
            Object3DInt objD = tmpPop.getFirstObject();
            objD.setLabel(obj.getLabel());
            donutPop.addObject(objD);
        }
        closeImages(imgCopy);
        imgBin.closeImagePlus();
        return(donutPop);
    }
    
   
    
    /**
     * Find dots
     */
    public Objects3DIntPopulation find_dots(ImagePlus img, int sig1, int sig2, String th) {
        ImagePlus imgDup = new Duplicator().run(img);
        //median_filter(imgDup, 1.5);
        IJ.run(imgDup, "Difference of Gaussians", "  sigma1="+sig1+" sigma2=" +sig2+" enhance stack");
        imgDup.setSlice(imgDup.getNSlices()/2);
        IJ.setAutoThreshold(imgDup, th+" dark");
        Prefs.blackBackground = false;
        IJ.run(imgDup, "Convert to Mask", "method="+th+" background=Dark");
        Objects3DIntPopulation pop = new Objects3DIntPopulationComputation(getPopFromImage(imgDup)).getFilterSize(minDotVol, maxDotVol);
        pop.resetLabels();
        pop.setVoxelSizeXY(cal.pixelWidth);
        pop.setVoxelSizeZ(cal.pixelDepth);
        closeImages(imgDup);
        return(pop);
    }
    
    /**
     * Find dots in cells 
     * check if dot centroid inside nucleus
     * @param allDots
     * @param nucObj
     * @return 
     */
     public Objects3DIntPopulation findDotsPop(Objects3DIntPopulation allDots, Object3DInt nucObj) {
        Objects3DIntPopulation nucPop = new Objects3DIntPopulation();
        nucPop.addObject(nucObj);
        Objects3DIntPopulation popDotsIn = new Objects3DIntPopulation();
        MeasurePopulationColocalisation coloc = new MeasurePopulationColocalisation(allDots, nucPop);
        for (Object3DInt dot : allDots.getObjects3DInt()) {
            if (coloc.getValueObjectsPair(dot, nucPop.getFirstObject()) > 0.5*dot.size()) {
                popDotsIn.addObject(dot);
            }
        }
        return(popDotsIn);
    }
    
    /**
     * Find objects sum of volume
     */
    private double findObjectsVolume(Objects3DIntPopulation pop) {
        double objsVol = 0;
        for (Object3DInt obj : pop.getObjects3DInt()) {
            objsVol += new MeasureVolume(obj).getVolumeUnit();
        }
        return(objsVol);
    }
    
    /**
     * Find objects sum of Intensity
     */
    private double findObjectsIntensity(Objects3DIntPopulation pop, ImageHandler imh) {
        double objsInt = 0;
        for (Object3DInt obj : pop.getObjects3DInt()) {
            objsInt += new MeasureIntensity(obj, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
        }
        return(objsInt);
    }
    
    /**
     * Label object
     * @param popObj
     * @param img 
     */
    public void labelsObject (Objects3DIntPopulation popObj, ImagePlus img, int fontSize) {
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
    
   /**
     * Tags cell with  parameters....
     * @param img
     * @param nucPop
     * @param innerNucPop
     * @param innerRingPop
     * @param cytoPop
     * @param allDots
     * @return 
     */
    
    public ArrayList<Nucleus> tagsNuclei(ImagePlus img, Objects3DIntPopulation nucPop, Objects3DIntPopulation innerNucPop, Objects3DIntPopulation innerRingPop,
            Objects3DIntPopulation outerNucPop, Objects3DIntPopulation cellsPop, Objects3DIntPopulation allDots) {
        if (allDots.getNbObjects() > 0)
            allDots.drawInImage(ImageHandler.wrap(img).createSameDimensions());
        ArrayList<Nucleus> nuclei = new ArrayList<>();
        ImageHandler imh = ImageHandler.wrap(img);
        for (Object3DInt nucObj : nucPop.getObjects3DInt()) {
            // nucleus
            double nucVol = new MeasureVolume(nucObj).getVolumeUnit();
            double nucComp = new MeasureCompactness(nucObj).getValueMeasurement(MeasureCompactness.COMP_CORRECTED);
            double nucSph = new MeasureCompactness(nucObj).getValueMeasurement(MeasureCompactness.SPHER_CORRECTED);
            double nucElongation = new MeasureEllipsoid(nucObj).getValueMeasurement(MeasureEllipsoid.ELL_ELONGATION);
            double nucFlatness = new MeasureEllipsoid(nucObj).getValueMeasurement(MeasureEllipsoid.ELL_FLATNESS);
            double nucInt = new MeasureIntensity(nucObj, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
            //Finding dots in nucleus
            Objects3DIntPopulation nucDotsPop = findDotsPop(allDots, nucObj);
            int nucDots = nucDotsPop.getNbObjects();
            double nucDotsVol = 0;
            double nucDotsInt = 0;
            if (nucDots != 0) {
                nucDotsVol = findObjectsVolume(nucDotsPop);
                nucDotsInt = findObjectsIntensity(nucDotsPop, imh);
            }
            
            // inner nucleus
            IJ.showStatus("Finding inner nucleus "+nucObj.getLabel()+"/"+nucPop.getNbObjects()+"  parameters ...");
            Object3DInt innerNucObj = innerNucPop.getObjectByLabel(nucObj.getLabel());
            double innerNucVol = new MeasureVolume(innerNucObj).getVolumeUnit();
            double innerNucInt = new MeasureIntensity(innerNucObj, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
            //Finding dots in inner nucleus
            Objects3DIntPopulation innerNucDotsPop = findDotsPop(allDots, innerNucObj);
            int innerNucDots = innerNucDotsPop.getNbObjects();
            double innerNucDotsVol = 0;
            double innerNucDotsInt = 0;
            if (innerNucDots != 0) {
                innerNucDotsVol = findObjectsVolume(innerNucDotsPop);
                innerNucDotsInt = findObjectsIntensity(innerNucDotsPop, imh);
            }
            
            //Finding outer nucleus
            Object3DInt outerNucObj = outerNucPop.getObjectByLabel(nucObj.getLabel());
            double outerNucVol = new MeasureVolume(outerNucObj).getVolumeUnit();
            double outerNucInt = new MeasureIntensity(outerNucObj, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
            IJ.showStatus("Finding dots in outer nucleus "+outerNucObj.getLabel()+"/"+outerNucPop.getNbObjects()+" ...");
            Objects3DIntPopulation outerNucDotsPop = findDotsPop(allDots, outerNucObj);
            int outerNucDots = outerNucDotsPop.getNbObjects();
            double outerNucDotsVol = 0;
            double outerNucDotsInt = 0;
            if (outerNucDots != 0) {
                outerNucDotsVol = findObjectsVolume(outerNucDotsPop);
                outerNucDotsInt = findObjectsIntensity(outerNucDotsPop, imh);
            }
            
            // inner nucleus Ring
            IJ.showStatus("Finding inner ring "+nucObj.getLabel()+"/"+nucPop.getNbObjects()+"  parameters ...");
            Object3DInt innerRingNucObj = innerRingPop.getObjectByLabel(nucObj.getLabel());
            double innerRingNucVol = new MeasureVolume(innerRingNucObj).getVolumeUnit();
            double innerRingNucInt = new MeasureIntensity(innerRingNucObj, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
            IJ.showStatus("Finding dots in innerRing nucleus "+innerRingNucObj.getLabel()+"/"+innerRingPop.getNbObjects()+" ...");
            Objects3DIntPopulation innerRingNucDotsPop = findDotsPop(allDots, innerRingNucObj);
            int innerRingNucDots = innerRingNucDotsPop.getNbObjects();
            double innerRingNucDotsVol = 0;
            double innerRingNucDotsInt = 0;
            if (innerRingNucDots != 0) {
                innerRingNucDotsVol = findObjectsVolume(innerRingNucDotsPop);
                innerRingNucDotsInt = findObjectsIntensity(innerRingNucDotsPop, imh);
            }
            
            // cell cyto
            IJ.showStatus("Finding cell cytoplasm "+nucObj.getLabel()+"/"+nucPop.getNbObjects()+"  parameters ...");
            Object3DInt cellObj = cellsPop.getObjectByLabel(nucObj.getLabel());
            double cytoVol = 0;
            double cytoInt = 0;
            double cytoDotsVol = 0;
            double cytoDotsInt = 0;
            int cytoDots = 0;
            if (cellObj != null) {
                cytoVol = new MeasureVolume(cellObj).getVolumeUnit();
                cytoInt = new MeasureIntensity(cellObj, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
                IJ.showStatus("Finding dots in cytoplasm "+cellObj.getLabel()+"/"+cellsPop.getNbObjects()+" ...");
                Objects3DIntPopulation cytoDotsPop = findDotsPop(allDots, cellObj);
                cytoDots = cytoDotsPop.getNbObjects();
                if (cytoDots != 0) {
                    cytoDotsVol = findObjectsVolume(cytoDotsPop);
                    cytoDotsInt = findObjectsIntensity(cytoDotsPop, imh);
                }
            }
            
            // add cell parameters
            Nucleus nucleus = new Nucleus((int)nucObj.getLabel(), nucVol, nucComp, nucSph, nucElongation, nucFlatness, nucInt, nucDots, nucDotsVol, nucDotsInt, innerNucVol, innerNucInt, innerNucDots, innerNucDotsVol, innerNucDotsInt,
            innerRingNucVol, innerRingNucInt, innerRingNucDots, innerRingNucDotsVol, innerRingNucDotsInt, outerNucVol, outerNucInt, outerNucDots, outerNucDotsVol, outerNucDotsInt,
                    cytoVol, cytoInt, cytoDots, cytoDotsVol, cytoDotsInt);
            nuclei.add(nucleus);
        }
        imh.closeImagePlus();
        return(nuclei);
    }
    
}