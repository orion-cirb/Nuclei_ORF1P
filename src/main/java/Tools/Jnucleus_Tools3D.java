package Tools;

import StardistOrion.StarDist2D;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.Duplicator;
import ij.plugin.GaussianBlur3D;
import ij.plugin.RGBStackMerge;
import ij.plugin.ZProjector;
import ij.plugin.filter.Analyzer;
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
import mcib3d.geom.Object3D;
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Point3D;
import mcib3d.image3d.ImageFloat;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import mcib3d.image3d.distanceMap3d.EDT;
import mcib3d.image3d.processing.FastFilters3D;
import mcib3d.image3d.regionGrowing.Watershed3D;
import org.apache.commons.io.FilenameUtils;



/**
 *
 * @author phm
 */

public class Jnucleus_Tools3D {
    

    public double minNucVol= 100;
    public double maxNucVol = 2000;
    public double minDotVol= 0.05;
    public double maxDotVol = 50;
    public float innerNucDil = 1;
    public float outerNucDil = 2;
    public int zMax = 1;
    public boolean zCrop = false;
    public boolean dotsDetect = true;
    public Calibration cal;
    
    // StarDist
    public Object syncObject = new Object();
    public final double stardistPercentileBottom = 0.2;
    public final double stardistPercentileTop = 99.8;
    public final double stardistProbThreshNuc = 0.65;
    public final double stardistOverlayThreshNuc = 0.01;
    public File modelsPath = new File(IJ.getDirectory("imagej")+File.separator+"models");
    public String stardistModel = "";
    public String stardistOutput = "Label Image"; 
    
    public String nucleusDetector = "";
    
    
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

    public  Objects3DPopulation getPopFromImage(ImagePlus img) {
        // label binary images first
        ImageLabeller labeller = new ImageLabeller();
        ImageInt labels = labeller.getLabels(ImageHandler.wrap(img));
        Objects3DPopulation pop = new Objects3DPopulation(labels);
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
     * Filters cells on sphericity
     */
    public void filterCells(Objects3DPopulation popPV, double sphCoef) {
        for (int i = 0; i < popPV.getNbObjects(); i++) {
            Object3D obj = popPV.getObject(i);
            double sph = obj.getSphericity(true);
            if (sph < sphCoef){
                popPV.removeObject(i);
                i--;
            }
        }
    }
  
    /**
     * Ask for parameters
     * @param channels
     * @param channelsName
     * @param dilate
     * @return 
     */
    
    public ArrayList dialog(List<String> channels, List<String> channelsName, boolean dilate) {
        ArrayList ch = new ArrayList();
        String[] models = findStardistModels();
        String[] nucleusDetectors = {"StarDist", "DOG"}; 
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 80, 0);
        gd.addImage(icon);
        gd.addMessage("Choose channels", Font.getFont("Monospace"), Color.blue);
        int index = 0;
        for (String chName : channelsName) {
            gd.addChoice(chName, channels.toArray(new String[0]), channels.get(index));
            index++;
        }
        gd.addMessage("Nucleus detection method", Font.getFont("Monospace"), Color.blue);
        gd.addChoice("Nucleus segmentation method :",nucleusDetectors, nucleusDetectors[0]);
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
        if (dilate) {
            gd.addMessage("Nucleus dounuts", Font.getFont("Monospace"), Color.blue);
            gd.addNumericField("Nucleus cyto ring (µm) :", outerNucDil);
            gd.addNumericField("Nucleus inner ring (µm) :", innerNucDil);
        }
        else {
            gd.addNumericField("Min dots vol. :", minDotVol);
            gd.addNumericField("Max dots vol. :", maxDotVol);
        }
        gd.addCheckbox("  Do zCrop", zCrop);
        gd.addCheckbox("  Do dots detection", dotsDetect);
        gd.showDialog();
        for (int i = 0; i < index; i++)
            ch.add(i, gd.getNextChoice());
        
        nucleusDetector = gd.getNextChoice();
        if (models.length > 0) {
            stardistModel = modelsPath+File.separator+gd.getNextChoice();
        }
        else {
            stardistModel = gd.getNextString();
        }
        if (nucleusDetector.equals("StarDist") && stardistModel.isEmpty()) {
            IJ.error("No model specify !!");
            return(null);
        }
        minNucVol = (float)gd.getNextNumber();
        maxNucVol = (float)gd.getNextNumber();
        if (dilate) {
            outerNucDil = (float)gd.getNextNumber();
            innerNucDil = (float)gd.getNextNumber();
        }
        else {
            minDotVol = (float)gd.getNextNumber();
            maxDotVol = (float)gd.getNextNumber();
        }
        zCrop = gd.getNextBoolean();
        dotsDetect = gd.getNextBoolean();
        if(gd.wasCanceled())
            ch = null;
        return(ch);
    }
    
    /**
     * Mask image draw obj with 0 
     * @param img
     * @param objPop
     * @return 
     */
    
    public ImagePlus maskImage(ImagePlus img, Objects3DPopulation objPop, String dir, String filename) {
        ImageHandler imh = ImageHandler.wrap(img).duplicate();
        // draw obj with 0
        for (int i = 0; i < objPop.getNbObjects(); i++) {
            Object3D obj = objPop.getObject(i);
            obj.draw(imh, 0);
        }
        
        // Save masked image
        imh.setTitle(filename);
        imh.save(dir);
        return(imh.getImagePlus());
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
        Calibration cal = new Calibration();  
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
     * Crop stack from max Z intensity to end stack
     */
    public ImagePlus cropZmax(ImagePlus img) {
        double meanInt = 0;
        for (int z = 1; z <= img.getNSlices(); z++) {
            img.setSlice(z);
            ImageProcessor imp = img.getProcessor();
            double mean = imp.getStatistics().mean;
            if (mean >= meanInt) {
                meanInt = mean;
                zMax = z;
            }
        }
        if (cal.pixelWidth < 0.103)
            zMax = (zMax == 1) ? zMax : zMax - 1;
        ImagePlus imgCroped = new Duplicator().run(img, zMax, img.getNSlices());
        return(imgCroped);
    }
    
    
    /**
     * Return dilated object restriced to image borders
     * @param img
     * @param obj
     * @return 
     */
    public Object3DVoxels dilCellObj(ImagePlus img, Object3D obj, double nucDil, boolean dil) {
        Calibration cal = img.getCalibration();
        float dilXY = (float)(nucDil/cal.pixelWidth);
        float dilZ = (float)(nucDil/cal.pixelHeight);
        Object3D objDil = null;
        if (dil) {
            objDil = obj.getDilatedObject(dilXY, dilXY, dilZ);
            // check if object go outside image
            if (objDil.getXmin() < 0 || objDil.getXmax() > img.getWidth() || objDil.getYmin() < 0 || objDil.getYmax() > img.getHeight()
                    || objDil.getZmin() < 0 || objDil.getZmax() > img.getNSlices()) {
                Object3DVoxels voxObj = new Object3DVoxels(objDil.listVoxels(ImageHandler.wrap(img)));
                return(voxObj);
            }
        }
        else
            objDil = obj.getErodedObject(dilXY, dilXY, 0);
        return(objDil.getObject3DVoxels());
    }
    
    /*
    Find cell cytoplasm
    */
    public Objects3DPopulation findCells (ImagePlus img, Objects3DPopulation nucPop) {
        Objects3DPopulation cellPop = new Objects3DPopulation();
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
        for (int i = 0; i < nucPop.getNbObjects(); i++) {
            Object3D nucObj = nucPop.getObject(i);
            Point3D pt = nucObj.getCenterAsPoint();
            imgCell.setSlice(pt.getRoundZ());
            PointRoi roi = new PointRoi(pt.getRoundX(), pt.getRoundY());
            imgCell.setRoi(roi, true);
            cell.setup("", imgCell);
            cell.run(imgCell.getProcessor());
            ImagePlus imgMask = cell.maskImp;
            imgMask.setCalibration(cal);
            nucObj.draw(imgMask.getStack(), 0);
            // get the cell cytoplasm
            Object3D cellCyto = getPopFromImage(imgMask).getObject(0);
            closeImages(imgMask);
            cellPop.addObject(cellCyto);
        }
        closeImages(imgCell);
        System.out.println(cellPop.getNbObjects()+" cells found");
        return(cellPop);
    }
    
    
    public Objects3DPopulation findNucleus(ImagePlus imgNuc, int blur1, int blur2, int radOut, String thMethod) {
        Objects3DPopulation nucPopOrg = find_nucleus2(imgNuc, blur1, blur2, radOut, thMethod);
        System.out.println("-- Total nucleus Population :"+nucPopOrg.getNbObjects());
        // size filter
        Objects3DPopulation nucPop = new Objects3DPopulation(nucPopOrg.getObjectsWithinVolume(minNucVol, maxNucVol, true));
        int nbNucPop = nucPop.getNbObjects();
        System.out.println("-- Total nucleus Population after size filter: "+ nbNucPop);
        return(nucPop);
    }
        
    /**
     * Nucleus segmentation 2
     * @param imgNuc
     * @return cellPop
     */
    public Objects3DPopulation find_nucleus2(ImagePlus imgNuc, int blur1, int blur2, int radOut, String thMethod) {
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
        imgStack.setCalibration(imgNuc.getCalibration());
        IJ.run(imgStack, "Watershed", "stack");
        Objects3DPopulation nucPop = getPopFromImage(imgStack);
        closeImages(img);
        closeImages(imgStack);
        return(nucPop);
    }
    
    
    /** Look for all nuclei
         Do z slice by slice stardist 
         * return nuclei population
         */
        public Objects3DPopulation stardistNucleiPop(ImagePlus imgNuc) throws IOException{
            // resize to be in a stardist-friendly scale
            ImagePlus img = null;
            int width = imgNuc.getWidth();
            int height = imgNuc.getHeight();
            float factor = 0.25f;
            boolean resized = false;
            if (imgNuc.getWidth() > 512) {
                img = imgNuc.resize((int)(width*factor), (int)(height*factor), 1, "none");
                resized = true;
            }
            else
                img = new Duplicator().run(imgNuc);int newWidth = width/2;
            
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
            Objects3DPopulation nucPop = new Objects3DPopulation(label3D);
            Objects3DPopulation nPop = new Objects3DPopulation(nucPop.getObjectsWithinVolume(minNucVol, maxNucVol, true));
            closeImages(nuclei);
            return(nPop);
        }
    
    /**
     * Read sum of pixel cell extensions
     * @param img
     * @return 
     */
    public double readSumIntensity(ImagePlus img, double bg, String name) {
        // Subtract background
        if (bg != 0) {
            median_filter(img, 2);
            IJ.run(img, "Subtract...", "value="+bg+" stack");
        }
        ImagePlus imgProj = doZProjection(img, 3);
        IJ.setAutoThreshold(imgProj, "Percentile dark");
        Prefs.blackBackground = false;
        ResultsTable rt = new ResultsTable();
        Analyzer ana = new Analyzer(imgProj, Measurements.INTEGRATED_DENSITY + Measurements.LIMIT, rt);
        ana.measure();
        double intValue = rt.getValue("IntDen", 0);
        
        if (name != null) {
            IJ.run(imgProj, "Convert to Mask", "method=Percentile background=Dark");
            FileSaver ImgObjectsFile = new FileSaver(imgProj);
            ImgObjectsFile.saveAsTiff(name);
        }
        closeImages(imgProj);
        return(intValue);
    }

    public ArrayList<Double> readIntensity(ImagePlus img, Objects3DPopulation objPop) {
        ArrayList<Double> intensity = new ArrayList();
        ImageHandler ima = ImageHandler.wrap(img);
        for (int i = 0; i < objPop.getNbObjects(); i++) {
            intensity.add(objPop.getObject(i).getIntegratedDensity(ima));
        }
        return(intensity);
    }
    
    
    private ImagePlus WatershedSplit(ImagePlus binaryMask, float rad) {
        float resXY = 1;
        float resZ = 1;
        float radXY = rad;
        float radZ = rad;
        Calibration cal = binaryMask.getCalibration();
        if (cal != null) {
            resXY = (float) cal.pixelWidth;
            resZ = (float) cal.pixelDepth;
            radZ = radXY * (resXY / resZ);
        }
        ImageInt imgMask = ImageInt.wrap(binaryMask);
        ImageFloat edt = EDT.run(imgMask, 0, resXY, resZ, false, 0);
        ImageHandler edt16 = edt.convertToShort(true);
        ImagePlus edt16Plus = edt16.getImagePlus();
        GaussianBlur3D.blur(edt16Plus, 4.0, 4.0, 4.0);
        edt16 = ImageInt.wrap(edt16Plus);
        edt16.intersectMask(imgMask);
        // seeds
        ImageHandler seedsImg = FastFilters3D.filterImage(edt16, FastFilters3D.MAXLOCAL, radXY, radXY, radZ, 0, false);
        Watershed3D water = new Watershed3D(edt16, seedsImg, 0, 0);
        water.setLabelSeeds(true);
        return(water.getWatershedImage3D().getImagePlus());
    }

   
    
   /**
    * Save image objects
    * nucPop blue cellPop red
    * @param nucPop
    * @param cellPop
    * @param imgCells
    * @param name 
    */
    
    public void saveImageObjects(Objects3DPopulation pop1, Objects3DPopulation pop2, Objects3DPopulation pop3, ImagePlus img, String name, int fontSize) {
        
        ImageHandler imgObj1 = ImageHandler.wrap(img).createSameDimensions();
        imgObj1.setCalibration(img.getCalibration());
        ImageHandler imgObj2 = imgObj1.duplicate();
        ImageHandler imgObj3 = imgObj1.duplicate();
        if (pop3 != null)
            pop3.draw(imgObj3, 255);
        // draw obj population
        pop1.draw(imgObj1, 255);
        labelsObject(pop1, imgObj1.getImagePlus(), fontSize);
        pop2.draw(imgObj2, 255);
        ImagePlus[] imgColors = {imgObj2.getImagePlus(), imgObj3.getImagePlus(), imgObj1.getImagePlus(), img};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, true);
        imgObjects.setCalibration(img.getCalibration());
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
    public double[] find_background(ImagePlus img) {
      double[] bg = new double[2];
      ImagePlus imgProj = doZProjection(img, ZProjector.MIN_METHOD);
      ImageProcessor imp = imgProj.getProcessor();
      bg[0] = imp.getStatistics().mean;
      bg[1] = imp.getStatistics().stdDev;
      System.out.println("Background =  " + bg[0] + "+-" + bg[1]);
      closeImages(imgProj);
      return(bg);
    }
    
    /**
     * Get inner nucleus
     * 
     */
    public Objects3DPopulation getInnerNucleus(Objects3DPopulation pop, ImagePlus img, float ringXY, boolean dil) {
        Objects3DPopulation innerNucleusPop = new Objects3DPopulation();
        for (int i = 0; i < pop.getNbObjects(); i++) {
            Object3D obj = pop.getObject(i);
            Object3D objDil = dilCellObj(img, obj, ringXY, false);
            innerNucleusPop.addObject(objDil);
        }
        innerNucleusPop.setCalibration(cal);
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
    public Objects3DPopulation createDonutPop(Objects3DPopulation pop, ImagePlus img, float ringXY, boolean dil) {
        Calibration cal = img.getCalibration();
        float dilXY = (float)(ringXY/cal.pixelWidth);
        float dilZ = (float)(ringXY / (cal.pixelDepth/cal.pixelWidth));
        ImagePlus imgCopy = new Duplicator().run(img);
        ImageInt imgBin = ImageInt.wrap(imgCopy);
        Objects3DPopulation donutPop = new Objects3DPopulation();
        for (int i = 0; i < pop.getNbObjects(); i++) {
            imgBin.fill(0);
            Object3D obj = pop.getObject(i);
            Object3D objDil = dilCellObj(img, obj, ringXY, dil);
            if (dil) {
                objDil.draw(imgBin, 255);
                obj.draw(imgBin, 0);
            }
            else {
                obj.draw(imgBin, 255);
                objDil.draw(imgBin, 0);
            }
            Objects3DPopulation tmpPop = getPopFromImage(imgBin.getImagePlus());
            donutPop.addObject(tmpPop.getObject(0));
        }
        closeImages(imgCopy);
        imgBin.closeImagePlus();
        return(donutPop);
    }
    
   
    
    /**
     * Find dots
     */
    public Objects3DPopulation find_dots(ImagePlus img, int sig1, int sig2, String th) {
        ImagePlus imgDup = new Duplicator().run(img);
        //median_filter(imgDup, 1.5);
        IJ.run(imgDup, "Difference of Gaussians", "  sigma1="+sig1+" sigma2=" +sig2+" enhance stack");
        imgDup.setSlice(imgDup.getNSlices()/2);
        IJ.setAutoThreshold(imgDup, th+" dark");
        Prefs.blackBackground = false;
        IJ.run(imgDup, "Convert to Mask", "method="+th+" background=Dark");
        Objects3DPopulation pop = new Objects3DPopulation(getPopFromImage(imgDup).getObjectsWithinVolume(minDotVol, maxDotVol, true));
        closeImages(imgDup);
        return(pop);
    }
    
    /**
     * Find dots in obj bounding box
     */
    public Objects3DPopulation findDotsPop(ImagePlus img, ImagePlus imgDots, Object3D obj) {
        Objects3DPopulation popDotsIn = new Objects3DPopulation();
        Roi roi = new Roi(obj.getXmin(), obj.getYmin(), obj.getXmax() - obj.getXmin() + 1, obj.getYmax() - obj.getYmin() + 1);
        // crop img dots
        imgDots.setRoi(roi, true);
        ImagePlus imgDotsObj = new Duplicator().run(imgDots, 1, img.getNSlices());
        Objects3DPopulation popDots = getPopFromImage(imgDotsObj);
        closeImages(imgDotsObj);
        System.out.println(popDots.getNbObjects()+" dots");
        img.setRoi(roi, true);
        ImagePlus imgObj = new Duplicator().run(img, 1, img.getNSlices());
        // crop img obj
        obj.setNewCenter(imgObj.getWidth()/2, imgObj.getHeight()/2, obj.getCenterZ());
        obj.draw(imgObj.getStack(), 255);
        closeImages(imgObj);
        for (int i = 0; i < popDots.getNbObjects(); i++) {
            Object3D obj2 = popDots.getObject(i);
            if (obj2.hasOneVoxelColoc(obj))
                popDotsIn.addObject(obj2);
        }
        return(popDotsIn);
    }
    
    /*
     * Find dots Volume
    */
    public double findDotsVolume(Objects3DPopulation pop) {
        double vol = 0;
        for (int i = 0; i <  pop.getNbObjects(); i++) {
            Object3D obj = pop.getObject(i);
            vol += obj.getVolumeUnit();
        }
        return(vol);
    }
    
    /*
     * Find dots intensity
    */
    public double findDotsIntensity(Objects3DPopulation pop, ImageHandler imh) {
        double sumInt = 0;
        for (int i = 0; i <  pop.getNbObjects(); i++) {
            Object3D obj = pop.getObject(i);
            sumInt += obj.getIntegratedDensity(imh);
        }
        
        return(sumInt);
    }
    
    
    /**
     * Label object
     * @param popObj
     * @param img 
     */
    public void labelsObject (Objects3DPopulation popObj, ImagePlus img, int fontSize) {
        if (IJ.isMacOSX())
            fontSize *= 3;
        Font tagFont = new Font("SansSerif", Font.PLAIN, fontSize);
        for (int n = 0; n < popObj.getNbObjects(); n++) {
            Object3D obj = popObj.getObject(n);
            int[] box = obj.getBoundingBox();
            int z = (int)obj.getCenterZ();
            int x = box[0] - 2;
            int y = box[2] - 2;
            img.setSlice(z+1);
            ImageProcessor ip = img.getProcessor();
            ip.setFont(tagFont);
            ip.setColor(255);
            ip.drawString(Integer.toString(n), x, y);
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
    
    public ArrayList<Nucleus> tagsNuclei(ImagePlus img, Objects3DPopulation nucPop, Objects3DPopulation innerNucPop, Objects3DPopulation innerRingPop,
            Objects3DPopulation outerPop, Objects3DPopulation cytoPop, Objects3DPopulation allDots) {
        ImagePlus imgDots = IJ.createImage("Dots",img.getWidth(), img.getHeight(),img.getNSlices(), 8);
        imgDots.setCalibration(cal);
        if (allDots.getNbObjects() > 0)
            allDots.draw(imgDots.getStack(), 255);
        ArrayList<Nucleus> nuclei = new ArrayList<>();
        ImageHandler imh = ImageHandler.wrap(img);
        int nucs = nucPop.getNbObjects();
        for (int i = 0; i < nucs; i++) {
            // nucleus
            Object3D nucObj = nucPop.getObject(i);
            double nucVol = nucObj.getVolumeUnit();
            double nucInt = nucObj.getIntegratedDensity(imh);
            System.out.println("Finding nucleus "+i+"/"+nucs+" dots ...");
            Objects3DPopulation nucDotsPop = findDotsPop(img, imgDots, nucObj);
            int nucDots = nucDotsPop.getNbObjects();
            double nucDotsVol = 0;
            double nucDotsInt = 0;
            if (nucDots != 0) {
                nucDotsVol = findDotsVolume(nucDotsPop);
                nucDotsInt = findDotsIntensity(nucDotsPop, imh);
            }
            
            // inner nucleus
            IJ.showStatus("Finding inner nucleus "+i+"/"+nucs+"  parameters ...");
            Object3D innerNucObj = innerNucPop.getObject(i);
            double innerNucVol = innerNucObj.getVolumeUnit();
            double innerNucInt = innerNucObj.getIntegratedDensity(imh);
            Objects3DPopulation innerNucDotsPop = findDotsPop(img, imgDots, innerNucObj);
            int innerNucDots = innerNucDotsPop.getNbObjects();
            double innerNucDotsVol = 0;
            double innerNucDotsInt = 0;
            if (innerNucDots != 0) {
                innerNucDotsVol = findDotsVolume(innerNucDotsPop);
                innerNucDotsInt = findDotsIntensity(innerNucDotsPop, imh);
            }
            
            // outer nucleus Ring
            IJ.showStatus("Finding outer ring  "+i+"/"+nucs+" parameters ...");
            Object3D outerObj = outerPop.getObject(i);
            double outerVol = outerObj.getVolumeUnit();
            double outerInt = outerObj.getIntegratedDensity(imh);
            Objects3DPopulation outerDotsPop = findDotsPop(img, imgDots, outerObj);
            int outerDots = outerDotsPop.getNbObjects();
            double outerDotsVol = 0;
            double outerDotsInt = 0;
            if (outerDots != 0) {
                outerDotsVol = findDotsVolume(outerDotsPop);
                outerDotsInt = findDotsIntensity(outerDotsPop, imh);
            }
            // inner nucleus Ring
            IJ.showStatus("Finding inner ring  "+i+"/"+nucs+" parameters ...");
            Object3D innerRingObj = innerRingPop.getObject(i);
            double innerRingVol = innerRingObj.getVolumeUnit();
            double innerRingInt = innerRingObj.getIntegratedDensity(imh);
            Objects3DPopulation innerRingDotsPop = findDotsPop(img, imgDots, innerRingObj);
            int innerRingDots = innerRingDotsPop.getNbObjects();
            double innerRingDotsVol = 0;
            double innerRingDotsInt = 0;
            if (innerRingDots != 0) {
                innerRingDotsVol = findDotsVolume(innerRingDotsPop);
                innerRingDotsInt = findDotsIntensity(innerRingDotsPop, imh);
            }
            // Cell cytoplasm
            IJ.showStatus("Finding cell cytoplasm  "+i+"/"+nucs+" parameters ...");
            Object3D cytoObj = cytoPop.getObject(i);
            double cytoVol = cytoObj.getVolumeUnit();
            double cytoInt = cytoObj.getIntegratedDensity(imh);
            Objects3DPopulation cytoDotsPop = findDotsPop(img, imgDots, cytoObj);
            int cytoDots = cytoDotsPop.getNbObjects();
            double cytoDotsVol = 0;
            double cytoDotsInt = 0;
            if (cytoDots != 0) {
                cytoDotsVol = findDotsVolume(cytoDotsPop);
                cytoDotsInt = findDotsIntensity(cytoDotsPop, imh);
            }
            // add cell parameters
            Nucleus nucleus = new Nucleus(i, nucVol, nucInt, nucDots, nucDotsVol, nucDotsInt, innerNucVol, innerNucInt, innerNucDots, innerNucDotsVol, innerNucDotsInt,
            innerRingVol, innerRingInt, innerRingDots, innerRingDotsVol, innerRingDotsInt, outerVol, outerInt, outerDots, outerDotsVol, outerDotsInt,
                    cytoVol, cytoInt, cytoDots, cytoDotsVol, cytoDotsInt);
            nuclei.add(nucleus);
        }
        closeImages(imgDots);
        return(nuclei);
    }
    
}