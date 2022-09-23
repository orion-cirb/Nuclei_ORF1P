package Tools;


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
import ij.plugin.frame.RoiManager;
import ij.process.ImageProcessor;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Object3D;
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.Object3D_IJUtils;
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

public class nucleus_Tools3D {
    

    public static double minNucVol= 100;
    public static double maxNucVol = 2000;
    public static double minDotVol= 0.1;
    public static double maxDotVol = 10;
    public static float nucDil = 3;
   

    /**
     *
     * @param img
     */
    public static void closeImages(ImagePlus img) {
        img.flush();
        img.close();
    }
    
    
    
    
  /**
     * return objects population in an binary image
     * @param img
     * @return pop objects population
     */

    public static  Objects3DPopulation getPopFromImage(ImagePlus img) {
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
    public static void median_filter(ImagePlus img, double size) {
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
    public static void filterCells(Objects3DPopulation popPV, double sphCoef) {
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
     * @return 
     */
    
    public static ArrayList dialog(List<String> channels, List<String> channelsName, boolean dilate) {
        ArrayList ch = new ArrayList();
        
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.addMessage("Choose channels");
        int index = 0;
        for (String chName : channelsName) {
            gd.addChoice(chName, channels.toArray(new String[0]), channels.get(index));
            index++;
        }
        if (dilate)
            gd.addNumericField("Nucleus dilatation (µm) :", nucDil);
        gd.showDialog();
        for (int i = 0; i < index; i++)
            ch.add(i, gd.getNextChoice());
        if (dilate)
            nucDil = (float)gd.getNextNumber();
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
    
    public static ImagePlus maskImage(ImagePlus img, Objects3DPopulation objPop, String dir, String filename) {
        ImageHandler imh = ImageHandler.wrap(img);
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
    public static ArrayList findImages(String imagesFolder, String imageExt) {
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
            if (fileExt.equals(imageExt))
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
    public static List<String> findChannels (String imageName) throws DependencyException, ServiceException, FormatException, IOException {
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
    public static Calibration findImageCalib(IMetadata meta) {
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
     * Return dilated object restriced to image borders
     * @param img
     * @param obj
     * @return 
     */
    public static Object3DVoxels dilCellObj(ImagePlus img, Object3D obj, double nucDil) {
        Calibration cal = img.getCalibration();
        Object3D objDil = obj.getDilatedObject((float)(nucDil/cal.pixelWidth), (float)(nucDil/cal.pixelHeight), 
                (float)(nucDil/cal.pixelDepth));
        // check if object go outside image
        if (objDil.getXmin() < 0 || objDil.getXmax() > img.getWidth() || objDil.getYmin() < 0 || objDil.getYmax() > img.getHeight()
                || objDil.getZmin() < 0 || objDil.getZmax() > img.getNSlices()) {
            Object3DVoxels voxObj = new Object3DVoxels(objDil.listVoxels(ImageHandler.wrap(img)));
            return(voxObj);
        }
        else
            return(objDil.getObject3DVoxels());
    }
    
    public static Objects3DPopulation findNucleus(ImagePlus imgNuc, String thMethod) {
        Objects3DPopulation nucPopOrg = find_nucleus2(imgNuc, thMethod);
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
    public static Objects3DPopulation find_nucleus2(ImagePlus imgNuc, String thMethod) {
        ImagePlus img = new Duplicator().run(imgNuc);
        ImageStack stack = new ImageStack(img.getWidth(), imgNuc.getHeight());
        for (int i = 1; i <= img.getStackSize(); i++) {
            IJ.showStatus("Finding nucleus section "+i+" / "+img.getStackSize());
            img.setZ(i);
            img.updateAndDraw();
            IJ.run(img, "Nuclei Outline", "blur=18 blur2=20 threshold_method="+thMethod+" outlier_radius=0 outlier_threshold=1 max_nucleus_size=500 "
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
       
//        IJ.showStatus("Starting watershed...");
//        ImagePlus imgWater = WatershedSplit(imgStack, 8);
//        closeImages(imgStack);
//        imgWater.setCalibration(img.getCalibration());
//        imgWater.show();
//        Objects3DPopulation nucPop = new Objects3DPopulation(imgWater);
//        nucPop.removeObjectsTouchingBorders(imgWater, false);
//        closeImages(imgWater);
        Objects3DPopulation nucPop = getPopFromImage(imgStack);
        closeImages(img);
        closeImages(imgStack);
        return(nucPop);
    }
    
    /**
     * Read sum of pixel in image zero value excluded
     * @param img
     * @return 
     */
    public static double readSumIntensity(ImagePlus img) {
        ImagePlus imgProj = doZProjection(img, 3);
        IJ.setThreshold(img, 1, img.getProcessor().getMax());
        ResultsTable rt = new ResultsTable();
        Analyzer ana = new Analyzer(img, Measurements.INTEGRATED_DENSITY + Measurements.LIMIT, rt);
        ana.measure();
        double intValue = rt.getValue("IntDen", 0);
        closeImages(imgProj);
        return(intValue);
    }

    public static ArrayList<Double> readIntensity(ImagePlus img, Objects3DPopulation objPop) {
        ArrayList<Double> intensity = new ArrayList();
        ImageHandler ima = ImageHandler.wrap(img);
        for (int i = 0; i < objPop.getNbObjects(); i++) {
            intensity.add(objPop.getObject(i).getIntegratedDensity(ima));
        }
        return(intensity);
    }
    
    
    private static ImagePlus WatershedSplit(ImagePlus binaryMask, float rad) {
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
    
    public static void saveImageObjects(Objects3DPopulation nucPop, Objects3DPopulation cellPop, ImagePlus img, String name, int fontSize) {
        
        ImageHandler nucImgObj = ImageHandler.wrap(img).createSameDimensions();
        nucImgObj.setCalibration(img.getCalibration());
        ImageHandler cellImgObj = nucImgObj.duplicate();
        // draw obj population
        nucPop.draw(nucImgObj, 64);
        labelsObject(nucPop, nucImgObj.getImagePlus(), fontSize);
        cellPop.draw(cellImgObj, 64);
        ImagePlus[] imgColors = {cellImgObj.getImagePlus(), null, nucImgObj.getImagePlus(), img};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, true);
        imgObjects.setCalibration(img.getCalibration());
        FileSaver ImgObjectsFile = new FileSaver(imgObjects);
        ImgObjectsFile.saveAsTiff(name); 
        nucImgObj.closeImagePlus(); 
        cellImgObj.closeImagePlus();
    }

    
    /**
     * Do Z projection
     * @param img
     * @param projection parameter
     * @return 
     */
    public static ImagePlus doZProjection(ImagePlus img, int param) {
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
    public static double[] find_background(ImagePlus img) {
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
     * Create donut object population
     * 
     */
    public static Objects3DPopulation createDonutPop(Objects3DPopulation pop, ImagePlus img, float dilateStepXY, float dilateStepZ) {
        Calibration cal = img.getCalibration();
        float dil = (float)(dilateStepXY/cal.pixelWidth);
        ImagePlus imgCopy = new Duplicator().run(img);
        ImageInt imgBin = ImageInt.wrap(imgCopy);
        Objects3DPopulation donutPop = new Objects3DPopulation();
        Object3D obj, objDil;
        for (int i = 0; i < pop.getNbObjects(); i++) {
            IJ.showStatus("Creating cell object "+i+"/"+pop.getNbObjects());
            imgBin.fill(0);
            obj = pop.getObject(i);
            objDil = obj.getDilatedObject(dil, dil, dilateStepZ);
            // check if object go outside image
            if (objDil.getXmin() < 0 || objDil.getXmax() > img.getWidth() || objDil.getYmin() < 0 || objDil.getYmax() > img.getHeight()
                    || objDil.getZmin() < 0 || objDil.getZmax() > img.getNSlices()) {
                Object3DVoxels voxObj = new Object3DVoxels(objDil.listVoxels(ImageHandler.wrap(img)));
                voxObj.draw(imgBin, 255);
            }
            else
                objDil.draw(imgBin, 255);
            obj.draw(imgBin, 0);
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
    public static Objects3DPopulation find_dots(ImagePlus img) {
        ImagePlus imgDup = new Duplicator().run(img);
        median_filter(imgDup, 1.5);
        IJ.run(imgDup, "Difference of Gaussians", "  sigma1=2 sigma2=1 enhance stack");
        imgDup.setSlice(imgDup.getNSlices()/2);
        IJ.setAutoThreshold(imgDup, "Triangle dark");
        Prefs.blackBackground = false;
        IJ.run(imgDup, "Convert to Mask", "method=Triangle background=Dark");
        Objects3DPopulation pop = new Objects3DPopulation(getPopFromImage(imgDup).getObjectsWithinVolume(minDotVol, maxDotVol, true));
        closeImages(imgDup);
        return(pop);
    }
                
    
    /**
     * Label object
     * @param popObj
     * @param img 
     */
    public static void labelsObject (Objects3DPopulation popObj, ImagePlus img, int fontSize) {
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
    
   
    
}