package Nuclei_ORF1P_Tools;

import java.util.HashMap;
import mcib3d.geom2.Object3DInt;

/**
 * @author hm
 */
public class Cell {
    
    public Object3DInt cell;
    public Object3DInt nucleus;
    public Object3DInt cytoplasm;
    public Object3DInt innerRing;
    public Object3DInt outerRing;
    public Object3DInt innerNucleus;
    public HashMap<String, Double> params;
    
    public Cell(Object3DInt cell, Object3DInt nucleus, Object3DInt cytoplasm) {
        this.cell = cell;
        this.nucleus = nucleus;
        this.cytoplasm = cytoplasm;
        this.params = new HashMap<>();
    }
    
    public void setInnerRing(Object3DInt innerRing) {
        this.innerRing = innerRing;
    }
    
    public void setOuterRing(Object3DInt outerRing) {
        this.outerRing = outerRing;
    }
    
    public void setInnerNucleus(Object3DInt innerNucleus) {
        this.innerNucleus = innerNucleus;
    }
    
    public void setParams(double index, double nucVol, double nucComp, double nucSph, double nucEllElong, double nucEllFlat, double nucInt, 
                double nucDots, double nucDotsVol, double nucDotsInt, double innerNucVol, double innerNucInt, double innerNucDots, 
                double innerNucDotsVol, double innerNucDotsInt, double innerRingVol, double innerRingInt, double innerRingDots, 
                double innerRingDotsVol, double innerRingDotsInt, double outerRingVol, double outerRingInt, double outerRingDots, 
                double outerRingDotsVol, double outerRingDotsInt, double cytoVol, double cytoInt, double cytoDots, double cytoDotsVol, 
                double cytoDotsInt) {
        params.put("index", index);
        
        // Nucleus
        params.put("nucVol", nucVol);
        params.put("nucComp", nucComp);
        params.put("nucSph", nucSph);
        params.put("nucEllElong", nucEllElong);
        params.put("nucEllFlat", nucEllFlat);
        params.put("nucInt", nucInt);
        params.put("nucDots", nucDots);
        params.put("nucDotsVol", nucDotsVol);
        params.put("nucDotsInt", nucDotsInt);
        
        // Inner nucleus
        params.put("innerNucVol", innerNucVol);
        params.put("innerNucInt", innerNucInt);
        params.put("innerNucDots", innerNucDots);
        params.put("innerNucDotsVol", innerNucDotsVol);
        params.put("innerNucDotsInt", innerNucDotsInt);
        
        // Inner ring
        params.put("innerRingVol", innerRingVol);
        params.put("innerRingInt", innerRingInt);
        params.put("innerRingDots", innerRingDots);
        params.put("innerRingDotsVol", innerRingDotsVol);
        params.put("innerRingDotsInt", innerRingDotsInt);
        
        // Outer ring
        params.put("outerRingVol", outerRingVol);
        params.put("outerRingInt", outerRingInt);
        params.put("outerRingDots", outerRingDots);
        params.put("outerRingDotsVol", outerRingDotsVol);
        params.put("outerRingDotsInt", outerRingDotsInt);
        
        // Cytoplasm
        params.put("cytoVol", cytoVol);
        params.put("cytoInt", cytoInt);
        params.put("cytoDots", cytoDots);
        params.put("cytoDotsVol", cytoDotsVol);
        params.put("cytoDotsInt", cytoDotsInt);        
    }
}
