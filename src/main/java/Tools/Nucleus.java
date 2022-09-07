/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Tools;

/**
 *
 * @author phm
 */
public class Nucleus {
    
    // Nucleus index
    private int index;
    // Nucleus volume
    private double nucVol;
    // Nucleus intensity
    private double nucInt;
    // Nucleus dots number
    private int nucDots;
    // Nucleus dots volume
    private double nucDotsVol;
    // Nucleus dots intensity
    private double nucDotsInt;
    // nucleus circularity
    private double nucCir;
    
    // Inner nucleus volume
    private double innerNucVol;
    // Inner nucleus intensity
    private double innerNucInt;
    // Inner nucleus dots number
    private int innerNucDots;
    // Inner nucleus dots volume
    private double innerNucDotsVol;
    // Inner nucleus dots intensity
    private double innerNucDotsInt;
    
    // Inner nucleus ring volume
    private double innerRingVol;
    // Inner nucleus intensity
    private double innerRingInt;
    // Inner nucleus dots number
    private int innerRingDots;
    // Inner nucleus dots volume
    private double innerRingDotsVol;
    // Inner nucleus dots intensity
    private double innerRingDotsInt;
    
    // Outer nucleus ring volume
    private double outerRingVol;
    // Inner nucleus intensity
    private double outerRingInt;
    // Inner nucleus dots number
    private int outerRingDots;
    // Inner nucleus dots volume
    private double outerRingDotsVol;
    // Inner nucleus dots intensity
    private double outerRingDotsInt;
    
    // Cell cytoplasm
    private double cytoVol;
    // Inner nucleus intensity
    private double cytoInt;
    // Inner nucleus dots number
    private int cytoDots;
    // Inner nucleus dots volume
    private double cytoDotsVol;
    // Inner nucleus dots intensity
    private double cytoDotsInt;

  
	
	public Nucleus(int index, double nucVol, double nucCir, double nucInt, int nucDots, double nucDotsVol, double nucDotsInt, double innerNucVol, double innerNucInt, int innerNucDots, double innerNucDotsVol,
                    double innerNucDotsInt, double innerRingVol, double innerRingInt, int innerRingDots, double innerRingDotsVol, double innerRingDotsInt, double outerRingVol, double outerRingInt,
                    int outerRingDots, double outerRingDotsVol, double outerRingDotsInt, double cytoVol, double cytoInt, int cytoDots, double cytoDotsVol, double cytoDotsInt) {
            this.index = index;
            this.nucVol = nucVol;
            this.nucCir =nucCir; 
            this.nucInt = nucInt;
            this.nucDots = nucDots;
            this.nucDotsVol = nucDotsVol;
            this.nucDotsInt = nucDotsInt;
            this.innerNucVol = innerNucVol;
            this.innerNucInt = innerNucInt;
            this.innerNucDots = innerNucDots;
            this.innerNucDotsVol = innerNucDotsVol;
            this.innerNucDotsInt = innerNucDotsInt;
            this.innerRingVol = innerRingVol;
            this.innerRingInt = innerRingInt;
            this.innerRingDots = innerRingDots;
            this.innerRingDotsVol = innerRingDotsVol;
            this.innerRingDotsInt = innerRingDotsInt;
            this.outerRingVol = outerRingVol;
            this.outerRingInt = outerRingInt;
            this.outerRingDots = outerRingDots;
            this.outerRingDotsVol = outerRingDotsVol;
            this.outerRingDotsInt = outerRingDotsInt;
            this.cytoVol = cytoVol;
            this.cytoInt = cytoInt;
            this.cytoDots = cytoDots;
            this.cytoDotsVol = cytoDotsVol;
            this.cytoDotsInt = cytoDotsInt;
	}
        
        public void setIndex(int index) {
            this.index = index;
	}
        
        public void setNucVol(double nucVol) {
            this.nucVol = nucVol;
	}
        
        public void setNucCir(double nucCir) {
            this.nucCir = nucCir;
	}
        
        public void setNucInt(double nucInt) {
            this.nucInt = nucInt;
	}
        
        public void setNucDots(int nucDots) {
            this.nucDots = nucDots;
	}
        
        public void setNucDotsVol(double nucDotsVol) {
            this.nucDotsVol = nucDotsVol;
        }
        
        public void setNucDotsInt(double nucDotsInt) {
            this.nucDotsInt = nucDotsInt;
        }
        
        public void setInnerNucVol(double innerNucVol) {
            this.innerNucVol = innerNucVol;
	}
        
        public void setInnerNucInt(double innerNucInt) {
            this.innerNucInt = innerNucInt;
	}
        
        public void setInnerNucDots(int innerNucDots) {
            this.innerNucDots = innerNucDots;
	}
        
        public void setInnerNucDotsVol(double innerNucDotsVol) {
            this.innerNucDotsVol = innerNucDotsVol;
        }
        
        public void setInnerNucDotsInt(double innerNucDotsInt) {
            this.innerNucDotsInt = innerNucDotsInt;
        }
        
        public void setInnerRingVol(double innerRingVol) {
            this.innerRingVol = innerRingVol;
	}
        
        public void setInnerRingInt(double innerRingInt) {
            this.innerRingInt = innerRingInt;
	}
        
        public void setInnerRingDots(int innerRingDots) {
            this.innerRingDots = innerRingDots;
	}
        
        public void setInnerRingDotsVol(double innerRingDotsVol) {
            this.innerRingDotsVol = innerRingDotsVol;
        }
        
        public void setInnerRingDotsInt(double innerRingDotsInt) {
            this.innerRingDotsInt = innerRingDotsInt;
        }
        
        public void setOuterRingVol(double outerRingVol) {
            this.outerRingVol = outerRingVol;
	}
        
        public void setOuterRingInt(double outerRingInt) {
            this.outerRingInt = outerRingInt;
	}
        
        public void setOuterRingDots(int outerRingDots) {
            this.outerRingDots = outerRingDots;
	}
        
        public void setOuterRingDotsVol(double outerRingDotsVol) {
            this.outerRingDotsVol = outerRingDotsVol;
        }
        
        public void setOuterRingDotsInt(double outerRingDotsInt) {
            this.outerRingDotsInt = outerRingDotsInt;
        }
        
        public void setCytoVol(double cytoVol) {
            this.cytoVol = cytoVol;
	}
        
        public void setCytoInt(double cytoInt) {
            this.cytoInt = cytoInt;
	}
        
        public void setCytoDots(int cytoDots) {
            this.cytoDots = cytoDots;
	}
        
        public void setCytoDotsVol(double cytoDotsVol) {
            this.cytoDotsVol = cytoDotsVol;
        }
        
        public void setCytoDotsInt(double cytoDotsInt) {
            this.cytoDotsInt = cytoDotsInt;
        }
        
        public int getIndex() {
            return index;
        }
        
        public double getNucVol() {
            return nucVol;
        }
                
        public double getNucCir() {
            return nucCir;
        }
        
         public double getNucInt() {
            return nucInt;
	}
        
        public int getNucDots() {
            return nucDots;
	}
        
        public double getNucDotsVol() {
            return nucDotsVol;
        }
        
        public double getNucDotsInt() {
            return nucDotsInt;
        }
        
        public double getInnerNucVol() {
            return innerNucVol;
        }
                
        public double getInnerNucInt() {
            return innerNucInt;
	}
        
        public int getInnerNucDots() {
            return innerNucDots;
	}
        
        public double getInnerNucDotsVol() {
            return innerNucDotsVol;
        }
        
        public double getInnerNucDotsInt() {
            return innerNucDotsInt;
        }
        
        public double getInnerRingVol() {
            return innerRingVol;
        }
                
        public double getInnerRingInt() {
            return innerRingInt;
	}
        
        public int getInnerRingDots() {
            return innerRingDots;
	}
        
        public double getInnerRingDotsVol() {
            return innerRingDotsVol;
        }
        
        public double getInnerRingDotsInt() {
            return innerRingDotsInt;
        }
        
        public double getOuterRingVol() {
            return outerRingVol;
        }
                
        public double getOuterRingInt() {
            return outerRingInt;
	}
        
        public int getOuterRingDots() {
            return outerRingDots;
	}
        
        public double getOuterRingDotsVol() {
            return outerRingDotsVol;
        }
        
        public double getOuterRingDotsInt() {
            return outerRingDotsInt;
        }
        
        public double getCytoVol() {
            return cytoVol;
        }
                
        public double getCytoInt() {
            return cytoInt;
	}
        
        public int getCytoDots() {
            return cytoDots;
	}
        
        public double getCytoDotsVol() {
            return cytoDotsVol;
        }
        
        public double getCytoDotsInt() {
            return cytoDotsInt;
        }
}
