package fdr;

public class FalseDiscoveryRate {
    private double FalseDiscoveryRate;
    private double Pvalue;
    private int PassingAnnotations;
    private int nullAnnotations;

    public FalseDiscoveryRate(double falseDiscoveryRate, double pvalue, int passingAnnotations, int nullAnnotat) {
        FalseDiscoveryRate = falseDiscoveryRate;
        Pvalue = pvalue;
        this.PassingAnnotations= passingAnnotations;
        this.nullAnnotations = nullAnnotat;
    }

    public double getFalseDiscoveryRate() {
        return FalseDiscoveryRate;
    }

    public void setFalseDiscoveryRate(double falseDiscoveryRate) {
        FalseDiscoveryRate = falseDiscoveryRate;
    }

    public double getPvalue() {
        return Pvalue;
    }

    public void setPvalue(double pvalue) {
        Pvalue = pvalue;
    }

    public int getPassingAnnotation() {
        return PassingAnnotations;
    }

    public void setPassingAnnotations(int passingAnnotations) {
        this.PassingAnnotations = passingAnnotations;
    }
    
    public int getNullAnnotations() {
    	return nullAnnotations;
    }

}


