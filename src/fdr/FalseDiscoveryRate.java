package fdr;

public class FalseDiscoveryRate {
    private double FalseDiscoveryRate;
    private double Pvalue;
    private int PassingAnnotations;

    public FalseDiscoveryRate(double falseDiscoveryRate, double pvalue, int passingAnnotations) {
        FalseDiscoveryRate = falseDiscoveryRate;
        Pvalue = pvalue;
        this.PassingAnnotations= passingAnnotations;
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
}
