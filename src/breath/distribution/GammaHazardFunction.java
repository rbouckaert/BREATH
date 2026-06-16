package breath.distribution;

import org.apache.commons.math.distribution.GammaDistributionImpl;

import java.text.DecimalFormat;
import java.util.Arrays;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.GammaDistribution;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Function.Constant;
import beast.base.core.Input.Validate;
import beast.base.util.Randomizer;

@Description("Defines a hazard function based on a gamma distribution")
public class GammaHazardFunction extends HazardFunction {
    final public Input<Function> scaleInput = new Input<>("scale", "scale of gamma hazard function");
    final public Input<Function> shapeInput = new Input<>("shape", "shape of gamma hazard function");
    final public Input<Function> rate_trInput = new Input<>("rate", "rate of gamma hazard function", Validate.XOR, scaleInput);
    final public Input<Function> constantInput = new Input<>("C", "constant of tranmission process", new Constant("1.0"));
    final public Input<Boolean> approxInput = new Input<>("approx", "approximate cumulative gamma distribution (faster, but with error < 0.00011)", false);

    private GammaDistribution samplingDist = new GammaDistributionImpl(1, 1);
    
    private Function scale;
	private Function shape;
	private Function rate;
	private Function constant;
	private boolean needsupdate;
	private double [] x;
	private double [] y;
	private boolean approx;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();
		scale = scaleInput.get();
		shape = shapeInput.get();
		rate = rate_trInput.get();
		constant = constantInput.get();
		approx = approxInput.get();
		if (approx) {
			x = new double[1000];
			y = new double[1000];
		}
		update();
	}
	
	DecimalFormat f = new DecimalFormat("#.####");
	DecimalFormat f4 = new DecimalFormat("#.######");
	
	@Override
	public double logS(double t, double d) {
		if (needsupdate) {
			update();
		}
		double logS = 0;
		try {
			logS = -constant.getArrayValue() * cumulativeProbability(t - d);
		} catch (MathException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
//System.err.println("logS" + getID() + "(" + f.format(t) + "-" + f.format(d) + "=" + f.format(t-d) + ") \t= " + f4.format(logS));		
//System.err.println("logS" + getID() + "(" + f.format(t-d) + ") # = " + f4.format(logS));		
		return logS;
	}

	private double cumulativeProbability(double d) throws MathException {
		if (approx) {
			int i = Arrays.binarySearch(x, d);
			if (i >= 0) {
				return y[i];
			}
			i = -2-i;
			double result = y[i] + (d-x[i])/(x[i+1]-x[i]) * (y[i+1] - y[i]);
			return result;
		} else {
			return samplingDist.cumulativeProbability(d);
		}
	}

	private void update() {
		samplingDist.setAlpha(shape.getArrayValue());
		if (rate != null) {
			samplingDist.setBeta(1.0/rate.getArrayValue());
		} else {
			samplingDist.setBeta(scale.getArrayValue());
		}
		
		if (approx) {
			try {
				double max = samplingDist.inverseCumulativeProbability(0.9999);
				for (int i = 1; i < 999; i++) {
					x[i] = i * max / 999;
					y[i] = samplingDist.cumulativeProbability(x[i]);
				}
				x[999] = 1e100;
				y[999] = 1.0;
			} catch (MathException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	
		needsupdate = false;
	}

	@Override
	public double logH(double t, double d) {
		if (needsupdate) {
			update();
		}
		final double logH = Math.log(constant.getArrayValue()) + samplingDist.logDensity(t - d);
//System.err.println("logh" + getID() + "(" + f.format(t) + "-" + f.format(d) + "=" + f.format(t-d) + ") = " + f4.format(logH));		
//System.err.println("logh" + getID() + "(" + f.format(t-d) + ") # = " + f4.format(logH));		
		return logH;
	}

	
	@Override
	public double simulate() throws MathException {
		double p = Randomizer.nextDouble();
		update();
		double t = samplingDist.inverseCumulativeProbability(p);
		return t;
	}
	
	public double getRate() {
		if (rate != null) {
			return rate.getArrayValue();
		} else {
			return 1.0/scale.getArrayValue();
		}
	}
	
	
	
	@Override
	protected void restore() {
		needsupdate = true;
		super.restore();
	}
	
	@Override
	protected boolean requiresRecalculation() {
		needsupdate = true;
		return super.requiresRecalculation();
	}
	
	
	public static void main(String[] args) throws MathException {
		// test difference between
		GammaHazardFunction h0 = new GammaHazardFunction();
		h0.initByName("shape","10.0","rate","6.5","C","0.75");
		GammaHazardFunction h1 = new GammaHazardFunction();
		h1.initByName("shape","10.0","rate","6.5","C","0.75","approx", true);
		double x = 0;
		while (x < 10) {
			double y1 = h1.logS(x, 0);
			double y0 = h0.logS(x , 0);
			System.err.println(y1 + " " + 0 + " " + (y0-y1));
			x += 0.01;
		}
	}
}
