package breath.operator;

import java.text.DecimalFormat;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.Operator;
import beast.base.inference.operator.kernel.KernelDistribution;
import beast.base.inference.parameter.RealParameter;

@Description("Move both origin and root height at the same time")
public class OriginAndRootHeightMover extends Operator {
    final public Input<RealParameter> originInput = new Input<>("origin", "time at which the study start above the root of tree. Assumed to be at root if not specified");
    final public Input<TreeInterface> treeInput = new Input<>("tree", "tree over which to calculate a prior or likelihood");
    final public Input<Double> scaleFactorInput = new Input<>("scaleFactor", "scaling factor: range from 0 to 1. Close to zero is very large jumps, close to 1.0 is very small jumps.", 0.75);

    final public Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);

    public final Input<KernelDistribution> kernelDistributionInput = new Input<>("kernelDistribution", "provides sample distribution for proposals", 
    		KernelDistribution.newDefaultKernelDistribution());

    protected KernelDistribution kernelDistribution;
    private double scaleFactor;

	protected double getScaler(int i, double value) {
    	return kernelDistribution.getScaler(i, value, getCoercableParameterValue());
	}
    
    
    private RealParameter origin;
    private TreeInterface tree;
    
	@Override
	public void initAndValidate() {
		tree = treeInput.get();
		origin = originInput.get();
    	kernelDistribution = kernelDistributionInput.get();
    	
    	setCoercableParameterValue(scaleFactorInput.get());
	}

	@Override
	public double proposal() {
        final Node root = tree.getRoot();                    
        final double scale = getScaler(root.getNr(), root.getHeight());
        final double newHeight = root.getHeight() * scale;

        if (newHeight < Math.max(root.getLeft().getHeight(), root.getRight().getHeight())) {
            return Double.NEGATIVE_INFINITY;
        }
        double delta = newHeight - root.getHeight();
        root.setHeight(newHeight);
        
        final double oldOrigin = origin.getValue();
        origin.setValue(oldOrigin + delta);
        
        return Math.log(scale) + Math.log(origin.getValue() / oldOrigin);
	}

    @Override
    public double getCoercableParameterValue() {
        return scaleFactor;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        scaleFactor = value; // Math.max(Math.min(value, upper), lower);
    }

    @Override
    public void optimize(double logAlpha) {
        // must be overridden by operator implementation to have an effect
    	if (optimiseInput.get()) {
	        double delta = calcDelta(logAlpha);
	        double scaleFactor = getCoercableParameterValue();
	        delta += Math.log(scaleFactor);
	        scaleFactor = Math.exp(delta);
	        setCoercableParameterValue(scaleFactor);
    	}
    }
    
    @Override
    public double getTargetAcceptanceProbability() {
    	return 0.3;
    }
    
    @Override
    public String getPerformanceSuggestion() {
        double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        double newWindowSize = getCoercableParameterValue() * ratio;

        DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10 || prob > 0.40) {
            return "Try setting scale factor to about " + formatter.format(newWindowSize);
        } else return "";
    }

}
