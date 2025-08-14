package breath.evolution;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.coalescent.PopulationFunction;


/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
@Description("Coalescent intervals for a population growing linearly from time t=0; all events" +
        "are assumed to occur in negative time")
public class LinearGrowth extends PopulationFunction.Abstract {
    final public Input<Function> rateParameter = new Input<>("rate",
            "linear growth rate", Validate.REQUIRED);

    //
    // Public stuff
    //

    @Override
    public void initAndValidate() {
    	super.initAndValidate();
    	rate  = getRate();
    }

    /**
     * @return initial population size.
     */
    public double getRate() {
        rate = rateParameter.get().getArrayValue();
        return rate;
    }

    /**
     * sets rate
     *
     * @param rate new rate
     */
    public void setRate(double rate) {
        this.rate = rate;
    }


    // Implementation of abstract methods

    @Override
    public List<String> getParameterIds() {
        List<String> ids = new ArrayList<>();
        if (rateParameter.get() instanceof BEASTInterface)
            ids.add(((BEASTInterface)rateParameter.get()).getID());
        return ids;
    }

    @Override
    public double getPopSize(double t) {
        if(t>0){
            throw new IllegalArgumentException("All timings should be negative");
        }
        rate = rateParameter.get().getArrayValue();
        return(-t*rate);

    }

    @Override
    public double getIntensity(double t) {
        throw new RuntimeException("Not implemented yet");
    }

    @Override
    public double getInverseIntensity(double x) {
        throw new RuntimeException("Not implemented yet");
    }
    
    /**
     * Calculates the integral 1/N(t) dt between start and finish.
     */
    @Override
	public double getIntegral(double start, double finish) {
		return /*1/getN0() * */ Math.log(start/finish);
	}


    private double rate = 1.0;
}
