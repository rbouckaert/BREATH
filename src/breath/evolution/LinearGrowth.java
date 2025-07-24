package breath.evolution;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.parameter.RealParameter;


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



    private double rate = 1.0;
}
