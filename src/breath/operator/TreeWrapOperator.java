package breath.operator;

import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.Operator;
import beast.base.inference.OperatorSchedule;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.IntegerParameter;

@Description("Operator that wraps around a tree operator, taking in account the volume change for blocks when branch lengths change")
public class TreeWrapOperator extends Operator {
	
	final public Input<Operator> operatorInput = new Input<>("operator","tree operator to be wrapped", Validate.REQUIRED);
    final public Input<IntegerParameter> blockCountInput = new Input<>("blockcount", "number of transitions inside a block", Validate.REQUIRED);

    private IntegerParameter blockCount;
    private Operator operator;
    private TreeInterface tree;
    
	@Override
	public void initAndValidate() {
    	blockCount = blockCountInput.get();

    	operator = operatorInput.get();
    	tree = (TreeInterface) operator.getInput("tree").get();
	}

	@Override
	public double proposal() {
		double logHR = 0;
		int n = tree.getNodeCount();
		for (int i = 0; i < n - 1; i++) {
			int bc = blockCount.getValue(i);
			if (bc == 0) {
				double l = tree.getNode(i).getLength();
				logHR += Math.log(l);
			} else if (bc > 0) {
				double l = tree.getNode(i).getLength();
				logHR += 2*Math.log(l);
			}
		}

		logHR += operator.proposal();
		
		if (Double.isNaN(logHR) || logHR == Double.NEGATIVE_INFINITY) {
			return Double.NEGATIVE_INFINITY;
		}
		
		for (int i = 0; i < n - 1; i++) {
			int bc = blockCount.getValue(i);
			if (bc == 0) {
				double l = tree.getNode(i).getLength();
				logHR -= Math.log(l);
			} else if (bc > 0) {
				double l = tree.getNode(i).getLength();
				logHR -= 2*Math.log(l);
			}
		}
		
		return logHR;
	}
	
	@Override
	public double getCoercableParameterValue() {
		return operator.getCoercableParameterValue();
	}
	
	@Override
	public void setCoercableParameterValue(double value) {
		operator.setCoercableParameterValue(value);
	}
	
	@Override
	public String getPerformanceSuggestion() {
		return operator.getPerformanceSuggestion();
	}
	
	@Override
	public void accept() {
		operator.accept();
		super.accept();
	}
	
	@Override
	public void reject() {
		operator.reject();
		super.reject();
	}
	
	@Override
	public void reject(int reason) {
		operator.reject(reason);
		super.reject(reason);
	}
	
	@Override
	public void optimize(double logAlpha) {
		operator.optimize(logAlpha);
	}

	
	@Override
	public double getTargetAcceptanceProbability() {
		return operator.getTargetAcceptanceProbability();
	}
	
	@Override
	public List<StateNode> listStateNodes() {
		return operator.listStateNodes();
	}
	
    // Added for coupled MCMC
    public int get_m_nNrAccepted(){
    	return operator.get_m_nNrAccepted();
    }
    public int get_m_nNrRejected(){
    	return operator.get_m_nNrRejected();
    }
    public int get_m_nNrAcceptedForCorrection(){
    	return operator.get_m_nNrAcceptedForCorrection();
    }
    public int get_m_nNrRejectedForCorrection(){
    	return operator.get_m_nNrRejectedForCorrection();
    }
    
    // Added for coupled MCMC
    public void setAcceptedRejected(int m_nNrAccepted, int m_nNrRejected, int m_nNrAcceptedForCorrection, int m_nNrRejectedForCorrection){
    	operator.setAcceptedRejected(m_nNrAccepted, m_nNrRejected, m_nNrAcceptedForCorrection, m_nNrRejectedForCorrection);
    }
    
    
    @Override
    public void setOperatorSchedule(OperatorSchedule operatorSchedule) {
    	operator.setOperatorSchedule(operatorSchedule);
    	super.setOperatorSchedule(operatorSchedule);
    }
}
