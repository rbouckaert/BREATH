package transmission2;



import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.Operator;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

@Description("Operator that adds or deletes infection")
public class AddOrDeleteInfectionOperator extends Operator {
    final public Input<TransmissionSet> transmissionsInput = new Input<>("transmissions", "represents set of transmission on the tree", Validate.REQUIRED);
    final public Input<TreeInterface> treeInput = new Input<>("tree", "tree over which to calculate a prior or likelihood");

    private TreeInterface tree;
    private TransmissionSet transmissions;
    private int eligbleInfectionCount;
	private int [] colourAtBase;

    @Override
	public void initAndValidate() {
    	tree = treeInput.get();
    	transmissions = transmissionsInput.get();
	}

    
    final static boolean debug = false;
        
	@Override
	public double proposal() {
		IntegerParameter nodeNrs = transmissions.nodeNrInput.get();
		RealParameter branchFraction = transmissions.branchFractionInput.get();
		final int n = nodeNrs.getDimension();

		length = 0;
		for (Node node : tree.getNodesAsArray()) {
			length += node.getLength();
		}

		if (Randomizer.nextBoolean()) {
			// add transmission
			
			// grow dimensions
			nodeNrs.setDimension(n + 1);
			branchFraction.setDimension(n+1);
			
			// randomly pick a branch (exclude root)
			int i = chooseBlockToInsert();
			double h = r / tree.getNode(i).getLength();
			
			nodeNrs.setValue(n, i);
			branchFraction.setValue(n, h);

			eligbleInfectionCount = calcEligableInfectionCount();
			return -Math.log(tree.getNode(i).getLength()/ length) + Math.log(1.0/eligbleInfectionCount);
			//return 0;
		} else {
			// delete transmission
			int i = chooseInfectionToRemove();
			if (i < 0) {
				// cannot find suitable candidate to remove
				return Double.NEGATIVE_INFINITY;
			}
			for (int j = i; j < n - 1; j++) {
				nodeNrs.setValue(j, nodeNrs.getValue(j+1));
				branchFraction.setValue(j, branchFraction.getValue(j+1));
			}
			nodeNrs.setDimension(n - 1);
			branchFraction.setDimension(n - 1);

			return Math.log(tree.getNode(i).getLength()/ length) - Math.log(1.0/eligbleInfectionCount);
			// return 0;
			//return Math.log(1.0/eligbleInfectionCount) - Math.log(1.0/ tree.getNodeCount());
		}
	}

	
	private double length = 0, r;
	private int chooseBlockToInsert() {
		r = Randomizer.nextDouble() * length;
		int i = 0;
		while (r > 0) {
			Node node = tree.getNode(i);
			if (r < node.getLength()) {
				return i;
			}
			r = r - node.getLength();
			i++;
		}
		throw new RuntimeException("Programmer error: should not get here");
	}

	
	private int calcEligableInfectionCount() {
		colourAtBase = new int[tree.getNodeCount()];
		transmissions.getColour(colourAtBase);
		int eligbleInfectionCount = 0;
		int n = tree.getLeafNodeCount();
		for (int i = 0; i < colourAtBase.length; i++) {
			if (transmissions.infectionCount(i) == 1) {
				if (!(colourAtBase[i] < n && colourAtBase[tree.getNode(i).getParent().getNr()] < n)) {
					eligbleInfectionCount += 1;
				}
			} else {
				eligbleInfectionCount += transmissions.infectionCount(i);
			}
		}
		return eligbleInfectionCount;
	}

	private int chooseInfectionToRemove() {
		int n = tree.getLeafNodeCount();
		
		eligbleInfectionCount = calcEligableInfectionCount();
		if (eligbleInfectionCount == 0) {
			return -1;
		}
		
		
		int k = Randomizer.nextInt(eligbleInfectionCount);
		for (int i = 0; i < colourAtBase.length; i++) {
			if (transmissions.infectionCount(i) == 1) {
				if (!(colourAtBase[i] < n && colourAtBase[tree.getNode(i).getParent().getNr()] < n)) {
					k--;
				}
			} else {
				k -= transmissions.infectionCount(i);
			}
			if (k < 0) {
				return i;
			}
		}
		throw new RuntimeException("Programmer error: should not get here");
	}

    @Override
    public List<StateNode> listStateNodes() {
        final List<StateNode> list = new ArrayList<>();
        list.add(transmissionsInput.get().branchFractionInput.get());
        list.add(transmissionsInput.get().nodeNrInput.get());
        return list;
    }
    
}
