package breath.operator;

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
import breath.distribution.ColourProvider;

@Description("Operator that moves block parameters of a transmission tree")
public class BlockOperator extends Operator {
	final public Input<RealParameter> blockStartFractionInput = new Input<>("blockstart", "start of block in fraction of branch length", Validate.REQUIRED);
    final public Input<RealParameter> blockEndFractionInput = new Input<>("blockend", "end of block in fraction of branch length", Validate.REQUIRED);
    final public Input<IntegerParameter> blockCountInput = new Input<>("blockcount", "number of transitions inside a block", Validate.REQUIRED);
    final public Input<Boolean> keepConstantCountInput = new Input<>("keepconstantcount", "if true, for every deleting there is an insertion to keep total sum of block counts constant", false);
    final public Input<TreeInterface> treeInput = new Input<>("tree", "tree over which to calculate a prior or likelihood", Validate.REQUIRED);

    private RealParameter blockStartFraction;
    private RealParameter blockEndFraction;
    private IntegerParameter blockCount;
    private TreeInterface tree;
    private double lowerStart, upperStart;
    private double lowerEnd, upperEnd;
    
    @Override
	public void initAndValidate() {
    	blockStartFraction = blockStartFractionInput.get();
    	blockEndFraction = blockEndFractionInput.get();
    	blockCount = blockCountInput.get();
    	tree = treeInput.get();
    	
    	lowerStart = blockStartFraction.getLower();
    	if (lowerStart < 0) {
    		lowerStart = 0;
    	}
    	if (lowerStart > 1) {
    		throw new IllegalArgumentException("lower bound of block start should be less than 1");
    	}
    	upperStart = blockStartFraction.getUpper();
    	if (upperStart < 0) {
    		upperStart = 0;
    	}
    	if (upperStart > 1) {
    		throw new IllegalArgumentException("upper bound of block start should be less than 1");
    	}
    	if (upperStart < lowerStart) {
    		throw new IllegalArgumentException("upper bound of block start should be higher than lower bound");
    	}

    	lowerEnd = blockEndFraction.getLower();
    	if (lowerEnd < 0) {
    		lowerEnd = 0;
    	}
    	if (lowerEnd > 1) {
    		throw new IllegalArgumentException("lower bound of block end should be less than 1");
    	}
    	upperEnd = blockEndFraction.getUpper();
    	if (upperEnd < 0) {
    		upperEnd = 0;
    	}
    	if (upperEnd > 1) {
    		throw new IllegalArgumentException("upper bound of block end should be less than 1");
    	}
    	if (upperEnd < lowerEnd) {
    		throw new IllegalArgumentException("upper bound of block end should be higher than lower bound");
    	}
    	
    	if (lowerStart > lowerEnd) {
    		throw new IllegalArgumentException("lower bound of block start should be lower than lower bound of block end");
    	}
    	if (upperStart > upperEnd) {
    		throw new IllegalArgumentException("upper bound of block start should be lower than upper bound of block end");
    	}
    	
    }

	@Override
	public double proposal() {
//		if (true)
//			if (Randomizer.nextBoolean()) {
//				int [] i = chooseInfectionToRemove();
//				if (i == null) {
//					return Double.NEGATIVE_INFINITY;
//				}
//				return removeInfection0(i);
//			} else {
//				int k = chooseBlockToInsert();
//				return insertInfection0(k);
//			}

		
		
		if (Randomizer.nextBoolean()) {
			// move block boundaries for a branch that has a non-empty block
			
			int i = Randomizer.nextInt(blockStartFraction.getDimension());

			int attempts = 0;
			while (blockCount.getValue(i) == -1 && attempts < 100) {
				i = Randomizer.nextInt(blockStartFraction.getDimension());
				attempts++;
			}
				
			// only move start and end fraction but not block count
			switch (blockCount.getValue(i)) {
			case -1:
				// nothing to do since start and end fractions are ignored
				break;
			case 0:
				// make sure start == end fraction after proposal
				double f = lowerStart + Randomizer.nextDouble() * (upperEnd - lowerStart);
				blockStartFraction.setValue(i, f);
				blockEndFraction.setValue(i, f);
				break;
			default:
				double blockStart = Randomizer.nextDouble();
				double blockEnd = Randomizer.nextDouble();
				if (blockEnd < blockStart) {
					double tmp = blockEnd; blockEnd = blockStart; blockStart = tmp;
				}
				blockStartFraction.setValue(i, blockStart);
				blockEndFraction.setValue(i, blockEnd);					
			}
			return 0;
		}
		
		if (keepConstantCountInput.get()) {
			// NB By deafult keepConstantCount == false, so should not get here
			// but if it does 
			// remove one infection safely (so that colouring is still valid) 
			// then add one infection
			int pre = blockCount.getValue(0);

			int [] i = chooseInfectionToRemove();
			double logHR = 0;
			if (i != null) {
				logHR += removeInfection(i);
			} else {
				return Double.NEGATIVE_INFINITY;
			}
			int k = chooseBlockToInsert();
			logHR += insertInfection(k);
			
			int post = blockCount.getValue(0);
			updateStats(pre, post);

			
			return 0*logHR;
		} else	if (Randomizer.nextBoolean()) {
			// remove one infection
			
			// first, find an infection to remove
			// possibly, no infection can be removed safely 
			// (i.e. such that the remaining colouring is valid)
			int [] i = chooseInfectionToRemove();
			if (i == null) {
				return Double.NEGATIVE_INFINITY;
			}
			
			// found a good candidate, so remove it
			return removeInfection(i);
		} else {
			// add one infection
			
			// it is always possible to add infections, so no special 
			// case here (unlike when removing infections)
			int k = chooseBlockToInsert();
			return insertInfection(k);
		}
		
	}

    int [][] stats = new int[3][3];
	private void updateStats(int pre, int post) {
		if (true) return;
		stats[1+pre][1+post]++;
		int total = stats[0][0] + stats[0][1] + stats[1][0] + stats[1][1] + stats[1][2] + stats[2][1] + stats[2][2];
		if (stats[0][0] > 0 && (total) % 1000000 == 0) {
			StringBuilder b = new StringBuilder();
			double t = stats[0][0] + stats[0][1];
			b.append(stats[0][0]/t + " " + stats[0][1]/t + ";");
			t = stats[1][0] + stats[1][1] + stats[1][2];
			b.append(stats[1][0]/t + " " + stats[1][1]/t + " " + stats[1][2]/t + ";");
			t = stats[2][1] + stats[2][2];
			b.append(stats[2][1]/t + " " + stats[2][2]/t + ";");

			t = total;
			b.append((stats[0][0] + stats[0][1])/t + " " );
			b.append((stats[1][0] + stats[1][1] + stats[1][2])/t+ " " );
			b.append((stats[2][1] + stats[2][2])/t);
			b.append("\n");
			
			System.out.println(b.toString());
		}
	}

	
	private int eligibleInfectionCount = 0;
	
	
	private int [] calcEligibleInfectionCount() {
		// colour the tree based on current infections
		int [] colourAtBase = new int[tree.getNodeCount()];
		int n = tree.getLeafNodeCount();
		ColourProvider.getColour(tree.getRoot(), blockCount, n, colourAtBase);
		
		// go through the whole tree, and for each branch
		// check if removing the infection results in a valid colouring
		// If so, add to eligbleInfectionCount
		eligibleInfectionCount = 0;
		for (int i = 0; i < blockCount.getDimension(); i++) {
			if (blockCount.getValue(i) == 0) {
				// one infection on this branch, that if removed, can lead to invalid colouring
				// 1. colour at base = a sampled host colour (if < n), and
				// 2. colour at parent = another sampled host colour (if < n)
				if (!(colourAtBase[i] < n && !tree.getNode(i).isRoot() && colourAtBase[tree.getNode(i).getParent().getNr()] < n)) {
					// otherwise, it can be removed, and the infection can be added to eligbleInfectionCount
					eligibleInfectionCount += 1;
				}
			} else if (blockCount.getValue(i) > 0) {
				// more than one infection on this branch, so we can safely remove 
				// start or end infection of the block: i.e. two possibilities
				eligibleInfectionCount += 2;
			//} else {
				// cannot remove infection and leave a valid colouring
				// so leave eligbleInfectionCount unchanged
			}
		}
		return colourAtBase;
	}
	
	private int[] chooseInfectionToRemove() {
		// choose infection to be removed such that the remaining infections still leave a valid infection history 
		// (i.e. there is no path between any pair of leaves that does not contain an infection)
		int [] colourAtBase = calcEligibleInfectionCount();
		int n = tree.getLeafNodeCount();
		if (eligibleInfectionCount == 0) {
			return null;
		}
		
		// randomly pick one of the eligible infections to remove
		int k = Randomizer.nextInt(eligibleInfectionCount);
		
		// loop through the eligible infections, till we find the k-th one
		for (int i = 0; i < blockCount.getDimension(); i++) {
			if (blockCount.getValue(i) == 0) {
				if (!(colourAtBase[i] < n && !tree.getNode(i).isRoot() && colourAtBase[tree.getNode(i).getParent().getNr()] < n)) {
					k--;
				}
			} else if (blockCount.getValue(i) > 0) {
				k -= 2;
			}
			if (k < 0) {
				// found the branch containing the eligble infection
				return new int[] {i, k};
			}
		}
		throw new RuntimeException("Programmer error: should not get here");
	}
	
	private int chooseBlockToInsert() {
		// choose random location on branch proportional to lengths of branches
		
		// first calculate length of tree
		double length = 0;
		for (Node node : tree.getNodesAsArray()) {
			length += node.getLength();
		}
		
		// random point on length
		double r = Randomizer.nextDouble() * length;
		
		// find the node associated with r
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

	
	/** insert infection on branch i **/
	private double insertInfection(int i) {
		
		switch (blockCount.getValue(i)) {
		case -1:
			// add infection on branch without any infection
			// a random location on the branch must be chosen to put it
			// start and end of block must be equal
			blockCount.setValue(i, 0);
			double f = Randomizer.nextDouble();
			blockStartFraction.setValue(i, f);
			blockEndFraction.setValue(i, f);
			break;
			
		case 0:
			// add infection to branch already containing an infection
			// since block start == block end, we need to choose new values for
			// start and end
			blockCount.setValue(i, 1);
			
			double blockStart = Randomizer.nextDouble();
			double blockEnd = Randomizer.nextDouble();
			if (blockEnd < blockStart) {
				double tmp = blockEnd; blockEnd = blockStart; blockStart = tmp;
			}
			blockStartFraction.setValue(i, blockStart);
			blockEndFraction.setValue(i, blockEnd);					
			break;
		
		default:
			// add infection to block already containing 2 infections
			// assume it goes inside the block, so no need to update block boundaries
			blockCount.setValue(i, blockCount.getValue(i)+1);
			
		}

		// calculate the number of infections that can be removed safely
		// after we added this infection
		calcEligibleInfectionCount();
		double length = 0;
		for (Node node : tree.getNodesAsArray()) {
			length += node.getLength();
		}
		//              probability this infection got selected for removal
		// HR = ----------------------------------------------------------------------------
		//      probability density the infection gets inserted at this branch at this point
		// return log(HR)		
		return Math.log(1.0/ eligibleInfectionCount)
			   - Math.log(tree.getNode(i).getLength() / length);
	} // insertInfection

	private double removeInfection(int [] infection) {
		int i = infection[0];
		switch (blockCount.getValue(i)) {
		case -1:
			// do nothing, should not get here
			return 0; 
			
		case 0:
			// remove infection: no infections left
			blockCount.setValue(i, -1);
			break;
			
		case 1:
			// remove infection: 1 infection left, so block start and end time becomes the same
			blockCount.setValue(i, 0);
			if (Randomizer.nextBoolean()) {
				blockStartFraction.setValue(i, blockEndFraction.getValue(i));
			} else {
				blockEndFraction.setValue(i, blockStartFraction.getValue(i));
			}
			// potentially, this should draw a new random value
			// for block start == block end to be symmetric with
			// adding an infection, like so:
//			double blockStart = Randomizer.nextDouble();
//			blockStartFraction.setValue(i, blockStart);
//			blockEndFraction.setValue(i, blockStart);					

			break;
		default:
			// remove infection and 2 or more infections left
			// assume we remove one inside the block, so no change to boundaries
			blockCount.setValue(i, blockCount.getValue(i)-1);
		}
		
		
		// calculate length of tree
		double length = 0;
		for (Node node : tree.getNodesAsArray()) {
			length += node.getLength();
		}
		
		//      probability density the infection gets inserted at this branch at this point
		// HR = ----------------------------------------------------------------------------
		//              probability this infection got selected for removal
		// return log(HR)		
		return Math.log(tree.getNode(i).getLength() / length) 
				- Math.log(1.0/ eligibleInfectionCount);
	} // removeInfection


	private double insertInfection0(int i) {
			// add infection
			blockCount.setValue(i, blockCount.getValue(i)+1);
			double k = blockCount.getValue(i) + 2;
			blockStartFraction.setValue(i, 1.0/k);
			blockEndFraction.setValue(i, (k-1.0)/k);
			return 0;
	} // insertInfection

	private double removeInfection0(int [] infection) {
		int i = infection[0];

		int bc = blockCount.getValue(i);
		if (bc == -1) {
			throw new RuntimeException("Programmer error");
		}
		// remove infection
		blockCount.setValue(i, blockCount.getValue(i)-1);

		double k = blockCount.getValue(i) + 2;
		blockStartFraction.setValue(i, 1.0/k);
		blockEndFraction.setValue(i, (k-1.0)/k);
		return 0;
	} // removeInfection

	@Override
	public List<StateNode> listStateNodes() {
        final List<StateNode> list = new ArrayList<>();
        list.add(blockCount);
        list.add(blockStartFraction);
        list.add(blockEndFraction);
		return list;
	}
}
