package breath.operator;


import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import breath.distribution.TransmissionTreeLikelihood;
import breath.distribution.Validator;

@Description("Operator that randomly pick a node and moves adjacent infections elsewhere")
public class InfectionMover3 extends Operator {
	final public Input<RealParameter> blockStartFractionInput = new Input<>("blockstart", "start of block in fraction of branch length", Validate.REQUIRED);
    final public Input<RealParameter> blockEndFractionInput = new Input<>("blockend", "end of block in fraction of branch length", Validate.REQUIRED);
    final public Input<IntegerParameter> blockCountInput = new Input<>("blockcount", "number of transitions inside a block", Validate.REQUIRED);
	final public Input<TransmissionTreeLikelihood> likelihoodInput = new Input<>("likelihood", "transmission treelikelihood containing the colouring", Validate.REQUIRED);

    private RealParameter blockStartFraction;
    private RealParameter blockEndFraction;
    private IntegerParameter blockCount;
    private TransmissionTreeLikelihood likelihood;
    private TreeInterface tree;
    private int [] colourAtBase;

    @Override
	public void initAndValidate() {
    	blockStartFraction = blockStartFractionInput.get();
    	blockEndFraction = blockEndFractionInput.get();
    	blockCount = blockCountInput.get();
    	likelihood = likelihoodInput.get();
    	tree = likelihood.treeInput.get();
	}

    
    final static boolean debug = false;
        
	@Override
	public double proposal() {
		double logHR = 0;

		// Pick an internal node.
		int i = tree.getLeafNodeCount() + Randomizer.nextInt(tree.getInternalNodeCount());
		Node node = tree.getNode(i);
		Node left = node.getLeft();
		Node right = node.getRight();

		// If there are no infections on any of its adjacent branches then the move fails.
		// Don’t count the root branch at any time, so the root node has only two adjacent branches while all others have three.
		boolean nodeHasInfections = (node.isRoot() ? false : blockCount.getValue(node.getNr()) + 1 > 0);
		boolean leftHasInfections = blockCount.getValue(left.getNr()) + 1>0;
		boolean rightHasInfections = blockCount.getValue(right.getNr()) + 1>0;
		int numBranchesWithInfections =  (nodeHasInfections ? 1 : 0) + 
				(leftHasInfections ? 1 : 0) + 
				(rightHasInfections ? 1 : 0);
		if (numBranchesWithInfections == 0) {
			return Double.NEGATIVE_INFINITY;
		}
				
		// Pick an infection from an adjacent branch
		// and move it to another one
		// Redraw the block boundaries if you need to.
		int infection = Randomizer.nextInt(numBranchesWithInfections);
		if (nodeHasInfections && infection == 0) {
			removeInfection(node);
			if (Randomizer.nextBoolean()) {
				insertInfection(left);
			} else {
				insertInfection(right);
			}
		} else {
			infection -= nodeHasInfections ? 1 : 0;
			if (leftHasInfections && infection == 0) {
				removeInfection(left);
				if (Randomizer.nextBoolean() && !node.isRoot()) {
					insertInfection(node);
				} else {
					insertInfection(right);
				}
			} else {
				removeInfection(right);
				if (Randomizer.nextBoolean() && !node.isRoot()) {
					insertInfection(node);
				} else {
					insertInfection(left);
				}
			}
		}

		// make sure the colouring is valid
		colourAtBase = likelihood.getFreshColouring();		
		Validator validator = new Validator((Tree)tree, colourAtBase, blockCount, blockStartFraction, blockEndFraction);
		if (!validator.isValid(colourAtBase)) {
			// System.err.println("x");
			return Double.NEGATIVE_INFINITY;
		}
		
		// Hastings ratios are simple, and origin and destination branches are uniquely specified by the choice of node. 		
		boolean nodeHasInfections2 = (node.isRoot() ? false : blockCount.getValue(node.getNr()) + 1 > 0);
		boolean leftHasInfections2 = blockCount.getValue(left.getNr()) + 1>0;
		boolean rightHasInfections2 = blockCount.getValue(right.getNr()) + 1>0;
		int numBranchesWithInfections2 =  (nodeHasInfections2 ? 1 : 0) + 
				(leftHasInfections2 ? 1 : 0) + 
				(rightHasInfections2 ? 1 : 0);
		logHR = Math.log(numBranchesWithInfections) - Math.log(numBranchesWithInfections2);

		return logHR;
	}
	
	
	// reduce infections by 1 on branch above given node
	private void removeInfection(Node node) {
		int nodeNr = node.getNr();
		blockCount.setValue(nodeNr, blockCount.getValue(nodeNr) - 1);
		if (blockCount.getValue(nodeNr) == -1) {
			return;
		}
		if (blockCount.getValue(nodeNr) == 0) {
			if (Randomizer.nextBoolean()) {
				blockStartFraction.setValue(nodeNr, blockEndFraction.getValue(nodeNr));
			} else {
				blockEndFraction.setValue(nodeNr, blockStartFraction.getValue(nodeNr));
				
			}
			return;
		}
		
		// shrink the block if it is on a boundary
		// which is just as well as positions 1 and 2
		//if (k == -1 || k == -2) { // blockCount.getValue(node.getNr()) - k == 1) {
		//	if (Randomizer.nextBoolean()) {
		double blockStart = Randomizer.nextDouble();
		double blockEnd = Randomizer.nextDouble();
		if (blockEnd < blockStart) {
			double tmp = blockEnd; blockEnd = blockStart; blockStart = tmp;
		}
		blockStartFraction.setValue(nodeNr, blockStart);
		blockEndFraction.setValue(nodeNr, blockEnd);					
	}
	
	private Node insertInfection(Node node) {
		int nodeNr = node.getNr();
		blockCount.setValue(nodeNr, blockCount.getValue(nodeNr) + 1);
		if (blockCount.getValue(nodeNr) == 0) {
			double f = Randomizer.nextDouble();
			blockStartFraction.setValue(nodeNr, f);
			blockEndFraction.setValue(nodeNr, f);
			return node;
		} else { // blockCount > 0
			double blockStart = Randomizer.nextDouble();
			double blockEnd = Randomizer.nextDouble();
			if (blockEnd < blockStart) {
				double tmp = blockEnd; blockEnd = blockStart; blockStart = tmp;
			}
			blockStartFraction.setValue(nodeNr, blockStart);
			blockEndFraction.setValue(nodeNr, blockEnd);					
			return node;
		}
	}


}
