package breath.operator;


import java.util.ArrayList;
import java.util.List;

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

@Description("Operator that randomly picks an infection and moves it elsewhere")
public class InfectionMover extends Operator {
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
		
		if (false)
		{
			// randomly pick internal node
			int nodeNr = Randomizer.nextInt(tree.getNodeCount()-1);
			if (blockCount.getArrayValue(nodeNr) < 0) {
				// immediate reject if there is no infection
				return Double.NEGATIVE_INFINITY;
			}
			// remove infection
			removeInfectionFromPath(new Node[] {tree.getNode(nodeNr)}, 0);

			// randomly pick any internal node
			int nodeNr2 = Randomizer.nextInt(tree.getNodeCount()-1);
			Node node2 = tree.getNode(nodeNr2);
			// insert infection
			insertInfectionToPath(new Node[] {node2}, node2.getLength());
			
			// make sure the colouring is valid
			colourAtBase = likelihood.getFreshColouring();		
			Validator validator = new Validator((Tree)tree, colourAtBase, blockCount, blockStartFraction, blockEndFraction);
			if (!validator.isValid(colourAtBase)) {
				// immediate reject if colouring is invalid
				return Double.NEGATIVE_INFINITY;
			}

			return 0;
		}
		

if (true)
		if (Randomizer.nextBoolean()) {
			// move infection to its sibling
			
			// randomly pick internal node
			int nodeNr = Randomizer.nextInt(tree.getNodeCount()-1);
			if (blockCount.getArrayValue(nodeNr) < 0) {
				// immediate reject if there is no infection
				return Double.NEGATIVE_INFINITY;
			}
			// remove infection
			Node node = removeInfectionFromPath(new Node[] {tree.getNode(nodeNr)}, 0);

			// put the infection on its sibling
			Node sibling = node.getParent().getLeft() == node ? node.getParent().getRight() : node.getParent().getLeft(); 

			// insert infection
			insertInfectionToPath(new Node[] {sibling}, sibling.getLength());

			// make sure the colouring is valid
			colourAtBase = likelihood.getFreshColouring();		
			Validator validator = new Validator((Tree)tree, colourAtBase, blockCount, blockStartFraction, blockEndFraction);
			if (!validator.isValid(colourAtBase)) {
				// immediate reject if colouring is invalid
				return Double.NEGATIVE_INFINITY;
			}

			return 0;
		} else {
			// move infection anywhere in the tree
			
			// randomly pick internal node
			int nodeNr = Randomizer.nextInt(tree.getNodeCount()-1);
			if (blockCount.getArrayValue(nodeNr) < 0) {
				// immediate reject if there is no infection
				return Double.NEGATIVE_INFINITY;
			}
			// remove infection
			Node node = removeInfectionFromPath(new Node[] {tree.getNode(nodeNr)}, 0);

			// randomly pick any internal node proportional to branch length
			double pathLength = 0;
			for (Node node0 : tree.getNodesAsArray()) {
				pathLength += node0.getLength();
			}

			// insert infection
			Node node2 = insertInfectionToPath(tree.getNodesAsArray(), pathLength);
			
			// make sure the colouring is valid
			colourAtBase = likelihood.getFreshColouring();		
			Validator validator = new Validator((Tree)tree, colourAtBase, blockCount, blockStartFraction, blockEndFraction);
			if (!validator.isValid(colourAtBase)) {
				// immediate reject if colouring is invalid
				return Double.NEGATIVE_INFINITY;
			}

			return Math.log(node.getLength()) - Math.log(node2.getLength());
		}
		

		
		
		
		
		double logHR = 0;
		//System.out.print(blockCount.getValue());
		int pre = blockCount.getValue(0);
		
		// get all nodes
		Node [] path = tree.getNodesAsArray();
		double pathLength = 0;
		for (Node node : path) {
			pathLength += node.getLength();
		}

		
		// 1. determine number of eligible infections 
		int eligbleInfectionCount = 0;
		for (Node node : path) {
			int bc = blockCount.getValue(node.getNr());
			if (bc == 0) {
				eligbleInfectionCount += 1;
			} else if (bc > 0) {
				eligbleInfectionCount += 2;
			}
		}

		// 2. pick one uniformly at random from eligible nodes
		int k = Randomizer.nextInt(eligbleInfectionCount);
		Node nodeWithInfectionRemoved = removeInfectionFromPath(path, k);
		
		if (nodeWithInfectionRemoved.isRoot()) {
			return Double.NEGATIVE_INFINITY;
		}
		
		logHR += Math.log(nodeWithInfectionRemoved.getLength()/pathLength) - Math.log(1.0/eligbleInfectionCount);
		
		Node insertionNode = insertInfectionToPath(path, pathLength);
		
		eligbleInfectionCount = 0;
		for (Node node : path) {
			int bc = blockCount.getValue(node.getNr());
			if (bc == 0) {
				eligbleInfectionCount += 1;
			} else if (bc > 0) {
				eligbleInfectionCount += 2;
			}
		}
		logHR += Math.log(1.0/eligbleInfectionCount) - Math.log(insertionNode.getLength()/pathLength);
		
		if (debug) {
			int post = blockCount.getValue(0);
			updateStats(pre, post);
		}
		
		// make sure the colouring is valid
		colourAtBase = likelihood.getFreshColouring();		
		Validator validator = new Validator((Tree)tree, colourAtBase, blockCount, blockStartFraction, blockEndFraction);
		if (!validator.isValid(colourAtBase)) {
			// System.err.println("x");
			return Double.NEGATIVE_INFINITY;
		}

		//System.out.print("=>"+blockCount.getValue());

		return logHR;
	}
	
	
    int [][] stats = new int[3][3];
	private void updateStats(int pre, int post) {
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

	private Node removeInfectionFromPath(Node [] path, int k) {
		//delta = -1;
		for (Node node : path) {
			int nodeNr = node.getNr();
			int bc = blockCount.getValue(nodeNr) + 1;
			k -= bc > 2 ? 2 : bc;
			if (k < 0) {
				blockCount.setValue(node.getNr(), blockCount.getValue(node.getNr()) - 1);
				if (blockCount.getValue(nodeNr) == -1) {
					return node;
				}
				if (blockCount.getValue(nodeNr) == 0) {
					if (Randomizer.nextBoolean()) {
						blockStartFraction.setValue(nodeNr, blockEndFraction.getValue(nodeNr));
					} else {
						blockEndFraction.setValue(nodeNr, blockStartFraction.getValue(nodeNr));
						
					}
					return node;
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
//				if (k == -1) {
//					delta = blockEndFraction.getValue(nodeNr);
//					double newBlockStartFraction = Randomizer.nextDouble() * delta; //blockStartFraction.getValue(nodeNr) + Randomizer.nextDouble() * (blockEndFraction.getValue(nodeNr) -blockStartFraction.getValue(nodeNr));
//					blockStartFraction.setValue(nodeNr, newBlockStartFraction);
//				} else if (k == -2) {
//					delta = 1.0 - blockStartFraction.getValue(nodeNr);
//					double newBlockEndFraction = 1.0 - Randomizer.nextDouble() * delta; //blockEndFraction.getValue(nodeNr) - Randomizer.nextDouble() * (blockEndFraction.getValue(nodeNr) - blockStartFraction.getValue(nodeNr));
//					blockEndFraction.setValue(nodeNr, newBlockEndFraction);						
//				}
				//}
				return node;
			}
			
//			if (bc == 0) {
//				k--;
//				if (k < 0) {
//					// remove from top of block
//					blockCount.setValue(nodeNr, -1);
//					return;
//				}
//			} else if (bc > 0) {
//				k -= 2;
//				if (k < -1) {
//					// remove from bottom of block
//					blockCount.setValue(node.getNr(), blockCount.getValue(node.getNr()) - 1);
//					if (blockCount.getValue(nodeNr) == 0) {
//						blockStartFraction.setValue(nodeNr, blockEndFraction.getValue(nodeNr));						
//					} else if (blockCount.getValue(nodeNr) > 0) {
//						double newBlockStartFraction = blockStartFraction.getValue(nodeNr) + Randomizer.nextDouble() * (blockEndFraction.getValue(nodeNr) -blockStartFraction.getValue(nodeNr));
//						blockStartFraction.setValue(nodeNr, newBlockStartFraction);
//					}
//					return;
//				}
//				if (k < 0) {
//					// remove from top of block
//					blockCount.setValue(node.getNr(), blockCount.getValue(node.getNr()) - 1);
//					if (blockCount.getValue(nodeNr) == 0) {
//						blockEndFraction.setValue(nodeNr, blockStartFraction.getValue(nodeNr));						
//					} else if (blockCount.getValue(nodeNr) > 0) {
//						double newBlockEndFraction = blockEndFraction.getValue(nodeNr) - Randomizer.nextDouble() * (blockEndFraction.getValue(nodeNr) - blockStartFraction.getValue(nodeNr));
//						blockEndFraction.setValue(nodeNr, newBlockEndFraction);
//					}
//					return;
//				}
//			}
		}
		throw new RuntimeException("Programmer error: should not get here");
	}	

	private Node insertInfectionToPath(Node [] path, double pathLength) {
		// insert infection uniform randomly on path 

		//delta = -1;
//		int k = Randomizer.nextInt(path.size());
//		Node node = path.get(k);
		
		double r = Randomizer.nextDouble() * pathLength;
		for (Node node : path) {
			if (node.getLength() > r) {
				int nodeNr = node.getNr();
				blockCount.setValue(nodeNr, blockCount.getValue(nodeNr) + 1);
				if (blockCount.getValue(nodeNr) == 0) {
					double f = Randomizer.nextDouble();
					blockStartFraction.setValue(nodeNr, f);// / node.getLength());
					blockEndFraction.setValue(nodeNr, f);// / node.getLength());
					return node;
				} else { // blockCount > 0
					double blockStart = Randomizer.nextDouble();
					double blockEnd = Randomizer.nextDouble();
					if (blockEnd < blockStart) {
						double tmp = blockEnd; blockEnd = blockStart; blockStart = tmp;
					}
					blockStartFraction.setValue(nodeNr, blockStart);
					blockEndFraction.setValue(nodeNr, blockEnd);					
//					if (Randomizer.nextBoolean()) {
//						delta = blockEndFraction.getValue(nodeNr);
//						double newBlockStartFraction = Randomizer.nextDouble() * delta; //blockStartFraction.getValue(nodeNr) + Randomizer.nextDouble() * (blockEndFraction.getValue(nodeNr) -blockStartFraction.getValue(nodeNr));
//						blockStartFraction.setValue(nodeNr, newBlockStartFraction);
//					} else if (k == -2) {
//						delta = 1.0 - blockStartFraction.getValue(nodeNr);
//						double newBlockEndFraction = 1.0 - Randomizer.nextDouble() * delta; //blockEndFraction.getValue(nodeNr) - Randomizer.nextDouble() * (blockEndFraction.getValue(nodeNr) - blockStartFraction.getValue(nodeNr));
//						blockEndFraction.setValue(nodeNr, newBlockEndFraction);						
//					}

//					if (blockStartFraction.getValue(nodeNr) * node.getLength()> r) {
//						blockStartFraction.setValue(nodeNr, r / node.getLength());
//						return;
//					}
//					if (blockEndFraction.getValue(nodeNr) * node.getLength() < r) {
//						blockEndFraction.setValue(nodeNr, r / node.getLength());						
//					}
					// (blockStartFraction.getValue(nodeNr) < r && r < blockEndFraction.getValue(nodeNr))
					return node;
				}
			}
			r = r - node.getLength();
		}
		throw new RuntimeException("Progammer error 2: should never get here");
	}

	


	
    protected List<Node> getPathExcludingMRCA(int nodeIndex1, int nodeIndex2) {
    	Node n1 = tree.getNode(nodeIndex1); 
    	Node n2 = tree.getNode(nodeIndex2);
    	List<Node> path = new ArrayList<>();
    	path.add(n1);
    	path.add(n2);
    	
        while (n1 != n2) {
	        double h1 = n1.getHeight();
	        double h2 = n2.getHeight();
	        if ( h1 < h2 ) {
	            n1 = n1.getParent();
	            path.add(n1);
	        } else if( h2 < h1 ) {
	            n2 = n2.getParent();
	            path.add(n2);
	        } else {
	            //zero length branches hell
	            Node n;
	            double b1 = n1.getLength();
	            double b2 = n2.getLength();
	            if( b1 > 0 ) {
	                n = n2;
	            } else { // b1 == 0
	                if( b2 > 0 ) {
	                    n = n1;
	                } else {
	                    // both 0
	                    n = n1;
	                    while( n != null && n != n2 ) {
	                        n = n.getParent();
	                    }
	                    if( n == n2 ) {
	                        // n2 is an ancestor of n1
	                        n = n1;
	                    } else {
	                        // always safe to advance n2
	                        n = n2;
	                    }
	                }
	            }
	            if( n == n1 ) {
                    n = n1 = n.getParent();
                } else {
                    n = n2 = n.getParent();
                }
	            path.add(n);
	        }
        }
        // remove MRCA itself
        while (path.remove(n1)) {}
        return path;
    }

}
