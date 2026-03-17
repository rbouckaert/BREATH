package breath.operator;


import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
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
    final public Input<Boolean> useBranchLengthInput = new Input<>("useBranchLength", "use Branch Length in HR", true);

    private RealParameter blockStartFraction;
    private RealParameter blockEndFraction;
    private IntegerParameter blockCount;
    private TransmissionTreeLikelihood likelihood;
    private TreeInterface tree;
    private int[] colourAtBase;
    private boolean useBranchLength;

    @Override
    public void initAndValidate() {
        blockStartFraction = blockStartFractionInput.get();
        blockEndFraction = blockEndFractionInput.get();
        blockCount = blockCountInput.get();
        likelihood = likelihoodInput.get();
        tree = likelihood.treeInput.get();
        useBranchLength = useBranchLengthInput.get();
    }


    final static boolean debug = false;

    @Override
    public double proposal() {

        if (false) {
            // randomly pick internal node
            int nodeNr = Randomizer.nextInt(tree.getNodeCount() - 1);
            if (blockCount.getArrayValue(nodeNr) < 0) {
                // immediate reject if there is no infection
                return Double.NEGATIVE_INFINITY;
            }
            // remove infection
            removeInfectionFromPath(new Node[]{tree.getNode(nodeNr)}, 0);

            // randomly pick any internal node
            int nodeNr2 = Randomizer.nextInt(tree.getNodeCount() - 1);
            Node node2 = tree.getNode(nodeNr2);
            // insert infection
            insertInfectionToPath(new Node[]{node2}, node2.getLength());

            // make sure the colouring is valid
            colourAtBase = likelihood.getFreshColouring();
            Validator validator = new Validator((Tree) tree, colourAtBase, blockCount, blockStartFraction, blockEndFraction);
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
                int eligibleNodesBefore = 0;
                for(int i = 0; i < tree.getNodeCount()-1; i++) {
                    if(blockCount.getArrayValue(i) >= 0){
                        eligibleNodesBefore++;
                    }
                }

                int nodeNr = Randomizer.nextInt(tree.getNodeCount()-1);
                if (blockCount.getArrayValue(nodeNr) < 0) {
                    // immediate reject if there is no infection
                    return Double.NEGATIVE_INFINITY;
                }
                // remove infection
                Node source = removeInfectionFromPath(new Node[]{tree.getNode(nodeNr)}, 0);

                // put the infection on its sibling
                Node sibling = source.getParent().getLeft() == source ? source.getParent().getRight() : source.getParent().getLeft();

                // insert infection
                Node target = insertInfectionToPath(new Node[]{sibling}, sibling.getLength());


                // make sure the colouring is valid
                colourAtBase = likelihood.getFreshColouring();
                Validator validator = new Validator((Tree) tree, colourAtBase, blockCount, blockStartFraction, blockEndFraction);
                if (!validator.isValid(colourAtBase)) {
                    // immediate reject if colouring is invalid
                    return Double.NEGATIVE_INFINITY;
                }

                int eligibleNodesAfter = 0;
                for(int i = 0; i < tree.getNodeCount()-1; i++) {
                    if(blockCount.getArrayValue(i) >= 0){
                        eligibleNodesAfter++;
                    }
                }

                double logHR = Math.log(eligibleNodesBefore) - Math.log(eligibleNodesAfter);
                if (blockCount.getArrayValue(source.getNr()) < 0) {
                    double l = useBranchLength ? source.getLength() : 1;
                    logHR += Math.log(1 / l);
                } else if (blockCount.getArrayValue(source.getNr()) == 0) {
                    double l = useBranchLength ? source.getLength() : 1;
                    logHR += Math.log(2 / l);
                }
                if (blockCount.getArrayValue(target.getNr()) == 0) {
                    double l = useBranchLength ? target.getLength() : 1;
                    logHR += -Math.log(1 / l);
                } else if (blockCount.getArrayValue(target.getNr()) == 1) {
                    double l = useBranchLength ? target.getLength() : 1;
                    logHR += -Math.log(2 / l);
                }


                double HR = Math.exp(logHR);
                Log.debug(
                            String.format("moveToSibling: node=%d, sourceLength=%.5f, targetLength=%.5f, logHR=%.5f, HR=%.5f",
                                    source.getNr(), source.getLength(), target.getLength(), logHR, HR)
                    );

                return logHR;
            } else {
                // move infection anywhere in the tree

                // randomly pick internal node
                int eligibleNodesBefore = 0;
                for(int i = 0; i < tree.getNodeCount()-1; i++) {
                    if(blockCount.getArrayValue(i) >= 0){
                        eligibleNodesBefore++;
                    }
                }

                int nodeNr = Randomizer.nextInt(tree.getNodeCount() - 1);
                if (blockCount.getArrayValue(nodeNr) < 0) {
                    // immediate reject if there is no infection
                    return Double.NEGATIVE_INFINITY;
                }
                // remove infection
                Node source = removeInfectionFromPath(new Node[]{tree.getNode(nodeNr)}, 0);

                // randomly pick any internal node proportional to branch length
                double pathLength = 0;
                for (Node node0 : tree.getNodesAsArray()) {
                    pathLength += node0.getLength();
                }

                // insert infection
                Node target = insertInfectionToPath(tree.getNodesAsArray(), pathLength);

                // make sure the colouring is valid
                colourAtBase = likelihood.getFreshColouring();
                Validator validator = new Validator((Tree) tree, colourAtBase, blockCount, blockStartFraction, blockEndFraction);
                if (!validator.isValid(colourAtBase)) {
                    // immediate reject if colouring is invalid
                    return Double.NEGATIVE_INFINITY;
                }

                int eligibleNodesAfter = 0;
                for(int i = 0; i < tree.getNodeCount()-1; i++) {
                    if(blockCount.getArrayValue(i) >= 0){
                        eligibleNodesAfter++;
                    }
                }

                double logHR = Math.log(eligibleNodesBefore) - Math.log(eligibleNodesAfter);
                if (blockCount.getArrayValue(source.getNr()) < 0) {
                    double l = useBranchLength ? source.getLength() : 1;
                    logHR += Math.log(1 / l);
                } else if (blockCount.getArrayValue(source.getNr()) == 0) {
                    double l = useBranchLength ? source.getLength() : 1;
                    logHR += Math.log(2 / l);
                }
                if (blockCount.getArrayValue(target.getNr()) == 0) {
                    double l = useBranchLength ? target.getLength() : 1;
                    logHR -= Math.log(1 / l);
                } else if (blockCount.getArrayValue(target.getNr()) == 1) {
                    double l = useBranchLength ? target.getLength() : 1;
                    logHR -= Math.log(2 / l);
                }
                logHR += Math.log(source.getLength()) - Math.log(target.getLength());

                double HR = Math.exp(logHR);
                Log.debug(String.format("moveAnywhere: source=%d, target=%d, sourceLength=%.5f, targetLength=%.5f, logHR=%.5f, HR=%.5f",
                            source.getNr(), target.getNr(), source.getLength(), target.getLength(), logHR, HR));

                return logHR;
            }

        double logHR = 0;
        //System.out.print(blockCount.getValue());
        int pre = blockCount.getValue(0);

        // get all nodes
        Node[] path = tree.getNodesAsArray();
        double pathLength = 0;
        for (Node node : path) {
            pathLength += node.getLength();
        }


        // 1. determine number of eligible infections
        int eligibleInfectionCount = 0;
        for (Node node : path) {
            int bc = blockCount.getValue(node.getNr());
            if (bc == 0) {
                eligibleInfectionCount += 1;
            } else if (bc > 0) {
                eligibleInfectionCount += 2;
            }
        }

        // 2. pick one uniformly at random from eligible nodes
        int k = Randomizer.nextInt(eligibleInfectionCount);
        Node nodeWithInfectionRemoved = removeInfectionFromPath(path, k);

        if (nodeWithInfectionRemoved.isRoot()) {
            return Double.NEGATIVE_INFINITY;
        }

        logHR += Math.log(nodeWithInfectionRemoved.getLength() / pathLength) - Math.log(1.0 / eligibleInfectionCount);

        Node insertionNode = insertInfectionToPath(path, pathLength);

        eligibleInfectionCount = 0;
        for (Node node : path) {
            int bc = blockCount.getValue(node.getNr());
            if (bc == 0) {
                eligibleInfectionCount += 1;
            } else if (bc > 0) {
                eligibleInfectionCount += 2;
            }
        }
        logHR += Math.log(1.0 / eligibleInfectionCount) - Math.log(insertionNode.getLength() / pathLength);

        if (debug) {
            int post = blockCount.getValue(0);
            updateStats(pre, post);
        }

        // make sure the colouring is valid
        colourAtBase = likelihood.getFreshColouring();
        Validator validator = new Validator((Tree) tree, colourAtBase, blockCount, blockStartFraction, blockEndFraction);
        if (!validator.isValid(colourAtBase)) {
            // System.err.println("x");
            return Double.NEGATIVE_INFINITY;
        }

        //System.out.print("=>"+blockCount.getValue());

        double HR = Math.exp(logHR);
        Log.debug(String.format("pickRandomInfection: removedNode=%d, insertedNode=%d, logHR=%.5f, HR=%.5f",
                    nodeWithInfectionRemoved.getNr(), insertionNode.getNr(), logHR, HR));

        return logHR;
    }


    int[][] stats = new int[3][3];

    private void updateStats(int pre, int post) {
        stats[1 + pre][1 + post]++;
        int total = stats[0][0] + stats[0][1] + stats[1][0] + stats[1][1] + stats[1][2] + stats[2][1] + stats[2][2];
        if (stats[0][0] > 0 && (total) % 1000000 == 0) {
            StringBuilder b = new StringBuilder();
            double t = stats[0][0] + stats[0][1];
            b.append(stats[0][0] / t + " " + stats[0][1] / t + ";");
            t = stats[1][0] + stats[1][1] + stats[1][2];
            b.append(stats[1][0] / t + " " + stats[1][1] / t + " " + stats[1][2] / t + ";");
            t = stats[2][1] + stats[2][2];
            b.append(stats[2][1] / t + " " + stats[2][2] / t + ";");

            t = total;
            b.append((stats[0][0] + stats[0][1]) / t + " ");
            b.append((stats[1][0] + stats[1][1] + stats[1][2]) / t + " ");
            b.append((stats[2][1] + stats[2][2]) / t);
            b.append("\n");

            System.out.println(b);
        }
    }

    private Node removeInfectionFromPath(Node[] path, int k) {
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
                    double f = Randomizer.nextDouble();
                    blockStartFraction.setValue(nodeNr, f);// / node.getLength());
                    blockEndFraction.setValue(nodeNr, f);// / node.getLength());
                    return node;
                }

                // shrink the block if it is on a boundary
                // which is just as well as positions 1 and 2
                //if (k == -1 || k == -2) { // blockCount.getValue(node.getNr()) - k == 1) {
                //	if (Randomizer.nextBoolean()) {
                double blockStart = Randomizer.nextDouble();
                double blockEnd = Randomizer.nextDouble();
                if (blockEnd < blockStart) {
                    double tmp = blockEnd;
                    blockEnd = blockStart;
                    blockStart = tmp;
                }
                blockStartFraction.setValue(nodeNr, blockStart);
                blockEndFraction.setValue(nodeNr, blockEnd);
                return node;
            }
        }
        throw new RuntimeException("Programmer error: should not get here");
    }

    private Node insertInfectionToPath(Node[] path, double pathLength) {
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
                        double tmp = blockEnd;
                        blockEnd = blockStart;
                        blockStart = tmp;
                    }
                    blockStartFraction.setValue(nodeNr, blockStart);
                    blockEndFraction.setValue(nodeNr, blockEnd);
                    return node;
                }
            }
            r = r - node.getLength();
        }
        throw new RuntimeException("Programmer error 2: should never get here");
    }
}
