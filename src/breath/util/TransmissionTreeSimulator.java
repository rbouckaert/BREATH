package breath.util;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.GammaDistribution;
import org.apache.commons.math.distribution.GammaDistributionImpl;
import org.apache.commons.math.distribution.PoissonDistribution;
import org.apache.commons.math.distribution.PoissonDistributionImpl;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Function.Constant;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeParser;
import beast.base.evolution.tree.coalescent.ConstantPopulation;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.evolution.tree.coalescent.PopulationFunction.Abstract;
import beast.base.inference.Runnable;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Binomial;
import beast.base.util.HeapSort;
import beast.base.util.Randomizer;
import beastfx.app.tools.Application;
import beastfx.app.util.OutFile;
import breath.distribution.ColourProvider;
import breath.distribution.GammaHazardFunction;
import breath.distribution.TransmissionTreeLikelihood;
import breath.distribution.TransmissionTreeLikelihood1;
import breath.evolution.LinearGrowth;

@Description("Simulates transmission tree with colouring and block counts")
public class TransmissionTreeSimulator extends Runnable {
	final public Input<Function> endTimeInput = new Input<>("endTime", "end time of the study", new Constant("1.0"));
	final public Input<Function> popSizeInput = new Input<>("popSize",
			"population size governing the coalescent process", new Constant("0.1"));
	final public Input<Function> linearGrowthRateInput = new Input<>("linearGrowthRate",
			"if specified, the popSize is ignored and a linear growth rate is assumed "
			+ "to govern the coalescent process");

	final public Input<Function> sampleShapeInput = new Input<>("sampleShape",
			"shape parameter of the sampling intensity function", new Constant("2.0"));
	final public Input<Function> sampleRateInput = new Input<>("sampleRate",
			"rate parameter of the sampling intensity function", new Constant("5.0"));
	final public Input<Function> sampleConstantInput = new Input<>("sampleConstant",
			"constant multiplier of the sampling intensity function", new Constant("0.75"));

	final public Input<Function> transmissionShapeInput = new Input<>("transmissionShape",
			"shape parameter of the transmission intensity function", new Constant("2.5"));
	final public Input<Function> transmissionRateInput = new Input<>("transmissionRate",
			"rate parameter of the transmission intensity function", new Constant("10.0"));
	final public Input<Function> transmissionConstantInput = new Input<>("transmissionConstant",
			"constant multiplier of the transmission intensity function", new Constant("1.5"));

	final public Input<OutFile> outputInput = new Input<>("out","output file. Print to stdout if not specified");
	final public Input<OutFile> traceOutputInput = new Input<>("trace", "trace output file with end time, tree heights and tree lengths, or stdout if not specified", new OutFile("[[none]]"));
	final public Input<Long> seedInput = new Input<>("seed","random number seed used to initialise the random number generator");
	final public Input<Integer> maxAttemptsInput = new Input<>("maxAttempts",
			"maximum number of attempts to generate coalescent sub-trees", 1000);
	final public Input<Integer> taxonCountInput = new Input<>("taxonCount", "generate tree with taxonCount number of taxa. Ignored if negative", -1);
	final public Input<Integer> maxTaxonCountInput = new Input<>("maxTaxonCount", "reject any tree with more than this number of taxa. Ignored if negative", -1);
	final public Input<Integer> treeCountInput = new Input<>("treeCount", "generate treeCount number of trees", 1);
	final public Input<Boolean> directOnlyInput = new Input<>("directOnly", "consider direct infections only, if false block counts are ignored", true);
	final public Input<Boolean> quietInput = new Input<>("quiet", "suppress some screen output", false);
	final public Input<Boolean> calcLogPInput = new Input<>("calcLogP", "calculate transmission likelihood for generated trees", false);

	private Node root;
	private Map<Node, Integer> colourMap;
	private int nodeCount;
	private double logP;
	
	private final boolean debug = false;
	private String newick0;
	
	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		int taxonCount = taxonCountInput.get();
		int maxTaxonCount = maxTaxonCountInput.get();
		if (maxTaxonCount <= 0) {
			maxTaxonCount = Integer.MAX_VALUE;
		}
		if (seedInput.get() != null) {
			Randomizer.setSeed(seedInput.get());
		}
		PrintStream out = System.out;
		if (outputInput.get() != null) {
			out = new PrintStream(outputInput.get());
		}
		PrintStream traceout = System.out;
		if (traceOutputInput.get() != null && !traceOutputInput.get().getName().equals("[[none]]")) {
			Log.warning("Writing to file " + traceOutputInput.get().getPath());
			traceout = new PrintStream(traceOutputInput.get());
		}
    	traceout.print("Sample\t");
    	boolean printTaxaInTrace = false;
    	if (printTaxaInTrace && taxonCount > 0) {
    		for (int i = 0; i < taxonCount; i++) {
    			traceout.print("t" + format(i+1) + "\t");
    		}
    	}
    	traceout.println("endTime\tTree.height\tTree.treeLength\torigin\tlogP" + (calcLogPInput.get()? "\tlogP2" : ""));
		
    	
    	Map<Integer, Integer> taxonCounts = new HashMap<>();
    	Map<Integer, Integer> infectionCounts = new HashMap<>();
		
		PrintStream out0 = debug? new PrintStream("/tmp/out0.nwk") : null;

    	
    	long attempts = 0;
		for (int i = 0; i < treeCountInput.get(); i++) {
			int k;
			double logP;
			do {
				this.logP = 0;
				logP = runOnce(maxTaxonCount);
				k = root.getAllLeafNodes().size();
				if (!taxonCounts.containsKey(k)) {
					taxonCounts.put(k, 0);
				}
				taxonCounts.put(k, taxonCounts.get(k) + 1);
				attempts++;
				if (attempts % 100000 == 0) {
					// reportAttempts(taxonCounts, );
				}
			} while (Double.isInfinite(logP) || taxonCount > 0 && taxonCount != k);

			System.out.println(i+"\t" + nodeCount + "\t" + logP);
			
			// convert to binary tree
			while (root.getChildCount() == 1) {
				root = root.getChild(0);
			}
			root.setParent(null);
			String newick = toNewick(root);
			
			int infectionCount = infectionCount(root);
			if (!infectionCounts.containsKey(k)) {
				infectionCounts.put(k, 0);
			}
			infectionCounts.put(k, infectionCounts.get(k) + infectionCount);

			// for debugging
			double h = root.getHeight();
			for (Node node : root.getAllLeafNodes()) {
				h = Math.min(h,  node.getHeight());
			}
			// System.err.println(toShortNewick(root, colourMap));
		
			if ((i+1) % 10 == 0) {
				if ((i+1) % 100 == 0) {
					System.err.print("|");
				} else {
					System.err.print(".");
				}
			}
			out.println(newick);
	    	traceout.print(i +"\t");
	    	

			if (printTaxaInTrace && taxonCount > 0) {
				TreeParser tree = new TreeParser(newick);
				IntegerParameter blockCount = new IntegerParameter();
				blockCount.initByName("dimension", taxonCount*2+1, "value", "-1");
				for (int j = 0; j < tree.getNodeCount(); j++) {
					Node node = tree.getNode(j);
					Object o = node.getMetaData("blockcount");
					if (o != null) {
						blockCount.setValue(j, (int)(double)o);
					}
				}
				int [] colourAtBase = new int[taxonCount*2-1];
				ColourProvider.getColour(tree.getRoot(), blockCount, tree.getNodeCount(), colourAtBase);
				
				int [] infectedBy = new int[taxonCount];
				Arrays.fill(infectedBy, -1);
				collectInfectedBy(root, infectedBy, taxonCount, colourAtBase, blockCount);
	    	
	    		for (int j = 0; j < taxonCount; j++) {
	    			traceout.printf(infectedBy[j] + "\t");
	    		}
			}
			double height = (root.getHeight() - h);
			double origin = endTimeInput.get().getArrayValue() - root.getHeight() + height;
			double logP2 = calcLogPInput.get() ? calcLogP(newick, h, origin) : 0;
	    	traceout.println((-h) + "\t" + height + "\t" + length(root) + "\t" + (endTimeInput.get().getArrayValue()-h) + "\t" + logP +
	    			(calcLogPInput.get() ? "\t" + logP2 : "" ));
			
	    	if (debug && Math.abs(logP-logP2) > 0.01) {
	    		out0.println(newick0);
	    	}
	    	
		}
		if (debug) {
			out0.close();
		}
		
		System.err.println();
		reportAttempts(taxonCounts, infectionCounts);
		
		if (traceOutputInput.get() != null && !traceOutputInput.get().getName().equals("[[none]]")) {
			traceout.close();
		}
		if (outputInput.get() != null) {
			out.close();
		}		
		Log.warning("Done");
	}
		
	private double calcLogP(String newick, double h, double origin) {
		TreeParser tree = new TreeParser(newick);

		int taxonCount = tree.getLeafNodeCount();
		
		IntegerParameter blockCount = new IntegerParameter();
		RealParameter blockStart = new RealParameter();
		RealParameter blockEnd = new RealParameter();
		blockCount.initByName("dimension", taxonCount*2-1, "value", "-1", "lower", -1, "upper", 1000);
		blockStart.initByName("dimension", taxonCount*2-2, "value", "0.5", "lower", 0.0, "upper", 1.0);
		blockEnd.initByName(  "dimension", taxonCount*2-2, "value", "0.5", "lower", 0.0, "upper", 1.0);
		
		for (int j = 0; j < tree.getNodeCount(); j++) {
			Node node = tree.getNode(j);
			Object o = node.getMetaData("blockcount");
			if (o != null) {
				blockCount.setValue(j, (int)(double)o);
			}
			o = node.getMetaData("blockstart");
			if (o != null) {
				blockStart.setValue(j, (double)o);
			}
			o = node.getMetaData("blockend");
			if (o != null) {
				blockEnd.setValue(j, (double)o);
			}
		}
		
		PopulationFunction.Abstract popFun = getPopFun();
		
		GammaHazardFunction transmissionHazard = new GammaHazardFunction();
		transmissionHazard.initByName("shape", transmissionShapeInput.get().getArrayValue() + "", 
				"rate", transmissionRateInput.get().getArrayValue() + "",
				"C", transmissionConstantInput.get().getArrayValue() + "");
		GammaHazardFunction sampleHazard = new GammaHazardFunction(); 
		sampleHazard.initByName("shape", sampleShapeInput.get().getArrayValue() + "", 
				"rate", sampleRateInput.get().getArrayValue() + "",
				"C", sampleConstantInput.get().getArrayValue() + "");
		
		TransmissionTreeLikelihood tl = new TransmissionTreeLikelihood();
		tl.initByName("tree", tree, 
				"blockstart", blockStart,
				"blockend", blockEnd,
				"blockcount", blockCount,
				"populationModel", popFun, 
				"endTime", h-endTimeInput.get().getArrayValue() + "",
				"origin", origin + "",
				"samplingHazard", sampleHazard,
				"transmissionHazard", transmissionHazard,
				"includeCoalescent", true
				);
		
		double logP = tl.calculateLogP();

		TransmissionTreeLikelihood1 tl1 = new TransmissionTreeLikelihood1();
		tl1.initByName("tree", tree, 
				"blockstart", blockStart,
				"blockend", blockEnd,
				"blockcount", blockCount,
				"populationModel", popFun, 
				"endTime", h-endTimeInput.get().getArrayValue() + "",
				"origin", origin + "",
				"samplingHazard", sampleHazard,
				"transmissionHazard", transmissionHazard,
				"includeCoalescent", true
				);
		
		logP = tl1.calculateLogP();
		return logP;
	}

	private Abstract getPopFun() {
		Abstract popFun;
		if (linearGrowthRateInput.get() != null) {
			popFun = new LinearGrowth();
			popFun.initByName("rate", linearGrowthRateInput.get().getArrayValue() + "");
		} else {
			popFun = new ConstantPopulation();
			popFun.initByName("popSize", popSizeInput.get().getArrayValue() + "");
		}
		return popFun;
	}

	private void reportAttempts(Map<Integer, Integer> taxonCounts, Map<Integer, Integer> infectionCounts) {
		DecimalFormat f = new DecimalFormat("##.###");
		Integer [] keys = taxonCounts.keySet().toArray(new Integer[] {});
		Arrays.sort(keys);
		long sum = 0;
		double mean = 0;
		if (maxTaxonCountInput.get() > 0) {
			Log.warning("#taxa\tpercentage of trees");
			for (int i : taxonCounts.keySet()) {
				if ( i <= maxTaxonCountInput.get()) {
					sum +=  taxonCounts.get(i);
				}
			}
			for (int i = 0; i <= maxTaxonCountInput.get(); i++) {
				if (taxonCounts.containsKey(i) && infectionCounts.containsKey(i)) {
					double percentage = 100.0*taxonCounts.get(i);
					Log.warning.print(i + "\t" + (percentage/sum < 10 ? " " : "") + f.format(((double)percentage)/sum));// + "\t" + taxonCounts.get(i));
					Log.warning.print("\t\t\t" + f.format((double)infectionCounts.get(i) / (i*taxonCounts.get(i))));
					Log.warning.println();
					mean += taxonCounts.get(i) * i;
				} else {
					Log.warning(i + "\t 0.000");
				}
			}
			int i = taxonCountInput.get() ;
			if (i > 0) {
				Log.warning.println("mean #transmissions per taxon: " + f.format((double)infectionCounts.get(i) / (i*treeCountInput.get())));
			}
		} else {
			Log.warning("#taxa\tpercentage of trees\tinfection count per taxon");
			for (int i : taxonCounts.values()) {
				sum +=  i;
			}
			for (Integer i : keys) {
				double percentage = 100.0*taxonCounts.get(i);
				Log.warning.print(i + "\t" + (percentage < 10 ? " " : "") + f.format(((double)percentage)/sum));
				Log.warning.println("\t\t\t" + f.format((double)infectionCounts.get(i) / (i*taxonCounts.get(i))));
				mean += taxonCounts.get(i) * i;
				sum +=  i;
			}
		}
		Log.warning("Mean #taxa = " + (mean / sum) + " based on " +sum + " observations");
	}

	private void collectInfectedBy(Node node, int[] infectedBy, int taxonCount, int [] colourAtBase, IntegerParameter blockCount) {
		if (!node.isRoot()) {
    		Node parent = node.getParent();
    		if (colourAtBase[node.getNr()] < taxonCount && colourAtBase[parent.getNr()] < taxonCount && 
    				colourAtBase[node.getNr()] != colourAtBase[parent.getNr()]) {
    			if (!directOnlyInput.get() || blockCount.getValue(node.getNr()) == 0) {
    				infectedBy[colourAtBase[node.getNr()]] = colourAtBase[parent.getNr()];
    			}
    		}
    	}
		for (Node child : node.getChildren()) {
			collectInfectedBy(child, infectedBy, taxonCount, colourAtBase, blockCount);
		}
	}

	private double length(Node node) {
		double length = 0;
		for (Node child : node.getChildren()) {
			length += length(child);
		}
		if (!node.isRoot()) {
			length += node.getLength();
		}
		return length;
	}

	private double runOnce(int maxTaxonCount) throws MathException {	
		double endTime = endTimeInput.get().getArrayValue();
		double popSize = popSizeInput.get().getArrayValue();

		double sampleShape = sampleShapeInput.get().getArrayValue();
		double sampleRate = sampleRateInput.get().getArrayValue();
		double sampleConstant = sampleConstantInput.get().getArrayValue();

		double transmissionShape = transmissionShapeInput.get().getArrayValue();
		double transmissionRate = transmissionRateInput.get().getArrayValue();
		double transmissionConstant = transmissionConstantInput.get().getArrayValue();

		root = new Node();
		root.setHeight(endTime);
		List<Node> nodes = new ArrayList<>();
		nodes.add(root);

		PoissonDistribution poisson = new PoissonDistributionImpl(transmissionConstant);
		GammaDistribution sampleIntensity = new GammaDistributionImpl(sampleShape, 1.0 / sampleRate);
		GammaDistribution transmissionIntensity = new GammaDistributionImpl(transmissionShape, 1.0 / transmissionRate);
		
		
		PopulationFunction popFun = getPopFun();

		List<Node> leafs = new ArrayList<>();
		int colour = 0;
		colourMap = new HashMap<>();
		colourMap.put(root, colour);
		
		this.logP = 0;
		while (nodes.size() > 0) {
			// continue with last node
			Node node = nodes.remove(nodes.size() - 1);

			// 1. draw number of events
			double r = Randomizer.nextDouble();
			int n = poisson.inverseCumulativeProbability(r);
			addToLogP("#events", Math.log(poisson.probability(n)));
					
			// 2. draw whether colour will be sampled
			boolean sample = (Randomizer.nextDouble() < sampleConstant);
			if (sample) {
				addToLogP("SampleP", Math.log(sampleConstant));
			} else {
				addToLogP("SampleP", Math.log(1.0 - sampleConstant));
			}
			
			
			// 3. simulate the time of sampling:
			// Note: do not need multiply by sampleConstant
			r = Randomizer.nextDouble();
			double sampletime = !sample ? 0
					: node.getHeight() - sampleIntensity.inverseCumulativeProbability(r);
			addToLogP("SampleTime", !sample ? 0 : sampleIntensity.logDensity(node.getHeight() - sampletime));
			if (sampletime < 0) {
				sample = false;
			}

			// 4. Simulate the times when node infects the new infectees
			double[] times = new double[n];
			for (int i = 0; i < n; i++) {
				// Note: do not need multiply by transmissionConstant
				r = Randomizer.nextDouble();
				times[i] = node.getHeight()
						- transmissionIntensity.inverseCumulativeProbability(r);
				addToLogP("TransTime", transmissionIntensity.logDensity(node.getHeight()-times[i]));
			}
			Arrays.sort(times);
			// remove times that are invalid: after study time, or after sample time (if any)
			while (n > 0 && (times[n - 1] < 0 || times[n - 1] < sampletime)) {
				n--;
			}
			
			List<Node> current = new ArrayList<>();
			// create leaf node
			if (sample) {
				Node leaf = new Node();
				colourMap.put(leaf, colour);
				leaf.setHeight(sampletime);
				leafs.add(leaf);
				current.add(leaf);
				leaf.setID("t" + format(leafs.size()));
				if (leafs.size() > maxTaxonCount) {
					if (!quietInput.get()) {
						System.err.print("x");
					}
					logP = Double.NEGATIVE_INFINITY;
					//runOnce(maxTaxonCount);
					addToLogP("\nrunOnce(too many taxa)", logP);
					return logP;
				}
			}
			// create internal (infection) nodes
			for (int i = 0; i < n; i++) {
				Node infectee = new Node();
				colourMap.put(infectee, colour);
				infectee.setHeight(times[i]);
				nodes.add(infectee);
				current.add(infectee);
			}

			
			// 5. Sampling a within-host phylogeny
			double currentHeight = sample ? sampletime : (n > 0 ? times[0] : 0);
			
			double [] logPCoalescent = new double[1];
			Node fragment = simulateCoalescent(current, popFun, currentHeight, node.getHeight(), maxAttemptsInput.get(), logPCoalescent);
			addToLogP("Coalescent:", logPCoalescent[0]);
			if (fragment == null) {
				if (!quietInput.get()) {
					System.err.print("c");
				}
				addToLogP("runOnce2", runOnce(maxTaxonCount));
				return logP;
			}
			// connect to node
			node.addChild(fragment);

			colourFragment(fragment, colour, colourMap);
			
			colour++;
		}
		
		
		// find all nodes to include in tree
		Set<Node> includedNodes = new HashSet<>();
		for (Node node : leafs) {
			while (node != root) {
				includedNodes.add(node);
				node = node.getParent();
			}
		}
		
		if (debug) {
			newick0 = root.toNewick();
		}
		// remove nodes not included from tree
		traverse(root, includedNodes);
		nodeCount = getNodeCount(root);
		return logP;
	}

	private void addToLogP(String caller, double log) {
//		System.err.println(caller + " " + log);
		logP += log;
	}

	private int getNodeCount(Node node) {
		int nodeCount = 1;
		for (Node child : node.getChildren()) {
			nodeCount += getNodeCount(child);
		}
		return nodeCount;
	}

	private String format(int i) {
		if (i >= 100) {
			return i + "";
		} else if (i >= 10) {
			return "0"+i;
		} else {
			return "00"+i;
		}
	}
	
	private Node simulateCoalescent(List<Node> current, PopulationFunction popFun, double currentHeight,
			double height, int maxAttemptCount, double [] logPCoalescent) {
		if (current.size() == 0) {
			Node node =  new Node();
			node.setHeight(currentHeight);
			return node;
		}
		
		List<Node> fragment;
		int attempt = 0;
		do {
			List<Node> currentCopy = new ArrayList<>(current);
			fragment = simulateCoalescent(currentCopy, popFun, currentHeight, height, logPCoalescent);
			if (fragment.size() == 1) {
				return fragment.get(0);
			}
			attempt++;
			logPCoalescent[0] = 0;
		} while (attempt < maxAttemptCount);
		Log.warning("Could not find a proper coalescent tree after " + maxAttemptCount + " attempts. "
				+ "Consider decreasing the population size or increasing maxAttempts.");
		return null;
	}

	private void colourFragment(Node node, int colour, Map<Node, Integer> colourMap) {
		colourMap.put(node,  colour);
		for (Node child : node.getChildren()) {
			colourFragment(child, colour, colourMap);
		}
	}

	private String toShortNewick(Node node, Map<Node, Integer> colourMap) {
		if (node.isLeaf()) {
			return node.getID();
		} else {
			String newick = "(";
			for (int i = 0; i < node.getChildCount(); i++) {
				Node child = node.getChild(i);
				newick += toShortNewick(child, colourMap);
				newick += "[&color=" + colourMap.get(child)+"]" + ":" + node.getLength();
				if (i < node.getChildCount() - 1) {
					newick += ",";
				}
			}
			return newick + ")";
		}
	}
	
	private int infectionCount(Node node) {
		switch(node.getChildCount()) {
		case 0: {// leaf
			Node p = node.getParent();
			int blockCount = -1;
			while (p != null && p.getChildCount() == 1) {
				p = p.getParent();
				if (colourMap.get(node) != colourMap.get(p)) {
					blockCount++;
				}
			}
			return blockCount+1;
		}
		case 1:
			return infectionCount(node.getChild(0));
		case 2:
			Node left = node.getLeft();
			int infectionCount = infectionCount(left);
			Node right = node.getRight();
			infectionCount += infectionCount(right);

			Node p = node.getParent();
			int blockCount = -1;
			while (p != null && p.getChildCount() == 1) {
				p = p.getParent();
				if (colourMap.get(node) != colourMap.get(p)) {
					blockCount++;
				}
			}
			return infectionCount + blockCount + 1;
		}
		return 0;		
	}
	
	private String toNewick(Node node) {
		switch(node.getChildCount()) {
		case 0: {// leaf
			double length = node.getLength();
			Node p = node.getParent();
			double blockStart = length;
			double blockEnd = length;
			int blockCount = -1;
			while (p != null && p.getChildCount() == 1) {
				blockEnd = length;
				length += p.getLength();
				p = p.getParent();
				if (colourMap.get(node) != colourMap.get(p)) {
					blockCount++;
				}
			}
			return node.getID() + "[&blockcount=" + blockCount + (blockCount >= 0 ? ",blockstart=" + (blockStart/length) + ",blockend=" + (blockEnd/length): "") +",color=" + colourMap.get(node) + "]:" + length;
		}
		case 1:
			return toNewick(node.getChild(0));
		case 2:
			Node left = node.getLeft();
			String leftNewick = toNewick(left);
			Node right = node.getRight();
			String rightNewick = toNewick(right);

			double length = node.getLength();
			Node p = node.getParent();
			double blockStart = length;
			double blockEnd = length;
			int blockCount = -1;
			while (p != null && p.getChildCount() == 1) {
				blockEnd = length;
				length += p.getLength();
				p = p.getParent();
				if (colourMap.get(node) != colourMap.get(p)) {
					blockCount++;
				}
			}
			if (p == null) {
				return "(" + leftNewick + "," + rightNewick + ")";
			}
			return "(" + leftNewick + "," + rightNewick + ")" + "[&blockcount=" + blockCount + (blockCount >= 0 ? ",blockstart=" + (blockStart/length) + ",blockend=" + (blockEnd/length): "") +",color=" + colourMap.get(node) + "]:" + length; 
		}
		return null;
	}

	private void traverse(Node node, Set<Node> includedNodes) {
		List<Node> children = node.getChildren();
		for (int i = children.size()-1; i >= 0; i--) {
			Node child = children.get(i);
			if (!includedNodes.contains(child)) {
				node.removeChild(child);
			}
		}
		for (int i = node.getChildren().size()-1; i >= 0; i--) {
			Node child = children.get(i);
			traverse(child, includedNodes);
		}		
	}
	
	public List<Node> simulateCoalescent(final List<Node> nodes, final PopulationFunction demographic,
			double currentHeight, final double maxHeight, double [] logPCoalescent) {
		// If only one node, return it
		// continuing results in an infinite loop
		if (nodes.size() == 1)
			return nodes;

		final double[] heights = new double[nodes.size()];
		for (int i = 0; i < nodes.size(); i++) {
			heights[i] = nodes.get(i).getHeight();
		}
		final int[] indices = new int[nodes.size()];
		HeapSort.sort(heights, indices);

		// node list
		List<Node> nodeList = new ArrayList<>();
		int activeNodeCount = 0;
		for (int i = 0; i < nodes.size(); i++) {
			nodeList.add(nodes.get(indices[i]));
		}
        while (getMinimumInactiveHeight(activeNodeCount, nodeList) <= currentHeight) {
            activeNodeCount += 1;
        }

		// get at least two tips
		while (activeNodeCount < 2) {
			currentHeight = getMinimumInactiveHeight(activeNodeCount, nodeList);
	        while (getMinimumInactiveHeight(activeNodeCount, nodeList) <= currentHeight) {
	            activeNodeCount += 1;
	        }
		}

		// simulate coalescent events
		double r = Randomizer.nextDouble();
		double lambda = demographic instanceof LinearGrowth ? ((LinearGrowth)demographic).getRate() : 0;
		double nextCoalescentHeight = currentHeight
				+ ( demographic instanceof LinearGrowth ?
//					+ (1.0/lambda)*(Math.pow(1-r, lambda/(activeNodeCount * (activeNodeCount-1)/2))) * lambda * (maxHeight - currentHeight);
					+ (Math.pow(1-r, lambda/(activeNodeCount * (activeNodeCount-1)/2))) * (maxHeight - currentHeight)
				: PopulationFunction.Utils.getSimulatedInterval(demographic, activeNodeCount, currentHeight));
        
        double kChoose2 = Binomial.choose2(activeNodeCount);
		double intervalArea = demographic.getIntegral(currentHeight, nextCoalescentHeight);
		logPCoalescent[0] += -kChoose2 * intervalArea;

		// while (nextCoalescentHeight < maxHeight && (getNodeCount() > 1)) {
		while (nextCoalescentHeight < maxHeight && (nodeList.size() > 1)) {

			if (nextCoalescentHeight >= getMinimumInactiveHeight(activeNodeCount, nodeList)) {
				currentHeight = getMinimumInactiveHeight(activeNodeCount, nodeList);
		        while (getMinimumInactiveHeight(activeNodeCount, nodeList) <= currentHeight) {
		            activeNodeCount += 1;
		        }
			} else {
				currentHeight = coalesceTwoActiveNodes(currentHeight, nextCoalescentHeight, nodeList, activeNodeCount, logPCoalescent);
				activeNodeCount--;
			}

			
			if (!(demographic instanceof LinearGrowth)) {
				final double demographicAtCoalPoint = demographic.getPopSize(currentHeight);
				logPCoalescent[0] -= Math.log(demographicAtCoalPoint);
			}

                // if (getNodeCount() > 1) {
			if (nodeList.size() > 1) {
				// get at least two tips
				while (activeNodeCount < 2) {
					currentHeight = getMinimumInactiveHeight(activeNodeCount, nodeList);
			        while (getMinimumInactiveHeight(activeNodeCount, nodeList) <= currentHeight) {
			            activeNodeCount += 1;
			        }
				}

				// nextCoalescentHeight = currentHeight +
				// DemographicFunction.Utils.getMedianInterval(demographic,
				// getActiveNodeCount(), currentHeight);
				r = Randomizer.nextDouble();
				nextCoalescentHeight = currentHeight
						+ ( demographic instanceof LinearGrowth ?
//							+ (1.0/lambda)*(Math.pow(1-r, lambda/(activeNodeCount * (activeNodeCount-1)/2))) * lambda * (maxHeight - currentHeight);
							+ (Math.pow(1-r, lambda/(activeNodeCount * (activeNodeCount-1)/2))) * (maxHeight - currentHeight)
						: PopulationFunction.Utils.getSimulatedInterval(demographic, activeNodeCount, currentHeight));

				kChoose2 = Binomial.choose2(activeNodeCount);
				intervalArea = demographic.getIntegral(currentHeight, nextCoalescentHeight);
				logPCoalescent[0] += -kChoose2 * intervalArea;
			}
		}

		return nodeList;
	}
	
    /**
     * Coalesce two nodes in the active list. This method removes the two
     * (randomly selected) active nodes and replaces them with the new node at
     * the top of the active list.
     * @param minHeight
     * @param height
     * @return
     */
    private double coalesceTwoActiveNodes(final double minHeight, double height, List<Node> nodeList, int activeNodeCount, double[]logPCoalescent){
        final int node1 = Randomizer.nextInt(activeNodeCount);
        int node2 = node1;
        while (node2 == node1) {
            node2 = Randomizer.nextInt(activeNodeCount);
        }
        logPCoalescent[0] += -Math.log(activeNodeCount * (activeNodeCount-1)/2);

        final Node left = nodeList.get(node1);
        final Node right = nodeList.get(node2);

        final Node newNode = new Node();
//		System.err.println(2 * m_taxa.get().getNrTaxa() - nodeList.size());
//        newNode.setNr(nextNodeNr++);   // multiple tries may generate an excess of nodes assert(nextNodeNr <= nrOfTaxa*2-1);
        newNode.setHeight(height);
        newNode.setLeft(left);
        left.setParent(newNode);
        newNode.setRight(right);
        right.setParent(newNode);

        nodeList.remove(left);
        nodeList.remove(right);

        activeNodeCount -= 2;

        nodeList.add(activeNodeCount, newNode);

        activeNodeCount += 1;

        
        if (getMinimumInactiveHeight(activeNodeCount, nodeList) < height) {
            throw new RuntimeException(
                    "This should never happen! Somehow the current active node is older than the next inactive node!\n"
            		+ "One possible solution you can try is to increase the population size of the population model.");
        }
        return height;
    }
    
    private double getMinimumInactiveHeight(int activeNodeCount, List<Node> nodeList) {
        if (activeNodeCount < nodeList.size()) {
            return (nodeList.get(activeNodeCount)).getHeight();
        } else
            return Double.POSITIVE_INFINITY;
    }


	public static void main(String[] args) throws Exception {
		new Application(new TransmissionTreeSimulator(), "SimulatedTransmissionTree", args);
	}

}
