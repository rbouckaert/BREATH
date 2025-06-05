package breath.util;


import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.imageio.ImageIO;

import org.graphper.api.Graphviz;
import org.graphper.api.Graphviz.GraphvizBuilder;
import org.graphper.api.Line;
import org.graphper.api.Node;
import org.graphper.api.attributes.Color;

import beastfx.app.tools.Application;
import beastfx.app.util.LogFile;
import beastfx.app.util.OutFile;
import beastfx.app.util.TreeFile;
import beastfx.app.util.Utils;
import breath.distribution.ColourProvider;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.util.HeapSort;
import beastfx.app.tools.LogAnalyser;
import beastfx.app.treeannotator.TreeAnnotator;
import beastfx.app.treeannotator.TreeAnnotator.MemoryFriendlyTreeSet;

@Description("Create SVG files to visualise who infected who")
public class WIWVisualiser extends beast.base.inference.Runnable {
	final public Input<TreeFile> treeFile = new Input<>("trees", "tree file file with transmission trees.", new TreeFile("[[none]]"));
	final public Input<LogFile> inFile = new Input<>("log", "trace file containing infectorOf log. Ignored if tree file is specified", new LogFile("[[none]]"));
	final public Input<Integer> burnInPercentageInput = new Input<>("burnin", "percentage of trees to used as burn-in (and will be ignored)", 10);
	final public Input<OutFile> outputInput = new Input<>("out", "output file, or stdout if not specified",
			new OutFile("/tmp/wiw.svg"));
	final public Input<OutFile> matrixOutputInput = new Input<>("matrix", "transition probability matrix output file, potentially useful for post-processing e.g. visualising as heat map.. Ignored if not specified.");
	final public Input<String> prefixInput = new Input<>("prefix", "prefix of infectorOf entry, e.g., infectorOf", "infectorOf");
	final public Input<Double> thresholdInput = new Input<>("threshold", "probability threshold below which edges will be ignored.", 0.1);
	final public Input<String> partitionInput = new Input<>("partition", "name of the partition appended to `blockcount, blockend and blockstart`");
	final public Input<Boolean> suppressSingletonInput = new Input<>("suppressSingleton", "do not show taxa that are not connected to any other taxa", true);
	final public Input<Boolean> colourByAgeInput = new Input<>("colourByAge", "colour nodes in output by node age. All blacks if false", true);
	final public Input<Boolean> widthByPosteriorInput = new Input<>("widthByPosterior", "draw line between nodes with widths proportional to posterior support", true);

	final public Input<Float> saturationInput = new Input<>("saturation", "saturation used when colouring nodes.", 0.7f);
	final public Input<Float> brightnessInput = new Input<>("brightness", "brightness used when colouring nodes.", 0.7f);
	final public Input<String> filterInput = new Input<>("filter", "search/replace regular expression for filtering labels. Should be of the form '/searchRegExp/replaceString/'. Ignored if not specified");

	final static String DIR_SEPARATOR = (Utils.isWindows() ? "\\\\" : "/");

	
	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		String [] nodeLabels;
		double [] age = null;
		int n;
		double upper = 0;
		
		String search = null, replace = null;
		if (filterInput.get() != null) {
			String filter = filterInput.get();
			int sep = filter.indexOf('/', 1);
			search = filter.substring(1, sep);
			replace = filter.substring(sep + 1, filter.length() - 1);
		}

		double [][] transitions;
		if (treeFile.get() != null && !treeFile.get().getName().equals("[[none]]")) {
	        MemoryFriendlyTreeSet trees = new TreeAnnotator().new MemoryFriendlyTreeSet(treeFile.get().getAbsolutePath(), burnInPercentageInput.get());
	        trees.reset();
	    	Tree tree = trees.next();
	    	n = tree.getLeafNodeCount();
	    	
			if (colourByAgeInput.get()) {
				age = new double[n];
				for (int i = 0; i < n; i++) {
					age[i] = tree.getNode(i).getHeight();
					upper = Math.max(upper, age[i]);
				}
			}

			nodeLabels = new String[n];
			for (int i = 0; i < n; i++) {
				nodeLabels[i] = tree.getNode(i).getID();
				if (search != null) {
					nodeLabels[i] = nodeLabels[i].replaceAll(search, replace);
				}
			}
			
	        trees.reset();
			transitions = new double[n][n+1];
			int sampleCount = 0;
	        while (trees.hasNext()) {
	        	tree = trees.next();
	        	
	        	// extract meta data from tree
	        	Integer [] count = new Integer[n*2-1];
	        	for (int i = 0; i < tree.getNodeCount(); i++) {
	        		beast.base.evolution.tree.Node node = tree.getNode(i);
	        		Object o = node.getMetaData("blockcount");
	        		if (o == null) {
	        			o = node.getMetaData("blockcount.t:" + partitionInput.get());
	        		}
	        		count[i] = o == null ? 0 : (int)(double) o;
	        	}
	            IntegerParameter blockCount = new IntegerParameter(count);

	            // calculate colouring
	        	int [] colourAtBase = new int[n*2-1];
	        	ColourProvider.getColour(tree.getRoot(), blockCount, tree.getLeafNodeCount(), colourAtBase);
	        	
	        	// determine who infected who
	        	int [] infectedBy = new int[tree.getLeafNodeCount()];
	        	Arrays.fill(infectedBy, -1);
	        	for (int i = 0; i < 2 * n - 2; i++) {
	        		beast.base.evolution.tree.Node node = tree.getNode(i);
	        		beast.base.evolution.tree.Node parent = node.getParent();
	        		if (colourAtBase[node.getNr()] < n && colourAtBase[parent.getNr()] < n && 
	        				colourAtBase[node.getNr()] != colourAtBase[parent.getNr()]) {
	        			if (blockCount.getValue(node.getNr()) == 0) {
	        				infectedBy[colourAtBase[node.getNr()]] = colourAtBase[parent.getNr()];
	        			}
	        		}
	        	}
	        	
	        	
	        	// log the result
	        	for (int i = 0; i < n; i++) {
		        	transitions[i][(infectedBy[i] + n + 1) % (n+1)]++;
	        	}
	        	sampleCount++;
	        }			
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n+1; j++) {
					transitions[i][j] /= sampleCount;
				}
			}
		} else {
			LogAnalyser trace = new LogAnalyser(inFile.get().getPath(), burnInPercentageInput.get(), true, false);
			
			// count number of taxa from trace file
			List<String> labels = trace.getLabels();
			n = 0;
			String prefix = prefixInput.get() + ".1";
			int offset = 0;
			for (int i = 0; i < labels.size(); i++) {
				String label = labels.get(i);
				if (label.equals(prefix)) {
					offset = i;
					prefix = prefixInput.get();
					while (labels.get(i+n).equals(prefix + "." + (n+1))) {
						n++;
					}
					break;
				}
			}
			
			// count transitions
			transitions = new double[n][n+1];
			int sampleCount = 0;
			for (int i = 0; i < n; i++) {
				Double [] infectorOf = trace.getTrace(offset + i + 1);
				sampleCount = infectorOf.length;
				for (Double d : infectorOf) {
					transitions[i][ (int)(d + n + 1) % (n+1)]++;
				}
			}
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n+1; j++) {
					transitions[i][j] /= sampleCount;
				}
			}
			nodeLabels = new String[n];
			for (int i = 0; i < n; i++) {
				nodeLabels[i] = (i+1) + "";
			}
		}		
		

		
		// build transition graph
		GraphvizBuilder dotty = Graphviz.digraph();
		//dotty = dotty.tempNode(Node.builder().shape(NodeShapeEnum.RECT).build());
		
		// add nodes
		DecimalFormat f = new DecimalFormat("#.##");
		Node [] nodes = new Node[n];

		boolean [] nodesInUse = new boolean[n];
		double threshold = thresholdInput.get();
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i!=j && transitions[i][j] >= threshold) {
					nodesInUse[i] = true;
					nodesInUse[j] = true;
				}
			}
		}

		
		
		
		int suppressedCount = 0;
		for (int i = 0; i < n; i++) {
			if (nodesInUse[i] || !suppressSingletonInput.get()) {
				String colourString;
				Node node;
				if (age != null) {
					int c = java.awt.Color.HSBtoRGB((float)(age[i]/upper), saturationInput.get(), brightnessInput.get());
					colourString = "#" + Integer.toHexString(c).substring(2);
					Color colour = Color.ofRGB(colourString);
					c = java.awt.Color.HSBtoRGB((float)(age[i]/upper), 0.5f, 0.9f);
					colourString = "#" + Integer.toHexString(c).substring(2);
					Color fillcolour = Color.ofRGB(colourString);
					node = Node.builder()
							.label(nodeLabels[i] + " (" + f.format(1.0-transitions[i][n]-transitions[i][i]) + ")")
							.color(colour)
							.fillColor(fillcolour)
							.build();
				} else {
					node = Node.builder()
						.label(nodeLabels[i] + " (" + f.format(1.0-transitions[i][n]-transitions[i][i]) + ")")
						.build();
				}
	//			System.err.println(transitions[i][i]);
				dotty = dotty.addNode(node);
				nodes[i] = node;
			} else {
				suppressedCount++;
			}
		}

		if (suppressSingletonInput.get()) {
			if (suppressedCount > 0) {
				Log.warning( "[Info] " + suppressedCount + " " + (suppressedCount == 1 ? "node" : "nodes") +
						" not shown due to being " + (suppressedCount == 1 ? "a singleton" : "singletons") + ".");
			}
			if (suppressedCount == n) {
				Log.warning(
						       "[Warning] No who-infected-who network will be created - all nodes have been removed.\n"
								+ "          Possible remedies:\n"
								+ "            o Ensure the partition name is correct.\n"
								+ "            o Try reducing the threshold value.\n"
								+ "            o Check if the hazard function parameters scale well with the tree size.\n"
				);
				return;
			}
		}

		// add edges
		double pen;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i!=j && transitions[i][j] >= threshold) {
					if (widthByPosteriorInput.get()) {
						pen = transitions[i][j]*15;
						dotty = dotty.addLine(Line.builder(nodes[j], nodes[i]).penWidth(pen)
								.label(f.format(transitions[i][j]))
								.build());
					} else if (colourByAgeInput.get()) {
						int c = java.awt.Color.HSBtoRGB((float)(transitions[i][j] * 0.5 + 0.5), saturationInput.get(), brightnessInput.get());
						String colourString = "#" + Integer.toHexString(c).substring(2);
						Color fillcolour = Color.ofRGB(colourString);
						dotty = dotty.addLine(Line.builder(nodes[j], nodes[i])
								.penWidth(2.0)
								.color(fillcolour)
								.fontColor(fillcolour)
								.label(f.format(transitions[i][j]))
								.build());
					} else {
						dotty = dotty.addLine(Line.builder(nodes[j], nodes[i])
								.label(f.format(transitions[i][j]))
								.build());
					}
				}
			}
		}
		
		
		// System.out.println(dotty.build().toString());
		// create SVG
		String svg = dotty.build().toSvgStr();

		// output graph
		try {
			File svgFile = outputInput.get();
			Log.warning("Writing network to file: " + svgFile.getPath());
			FileWriter outfile = new FileWriter(svgFile);
			outfile.write(svg);
			outfile.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		if (colourByAgeInput.get()) {
			String parent = outputInput.get().getParent();
			String legendPath = (parent != null ? parent + DIR_SEPARATOR : "") + "legend.png";
			Log.warning("Writing legend to file: " + legendPath);

			BufferedImage image = new BufferedImage(200, 200, BufferedImage.TYPE_INT_RGB);
	        Graphics g = image.getGraphics();
	        g.setColor(java.awt.Color.white);
	        g.fillRect(0, 0, 200, 200);
	        for (int i = 0; i < 200; i++) {
	        	double x = (200 - i) / 200.0;
				g.setColor(java.awt.Color.getHSBColor((float) x, saturationInput.get(), brightnessInput.get()));
	        	g.drawLine(25, i, 75, i);
	        }
			g.setColor(java.awt.Color.black);
	        for (int i = 0; i < 10; i++) {
	        	g.drawString(f.format(upper * (10 - i)/ 10), 100, 200*i/10);
	        }
			ImageIO.write(image, "png", new File(legendPath));
		}
		outputMatrix(nodeLabels, transitions);
	}

	/*
	 * output transition matrix to tab separated file
	 */
	private void outputMatrix(String[] nodeLabels, double[][] transitions) throws FileNotFoundException {
		if (matrixOutputInput.get() != null) {
			System.err.println("Writing transition matrix " + matrixOutputInput.get().getPath());
			int [] order = new int[nodeLabels.length];
			
			// labels are numeric?
			boolean isNumeric = true;
			for (String str : nodeLabels) {
				try {
					int h = Integer.valueOf(str);
				} catch (NumberFormatException e) {
					isNumeric = false;
					break;
				}
			}
			if (isNumeric) {
				List<Integer> labels = new ArrayList<>();
				for (String label : nodeLabels) {
					labels.add(Integer.valueOf(label));
				}
				HeapSort.sort(labels, order);

			} else {
				List<String> labels = new ArrayList<>();
				for (String label : nodeLabels) {
					labels.add(label);
				}
				HeapSort.sort(labels, order);
			}

			PrintStream out = new PrintStream(matrixOutputInput.get());
			out.print("\t");
			for (int i = 0; i < nodeLabels.length; i++) {
				out.print(nodeLabels[order[i]] + "\t");
			}
			out.println();
			for (int i = 0; i < nodeLabels.length; i++) {
				out.print(nodeLabels[order[i]] + "\t");
				for (int j = 0; j < nodeLabels.length; j++) {
					out.print(transitions[order[i]][order[j]] + "\t");
				}
				out.println();
			}

			out.close();
		}
	}

	public static void main(String[] args) throws Exception {
		new Application(new WIWVisualiser(), "Who-Infected-Who Visualiser", args);
	}
}
