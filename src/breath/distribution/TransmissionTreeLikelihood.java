package breath.distribution;




import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import breath.evolution.LinearGrowth;
import org.apache.commons.math.MathException;
import org.apache.commons.math.special.Gamma;
import org.apache.commons.math3.util.FastMath;

import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.evolution.tree.IntervalList;
import beast.base.evolution.tree.IntervalType;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeDistribution;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.State;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Binomial;

@Description("Likelihood of a transmission tree")
@Citation(value="Caroline Colijn, Matthew David Hall, Remco Bouckaert.\n"
		+ "Taking a BREATH (Bayesian Reconstruction and Evolutionary Analysis of Transmission Histories) to simultaneously infer phylogenetic and transmission trees for partially sampled outbreaks.\n"
		+ "BioRxiv, 2024", year=2024, DOI="10.1101/2024.07.11.603095)")
public class TransmissionTreeLikelihood extends TreeDistribution {
	final public Input<RealParameter> blockStartFractionInput = new Input<>("blockstart", "start of block in fraction of branch length", Validate.REQUIRED);
	final public Input<RealParameter> blockEndFractionInput = new Input<>("blockend", "end of block in fraction of branch length", Validate.REQUIRED);
	final public Input<IntegerParameter> blockCountInput = new Input<>("blockcount", "number of transitions inside a block", Validate.REQUIRED);
	//    final public Input<IntegerParameter> colourInput = new Input<>("colour", "colour of the base of the branch", Validate.REQUIRED);
	final public Input<PopulationFunction> popSizeInput = new Input<>("populationModel", "A population size model", Validate.REQUIRED);

	final public Input<Function> originInput = new Input<>("origin", "time at which the study start above the root of tree. Assumed to be at root if not specified");
	final public Input<RealParameter> endTimeInput = new Input<>("endTime", "time at which the study finished", Validate.REQUIRED);
	//final public Input<RealParameter> lambdaTrInput = new Input<>("lambda", "lambda parameter of Poisson process", Validate.REQUIRED);

	final public Input<GammaHazardFunction> samplingHazardInput = new Input<>("samplingHazard", "determines the hazard of being sampled", Validate.REQUIRED);
	final public Input<GammaHazardFunction> transmissionHazardInput = new Input<>("transmissionHazard", "determines the hazard of transmitting an infection", Validate.REQUIRED);

	final public Input<Boolean> colourOnlyInput = new Input<>("colourOnly", "flag for debugging that calculates colour at base only, but does not contribute to posterior otherwise", false);
	final public Input<Boolean> includeCoalescentInput = new Input<>("includeCoalescent", "flag for debugging that includes contribution from coalescent to posterior if true", true);
	final public Input<Boolean> conditionOnInfectionTimeInput = new Input<>("conditionOnInfectionTime", "flag whether to condition the coalescent for subtrees on infection time. "
			+ "If true, population size is assumed to be constant within a host. "
			+ "If false, any population function can be used.", true);

	final public Input<Boolean> allowTransmissionsAfterSamplingInput = new Input<>("allowTransmissionsAfterSampling", "flag to indicate sampling does not affect the probability of onwards transmissions. "
			+ "If false, no onwards transmissions are allowed (not clear how this affects the unknown unknowns though).", true);

	final public Input<Double> branchLengthThresholdInput = new Input<>("branchLengthThreshold", "minimal branch length for which penalty applies (to prevent very samll branch lengths)", 1e-4);


	private Tree tree;
	private RealParameter blockStartFraction;
	private RealParameter blockEndFraction;
	private IntegerParameter blockCount;
	private int [] colourAtBase;
	private PopulationFunction popSizeFunction;
	private Validator validator;
	private Function origin;
	private double branchLengthThreshold;

	// hazard functions for sampling and transmission respectively

	private RealParameter endTime; // end time of study
	private double lambda;
	private GammaHazardFunction samplingHazard;
	private GammaHazardFunction transmissionHazard;

	private double Cs;
	private double Ctr;
	private double p0;
	private double phi;
	private double rho;
	private double atr;
	private double btr;

	//private double a, b;

	private boolean updateColours = true;
	private boolean allowTransmissionsAfterSampling;
	private boolean initialCalculation = true;
	private boolean conditionOnInfectionTime = true;

	@Override
	public void initAndValidate() {
		tree = (Tree) treeInput.get();
		if (tree == null) {
			tree = treeIntervalsInput.get().treeInput.get();
		}

		int n = tree.getNodeCount();
		blockStartFraction = blockStartFractionInput.get();
		blockEndFraction = blockEndFractionInput.get();
		blockCount = blockCountInput.get();
		colourAtBase = new int[n];

		sanityCheck(blockStartFraction, n-1 , "blockStart");
		sanityCheck(blockEndFraction, n-1, "blockEnd");

		if (blockCount.getDimension() != n) {
			blockCount.setDimension(n);
			// no transmission at root
			blockCount.setValue(n-1, -1);
			Log.warning("WARNING: Setting dimension of blockCount parameter " + blockCount.getID() + " to " + n);
		}
		if (blockCount.getLower() < -1) {
			blockCount.setBounds(-1, blockCount.getUpper());
			Log.warning("WARNING: Setting lower bound of blockCount parameter " + blockCount.getID() + " to -1");
		}

		popSizeFunction = popSizeInput.get();

		validator = new Validator(tree, colourAtBase, blockCount, blockStartFraction, blockEndFraction);

		origin = originInput.get();



		endTime = endTimeInput.get();
		// lambda_tr = lambdaTrInput.get();
		samplingHazard = samplingHazardInput.get();
		transmissionHazard = transmissionHazardInput.get();

		Cs = samplingHazard.constantInput.get().getArrayValue();
		Ctr = transmissionHazard.constantInput.get().getArrayValue();
		atr = transmissionHazard.shapeInput.get().getArrayValue();
		btr = transmissionHazard.getRate();

		double f = getRetainedFrac(50000);
		lambda = (Cs*f*Ctr + (1-Cs)*Ctr) ;

		p0 = getp0(Cs, lambda, 0.1);
		phi = getPhi(Cs, lambda, p0);
		rho = getRho(phi);

		allowTransmissionsAfterSampling = allowTransmissionsAfterSamplingInput.get();
		conditionOnInfectionTime = conditionOnInfectionTimeInput.get();
		branchLengthThreshold = branchLengthThresholdInput.get();
	}

	private double getRetainedFrac(int numSamps) {
		double f = 0;

		double mean = 0;
		double x2 = 0;

		for (int i = 0; i < numSamps; i++) {
			double tInf;
			try {
				tInf = transmissionHazard.simulate();
				double tSam = samplingHazard.simulate();
				if ( tInf < tSam) {
					f++;
				} else {
					mean += tInf;
					x2 += tInf * tInf;
				}
			} catch (MathException e) {
				e.printStackTrace();
			}
		}
		mean /= (numSamps-f);
		x2 /= (numSamps-f);
		double var = x2 - mean*mean;
//    	b = mean / var;
//    	a = mean * b;

		return f / numSamps;
	}

	private void sanityCheck(RealParameter blockFraction, int n, String paramName) {
		if (blockFraction.getDimension() != n) {
			blockFraction.setDimension(n);
			Log.warning("WARNING: Setting dimension of "  + paramName + " parameter " + blockFraction.getID() + " to " + n);
		}
		if (blockFraction.getLower() < 0) {
			blockFraction.setLower(0.0);
			Log.warning("WARNING: Setting lower bound of "  + paramName + " parameter " + blockFraction.getID() + " to 0");
		}
		if (blockFraction.getUpper() > 1) {
			blockFraction.setUpper(1.0);
			Log.warning("WARNING: Setting upper bound of  "  + paramName + " parameter " + blockFraction.getID() + " to 1");
		}
	}


	@Override
	public double calculateLogP() {
		if (initialCalculation) {
			double height = tree.getRoot().getHeight();
			if (height > origin.getArrayValue()) {
				Log.warning("\nWARNING: origin (" + origin.getArrayValue() +") is less than tree height (" + height +"). Consider increasing the origin start value to get the chain started.\n");
			}
			if (origin.getArrayValue() - height > 75) {
				Log.warning("\nWARNING: origin (" + origin.getArrayValue() +") is more than 75 over tree height (" + height +"). Consider decreasing the origin start value to prevent numerical issues.\n");
			}
		}
		initialCalculation = false;

		logP = 0;

		if (origin.getArrayValue() < tree.getRoot().getHeight()) {
			logP = Double.NEGATIVE_INFINITY;
			return logP;
		}

		if (!calcColourAtBase()) {
			logP = Double.NEGATIVE_INFINITY;
			return logP;
		}

		if (!validator.isValid(colourAtBase)) {
			logP = Double.NEGATIVE_INFINITY;
			return logP;
		}

		if (branchLengthThreshold > 0) {
			for (Node node : tree.getNodesAsArray()) {
				if (node.getLength() < branchLengthThreshold && !node.isRoot()) {
					logP += -10000;
				}
			}
		}

		if (colourOnlyInput.get()) {
			return logP;
		}

		segments = collectSegments();

		if (includeCoalescentInput.get()) {
			logP += calculateCoalescent();
		}

		logP += calcTransmissionLikelihood();
		if (Double.isInfinite(logP)) {
			logP = Double.NEGATIVE_INFINITY;
		}
		return logP;
	}




	// initialise colourAtBase
	// return true if a valid colouring can be found, 
	// return false if there is a path between leafs without a transmission
	public boolean calcColourAtBase() {
		updateColours = ColourProvider.getColour(tree.getRoot(), blockCount, tree.getLeafNodeCount(), colourAtBase);
		return updateColours;
	}

	public double calcTransmissionLikelihood() {
		double d = endTime.getArrayValue();
		double logP = 0;
		int n = tree.getLeafNodeCount();
		Node [] nodes = tree.getNodesAsArray();
		if (segments == null) {
			segments = collectSegments();
		}

		if (origin != null) {
			if (origin.getArrayValue() < tree.getRoot().getHeight()) {
				return Double.NEGATIVE_INFINITY;
			}
		}

//System.err.println("\n#contribution of sampled cases");
		// contribution of sampled cases
		for (int i = 0; i < n; i++) {
//System.err.println("#node " + (i+1));
			double logP1 = 0;
			// contribution of not being sampled
			SegmentIntervalList intervals = segments.get(i);
			double start = intervals.birthTime;
			double end = intervals.times.get(0);
			logP1 += logh_s(start, end) + logS_s(start, end);
			// contribution of causing infections
			if (allowTransmissionsAfterSampling) {
				logP1 +=  logS_tr(start, d); // further contribution below
			} else {
				logP1 +=  logS_tr(start, end); // further contribution below
			}
			logP1 -= logGetIndivCondition(p0, start, d);
			if (Double.isInfinite(logP1) && logP1 > 0) {
				System.err.println("Numerical instability encountered: ");
				System.err.println(start + " " + d + " " + end + " " + p0);
				System.err.println(logS_tr(start, d));
				System.err.println(logS_tr(start, end));
				System.err.println(logGetIndivCondition(p0, start, d));
			}
//System.err.println("#node " + (i+1) + " " + logP1);
			logP += logP1;
		}

//System.err.println("\n#transmissions from sampled cases");
		// further contribution of causing infections
		for (int i = 0; i < tree.getNodeCount() - 1; i++) {
			int baseColour = colourAtBase[i];
			int parent = nodes[i].getParent().getNr();
			int parentColour = colourAtBase[parent];
			if (baseColour != parentColour && parentColour < n) {
//System.err.println("#node " + (i+1));
				double tInf0 = segments.get(parentColour).birthTime;
				Node node = nodes[i];
				double tInf1 = node.getHeight() + node.getLength() * blockEndFraction.getArrayValue(node.getNr());
				double logP1 = logh_tr(tInf0, tInf1);
				logP += logP1;
			}
		}

//System.err.println("\n#contribution of unsampled cases");
		// contribution of unsampled cases
		for (int i = n; i < tree.getNodeCount(); i++) {
			if (colourAtBase[i] >= n) {
//System.err.println("#node " + (i+1));

				// contribution of not being sampled
				SegmentIntervalList intervals = segments.get(i);
				if (intervals != null) {
					double start = intervals.birthTime;
					Double logP1 = logS_s(start, d);
					// contribution of causing infections
					logP += logS_tr(start, d); // further contribution below
					logP -= logGetIndivCondition(p0, start, d);

					logP += logP1;
				}
			}
		}

//System.err.println("\n#transmissions from unsampled cases");
		// further contribution of causing infections
		for (int i = 0; i < tree.getNodeCount() - 1; i++) {
			int baseColour = colourAtBase[i];
			int parent = nodes[i].getParent().getNr();
			int parentColour = colourAtBase[parent];
			if (baseColour != parentColour && parentColour >= n) {
//System.err.println("#node " + (i+1));

				double tInf0 = segments.get(parentColour).birthTime;
				Node node = nodes[i];
				double tInf1 = node.getHeight() + node.getLength() * blockEndFraction.getArrayValue(node.getNr());
				double logP1 = logh_tr(tInf0, tInf1);
				logP += logP1;
			}
		}

//System.err.println("\n#contribution of cases in blocks");
		// contribution of cases in blocks
		for (int i = 0; i < tree.getNodeCount() - 1; i++) {
			if (blockCount.getValue(i) > 0) {
				double branchlength = nodes[i].getLength();
				double start = nodes[i].getHeight() + branchlength * blockStartFraction.getValue(i);
				double end   = nodes[i].getHeight() + branchlength * blockEndFraction.getValue(i);
				int blocks = blockCount.getValue(i);

				double logPBlock = getLogBlockLike(end - start, blocks, end - d);

//    			System.err.println("#node " + (i+1) + " " + logPBlock);
//    			System.err.println((tree.getRoot().getHeight() - end) + " - " + (tree.getRoot().getHeight() - start) + " = " + tau + " logPBlock=" + logPBlock);    			

				logP += logPBlock;
			} else  if (colourAtBase[i] != colourAtBase[nodes[i].getParent().getNr()]) {
				// blockCount[i] == 0 but parent colour differs from base colour
				// TODO: confirm there is no contribution ???
			}
		}

		return logP;
	}

	private double logS_tr(double t, double d) {
		return transmissionHazard.logS(tree.getRoot().getHeight() - d, tree.getRoot().getHeight() - t);
		// return transmissionHazard.logS(t, d);
	}

	private double logS_s(double t, double d) {
		return samplingHazard.logS(tree.getRoot().getHeight() - d, tree.getRoot().getHeight() - t);
		//return samplingHazard.logS(t, d);
	}

	private double logh_s(double t, double d) {
		return samplingHazard.logH(tree.getRoot().getHeight() - d, tree.getRoot().getHeight() - t);
		//return samplingHazard.logH(t, d);
	}

	private double logh_tr(double t, double d) {
		return transmissionHazard.logH(tree.getRoot().getHeight() - d, tree.getRoot().getHeight() - t);
		// return transmissionHazard.logH(t, d);
	}

	private List<SegmentIntervalList> segments;

	public double calculateCoalescent() {
		double logP = 0;
		for (SegmentIntervalList intervals : segments) {
			if (intervals != null) {
				if (conditionOnInfectionTime) {
					logP += calculateCoalescent(intervals, 0.0);
				} else {
					logP += calculateCoalescentUnconditioned(intervals, 0.0);
				}
			}
		}
		return logP;
	}

	public List<Double> calculateCoalescents() {
		segments = collectSegments();
		List<Double> logP = new ArrayList<>();
		for (SegmentIntervalList intervals : segments) {
			if (intervals != null) {
				if (conditionOnInfectionTime) {
					logP.add(calculateCoalescent(intervals, 0.0));
				} else {
					logP.add(calculateCoalescentUnconditioned(intervals, 0.0));
				}
			}
		}
		return logP;
	}

	class SegmentIntervalList implements IntervalList  {

		double birthTime;
		private List<Double> times = new ArrayList<>();
		private List<IntervalType> events = new ArrayList<>();

		private int intervalCount = 0;
		/**
		 * The widths of the intervals.
		 */
		private double[] intervals;
		/**
		 * The number of uncoalesced lineages within a particular interval.
		 */
		private int[] lineageCounts;


		@Override
		public int getIntervalCount() {
			return intervalCount;
		}

		@Override
		public int getSampleCount() {
			throw new RuntimeException("Not implemented yet");
		}

		@Override
		public double getInterval(int i) {
			if (i < 0 || i >= intervalCount) throw new IllegalArgumentException();
			return intervals[i];
		}

		@Override
		public int getLineageCount(int i) {
			if (i >= intervalCount) throw new IllegalArgumentException();
			return lineageCounts[i];
		}

		@Override
		public int getCoalescentEvents(int i) {
			if (i >= intervalCount) throw new IllegalArgumentException();
			if (i < intervalCount - 1) {
				return lineageCounts[i] - lineageCounts[i + 1];
			} else {
				return lineageCounts[i] - 1;
			}
		}

		@Override
		public IntervalType getIntervalType(int i) {
			if (i >= intervalCount) throw new IllegalArgumentException();
			int numEvents = getCoalescentEvents(i);

			if (numEvents > 0) return IntervalType.COALESCENT;
			else if (numEvents < 0) return IntervalType.SAMPLE;
			else return IntervalType.NOTHING;
		}

		@Override
		public double getTotalDuration() {
			return times.get(times.size() - 1) - times.get(0);
		}

		@Override
		public boolean isBinaryCoalescent() {
			throw new RuntimeException("Not implemented yet");
		}

		@Override
		public boolean isCoalescentOnly() {
			throw new RuntimeException("Not implemented yet");
		}


		public void calculateIntervals() {
			double multifurcationLimit = 0.0;
			int nodeCount = events.size();

			if (intervals == null || intervals.length != nodeCount) {
				intervals = new double[nodeCount];
				lineageCounts = new int[nodeCount];
			}

			// start is the time of the first tip
			double start = times.get(0);
			int numLines = 0;
			int nodeNo = 0;
			intervalCount = 0;
			while (nodeNo < nodeCount) {

				int lineagesRemoved = 0;
				int lineagesAdded = 0;

				double finish = times.get(nodeNo);
				double next;

				do {
					final int childIndex = nodeNo;
					final IntervalType type = events.get(childIndex);
					// don't use nodeNo from here on in do loop
					nodeNo += 1;
					if (type == IntervalType.SAMPLE) {
						lineagesAdded++;
					} else {
						lineagesRemoved++;

						// no mix of removed lineages when 0 th
						if (multifurcationLimit == 0.0) {
							break;
						}
					}

					if (nodeNo < nodeCount) {
						next = times.get(nodeNo);
					} else break;
				} while (Math.abs(next - finish) <= multifurcationLimit);

				if (lineagesAdded > 0) {

					if (intervalCount > 0 || ((finish - start) > multifurcationLimit)) {
						intervals[intervalCount] = finish - start;
						lineageCounts[intervalCount] = numLines;
						intervalCount += 1;
					}

					start = finish;
				}

				// add sample event
				numLines += lineagesAdded;

				if (lineagesRemoved > 0) {

					intervals[intervalCount] = finish - start;
					lineageCounts[intervalCount] = numLines;
					intervalCount += 1;
					start = finish;
				}
				// coalescent event
				numLines -= lineagesRemoved;
			}
		}

		public void addEvent(double time, IntervalType type) {
			if (times.size() == 0 || times.get(times.size()-1) < time) {
				times.add(time);
				events.add(type);
			} else {
				int index = Collections.binarySearch(times, time);
				if (index < 0) {
					index = -index - 1;
				}
				times.add(index, time);
				events.add(index, type);
			}
		}

		@Override
		public String toString() {
			if (times == null) {
				return "empty SegmentIntervalList";
			}
			String str = "";
			if (lineageCounts == null) {
				for (int i = 0;i < times.size(); i++) {
					//str += "(" + (events.get(i) == IntervalType.SAMPLE ? "S": "C") + " " + times.get(i) + ") ";
					str += "(" + (events.get(i) == IntervalType.SAMPLE ? "S": "C") + ") ";
				}
			} else {
				for (int i = 0;i < times.size(); i++) {
					//str += "(" + lineageCounts[i] + " " + (events.get(i) == IntervalType.SAMPLE ? "S": "C") + " " + times.get(i) + ") ";
					str += "(" + lineageCounts[i] + " " + (events.get(i) == IntervalType.SAMPLE ? "S": "C") + ") ";
				}
			}
			return str;
		}

	}

	private List<SegmentIntervalList> collectSegments() {
		List<SegmentIntervalList> segments = new ArrayList<>();
		int nodeCount = tree.getNodeCount();
		for (int i = 0; i < nodeCount; i++) {
			segments.add(null);
		}

		for (int i =  0; i < nodeCount; i++) {
			int colour = colourAtBase[i];
			if (segments.get(colour) == null) {
				segments.set(colour, new SegmentIntervalList());
			}
		}

		for (int i =  0; i < nodeCount; i++) {
			int colour = colourAtBase[i];
			Node node = tree.getNode(i);
			SegmentIntervalList intervals = (SegmentIntervalList) segments.get(colour);

			// add node event
			intervals.addEvent(node.getHeight(), node.isLeaf() ? IntervalType.SAMPLE : IntervalType.COALESCENT);
			if (!node.isRoot()) {
				int parentNr = node.getParent().getNr();
				int parentColour = colourAtBase[parentNr];
				if (colour != parentColour) {
					// add sampling event at top of block
					intervals = (SegmentIntervalList) segments.get(parentColour);
					double h = node.getHeight() + blockEndFraction.getValue(node.getNr()) * node.getLength();
					intervals.addEvent(h, IntervalType.SAMPLE);
					// set start of colour
					h = node.getHeight() + blockStartFraction.getValue(node.getNr()) * node.getLength();
					((SegmentIntervalList) segments.get(colour)).birthTime = h;
				}
			} else {
				((SegmentIntervalList) segments.get(colour)).birthTime = node.getHeight();
			}
		}

		// set origin in segment at root
		int colour = colourAtBase[tree.getNodeCount()-1];
		((SegmentIntervalList) segments.get(colour)).birthTime =
				origin != null ? origin.getArrayValue() : tree.getRoot().getHeight();

		for (IntervalList intervals : segments) {
			if (intervals != null) {
				((SegmentIntervalList)intervals).calculateIntervals();
			}
		}


		return segments;
	}

	/**
	 * Calculate contribution of coalescent NOT conditioned on infection time being before all coalescent events
	 * DOES NOT assume constant population size inside a host
	 */
	private double calculateCoalescentUnconditioned(SegmentIntervalList intervals, double threshold) {
		//private double calculateCoalescent(SegmentIntervalList intervals, double threshold) {
		double logL = 0.0;

		double startTime = 0.0;
		final int n = intervals.getIntervalCount();
		for (int i = 0; i < n; i++) {

			final double duration = intervals.getInterval(i);
			final double finishTime = startTime + duration;

			final double intervalArea = popSizeFunction.getIntegral(startTime, finishTime);
			if (intervalArea == 0 && duration > 1e-10) {
				/* the above test used to be duration != 0, but that leads to numerical issues on resume
				 * (https://github.com/CompEvol/beast2/issues/329) */
				return Double.NEGATIVE_INFINITY;
			}
			final int lineageCount = intervals.getLineageCount(i);

			final double kChoose2 = Binomial.choose2(lineageCount);
			// common part
			logL += -kChoose2 * intervalArea;

			if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {

				final double demographicAtCoalPoint = popSizeFunction.getPopSize(finishTime);

				// if value at end is many orders of magnitude different than mean over interval reject the interval
				// This is protection against cases where ridiculous infinitesimal
				// population size at the end of a linear interval drive coalescent values to infinity.

				if (duration == 0.0 || demographicAtCoalPoint * (intervalArea / duration) >= threshold) {
					//                if( duration == 0.0 || demographicAtCoalPoint >= threshold * (duration/intervalArea) ) {
					logL -= Math.log(demographicAtCoalPoint);
				} else {
					// remove this at some stage
					//  System.err.println("Warning: " + i + " " + demographicAtCoalPoint + " " + (intervalArea/duration) );
					return Double.NEGATIVE_INFINITY;
				}
			}
			startTime = finishTime;
		}

		return logL;
	}

	/**
	 * Calculate contribution of coalescent conditioned on infection time being before all coalescent events
	 * Assumes constant population size inside a host
	 */
	private double calculateCoalescent(SegmentIntervalList intervals, double threshold) {
		double logL = 0.0;

		double t0 = intervals.times.get(0);
		final int n = intervals.getIntervalCount();

		double tmax = intervals.birthTime;

		if(popSizeFunction instanceof LinearGrowth){

			t0 = t0 - tmax;
			double rate = ((LinearGrowth) popSizeFunction).getRate();

			for (int i = 0; i < n; i++) {

				final double duration = intervals.getInterval(i);
				final double t = t0 + duration;

				final int k = intervals.getLineageCount(i);

				if (t - 1e-6 > tmax) { // sanity check for debugging
					throw new RuntimeException("Programmer error: finish of interval > time of infection");
				}

				if (k > 1 && duration > 0) {
					logL -= (k*(k-1)/(2*rate))*(Math.log(-t0)-Math.log(-t));
					if(intervals.getIntervalType(i) == IntervalType.COALESCENT) {
						logL -= Math.log(rate);
						logL -= Math.log(-t);
					}

				}
				t0 = t;
			}

		} else {

			double N = popSizeFunction.getPopSize(0);
			for (int i = 0; i < n; i++) {

				final double duration = intervals.getInterval(i);
				final double t = t0 + duration;

				final int k = intervals.getLineageCount(i);

				if (t - 1e-6 > tmax) { // sanity check for debugging
					throw new RuntimeException("Programmer error: finish of interval > time of infection");
				}

				if (k > 1 && duration > 0) {
					final double rate = k * (k - 1) / (2.0 * N);

					switch (intervals.getIntervalType(i)) {
						case COALESCENT:
							// logL += log(rate*exp(-rate*duration))
							//       = -rate*duration + log(rate)
							logL += -rate * duration + Math.log(rate);
							// logL -= log(1-exp(-rate*(tmax-t0)))
							logL -= Math.log(1.0 - Math.exp(-rate * (tmax - t0)));
							break;
						case SAMPLE:
							// t0 = start of interval
							// tmax = time of infection
							// t = end of interval
							// logL += log(exp(-rate*(t-t0) - exp(-rate*(tmax-t0)))
							//       = log(exp(-rate*(t-t0) - exp(-rate*(tmax-t + t - t0)))
							//       = log(exp(-rate*(t-t0) - exp(-rate*(tmax-t)) * exp(-rate*(t-t0)))
							//       = log(exp(-rate*(t-t0) * (1 - exp(-rate*(tmax-t)))
							//       = -rate*(t-t0) + log(1 - exp(-rate*(tmax-t)))
							logL += -rate * duration + Math.log(1.0 - Math.exp(-rate * (tmax - t)));
							// logL -= log(1-exp(-rate*(tmax-t0)))
							logL -= Math.log(1.0 - Math.exp(-rate * (tmax - t0)));
							break;
						default:
							throw new RuntimeException("Don't know how to process " + intervals.getIntervalType(i));
					}
				}
				t0 = t;
			}
		}
		return logL;
	}



	public int getColour(int i) {
		if (updateColours) {
			calcColourAtBase();
		}
		return colourAtBase[i];
	}

	public int [] getColouring() {
		if (updateColours) {
			calcColourAtBase();
		}
		return colourAtBase;
	}

	public int [] getFreshColouring() {
		calcColourAtBase();
		return colourAtBase;
	}

	@Override
	public List<String> getConditions() {
		List<String> conditions = new ArrayList<>();
		conditions.add(blockStartFractionInput.get().getID());
		conditions.add(blockEndFractionInput.get().getID());
		conditions.add(blockCountInput.get().getID());
		//conditions.add(colourInput.get().getID());
		return conditions;
	}

	@Override
	public List<String> getArguments() {
		List<String> arguments = new ArrayList<>();
		arguments.add(treeInput.get().getID());
		return arguments;
	}

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub
	}



	@Override
	public void restore() {
		updateColours = true;
		super.restore();
	}

	@Override
	protected boolean requiresRecalculation() {
		updateColours = true;
		return true;
	}



	// recall newton's method: x_n+1 = x_n - f(x_n)/f'(x_n) to find a root of f 
	// here f = x- (1-C)*exp(lambda(x-1)) and f' is 1- (1-C)*lambda*exp(lambda(x-1))
	final static int  maxsteps = 1000;
	final static double tol=1e-6;
	private int n=0;

	private double f(double x, double Cs, double lambda) {
		return x - (1-Cs)*FastMath.exp(lambda*(x-1));
	}

//	private double fprime(double x, double Cs, double Ctr) {
//    	return 1 - (1-Cs)*Ctr*Math.exp(Ctr*(x-1));
//    }

	public double getp0(double Cs, double lambda, double x0) {
		n=0;
		double f = f(x0, Cs, lambda);
		while(Math.abs(f) > tol && n < maxsteps) {
			//x0 = x0 - f(x0, Cs, Ctr) / fprime(x0, Cs, Ctr); // Newton's method formula
			double tmp = (1-Cs)*FastMath.exp(lambda*(x0-1));
			x0 = x0 -(x0-tmp)/(1-tmp*lambda);
			f = f(x0, Cs, lambda);
			n=n+1;
		}
		if(n < maxsteps) {
			return x0;
		}
		throw new RuntimeException("The p0 algorithm did not converge after " + n + " iterations");
	}

	public double logGetIndivCondition(double p0, double t, double d) {
		final double TT = 1 - FastMath.exp(logS_tr(t, d)*(1-p0) + logS_s(t, d));
		if (TT == 0) {
			// something is wrong -- make sure the likelihood becomes NEGATIVE_INFINITY by returning POSITIVE_INFINITY
			return Double.POSITIVE_INFINITY;
		}
		final double logIndivCond = FastMath.log(TT);
//	    System.err.println("logGetIndivCondition(" +p0+"," + (t - d)+") = " + logIndivCond);
		return logIndivCond;
	}


	// recall p0 can be obtained with getp0(Cs,Ctr) 

	private double getPhi(double Cs,double lambda,double p0) {
		return(1 - p0*(1+ lambda*(1-p0)/(1-Cs)));
	}

	private double getRho(double phi) {
		return (1 - FastMath.exp(logS_tr(100, 0)*phi + logS_s(100, 0)));
	}

	final static double tol2=1e-7;
	final static int maxn = 1000000;
	private double getBlockCondition(double p0, double rho, double atr,double btr, double Yr) {
		double Z =0;// # init for the sum
		int n=1;
		double term = 1; // initialize the term at something > tol
		while ((term > tol2) & (n < maxn)) {
			term = FastMath.pow(1.0 - rho, n) * pgamma(Yr, n * atr, btr);
			Z = Z+term;
			n=n+1;
		}
//	    if (output==TRUE) {
//	        cat("Approximate sum:", Z, "\n")
//	        cat("Number of terms used:", n, "\n") 
//	        } 
		return Z;
	}

	private double getLogBlockLike(double tblock, int n, double Yr) {
		double blockLike = (1-FastMath.pow(rho,n)) * dgamma(tblock, n*atr, btr) / getBlockCondition(p0,rho, atr, btr, Yr);
//	    double blockLike = (1-FastMath.pow(rho,n)) * dgamma(tblock, n*a, b) / getBlockCondition(p0,rho, a, b, Yr);
		double logBlockLike = FastMath.log(blockLike);
//	    System.err.println("blockLike(" +tblock+"," + n +"," + Yr+") = " + blockLike);
		return logBlockLike;
	}

	// gives the density
	double dgamma(double x, double alpha, double rate) {
		if (x < 0) {
			throw new IllegalArgumentException("x should be non-negative");
		}
		return FastMath.pow(x * rate, alpha - 1) * rate * FastMath.exp(-x * rate) / FastMath.exp(Gamma.logGamma(alpha));
	}

	//  gives the cumulative distribution function
	double pgamma(double x, double shape, double rate) {
		if (x <= 0) {
			return 0;
		}
		try {
			return Gamma.regularizedGammaP(shape, x * rate);
		} catch (MathException e) {
			e.printStackTrace();
			return 0;
		}
	}

	public double calculateSampledHostContribution() {
		double d = endTime.getArrayValue();
		double logP = 0;
		int n = tree.getLeafNodeCount();
		Node [] nodes = tree.getNodesAsArray();
		calcColourAtBase();
		segments = collectSegments();

		if (origin != null) {
			if (origin.getArrayValue() < tree.getRoot().getHeight()) {
				return Double.NEGATIVE_INFINITY;
			}
		}

		// contribution of sampled cases
		for (int i = 0; i < n; i++) {
			double logP1 = 0;
			// contribution of not being sampled
			SegmentIntervalList intervals = segments.get(i);
			double start = intervals.birthTime;
			double end = intervals.times.get(0);
			logP1 += logh_s(start, end) + logS_s(start, end);
			// contribution of causing infections
			if (allowTransmissionsAfterSampling) {
				logP1 +=  logS_tr(start, d); // further contribution below
			} else {
				logP1 +=  logS_tr(start, end); // further contribution below
			}
			logP1 -= logGetIndivCondition(p0, start, d);
			if (Double.isInfinite(logP1) && logP1 > 0) {
				System.err.println("Numerical instability encountered: ");
				System.err.println(start + " " + d + " " + end + " " + p0);
				System.err.println(logS_tr(start, d));
				System.err.println(logS_tr(start, end));
				System.err.println(logGetIndivCondition(p0, start, d));
			}
			logP += logP1;
		}

		// further contribution of causing infections
		for (int i = 0; i < tree.getNodeCount() - 1; i++) {
			int baseColour = colourAtBase[i];
			int parent = nodes[i].getParent().getNr();
			int parentColour = colourAtBase[parent];
			if (baseColour != parentColour && parentColour < n) {
				double tInf0 = segments.get(parentColour).birthTime;
				Node node = nodes[i];
				double tInf1 = node.getHeight() + node.getLength() * blockEndFraction.getArrayValue(node.getNr());
				double logP1 = logh_tr(tInf0, tInf1);
				logP += logP1;
			}
		}

		return logP;
	}

	public double calculateUnsampledHostContribution() {
		double d = endTime.getArrayValue();
		double logP = 0;
		int n = tree.getLeafNodeCount();
		Node [] nodes = tree.getNodesAsArray();
		calcColourAtBase();
		segments = collectSegments();

		if (origin != null) {
			if (origin.getArrayValue() < tree.getRoot().getHeight()) {
				return Double.NEGATIVE_INFINITY;
			}
		}

		// contribution of unsampled cases
		for (int i = n; i < tree.getNodeCount(); i++) {
			if (colourAtBase[i] >= n) {
				// contribution of not being sampled
				SegmentIntervalList intervals = segments.get(i);
				if (intervals != null) {
					double start = intervals.birthTime;
					Double logP1 = logS_s(start, d);
					// contribution of causing infections
					logP += logS_tr(start, d); // further contribution below
					logP -= logGetIndivCondition(p0, start, d);

					logP += logP1;
				}
			}
		}

		// further contribution of causing infections
		for (int i = 0; i < tree.getNodeCount() - 1; i++) {
			int baseColour = colourAtBase[i];
			int parent = nodes[i].getParent().getNr();
			int parentColour = colourAtBase[parent];
			if (baseColour != parentColour && parentColour >= n) {

				double tInf0 = segments.get(parentColour).birthTime;
				Node node = nodes[i];
				double tInf1 = node.getHeight() + node.getLength() * blockEndFraction.getArrayValue(node.getNr());
				double logP1 = logh_tr(tInf0, tInf1);
				logP += logP1;
			}
		}
		return logP;
	}

	public double calculateBlockContribution() {
		double d = endTime.getArrayValue();
		double logP = 0;
		int n = tree.getLeafNodeCount();
		Node [] nodes = tree.getNodesAsArray();
		calcColourAtBase();
		segments = collectSegments();

		if (origin != null) {
			if (origin.getArrayValue() < tree.getRoot().getHeight()) {
				return Double.NEGATIVE_INFINITY;
			}
		}

		// contribution of cases in blocks
		for (int i = 0; i < tree.getNodeCount() - 1; i++) {
			if (blockCount.getValue(i) > 0) {
				double branchlength = nodes[i].getLength();
				double start = nodes[i].getHeight() + branchlength * blockStartFraction.getValue(i);
				double end   = nodes[i].getHeight() + branchlength * blockEndFraction.getValue(i);
				int blocks = blockCount.getValue(i);

				double logPBlock = getLogBlockLike(end - start, blocks, end - d);

				logP += logPBlock;
			} else  if (colourAtBase[i] != colourAtBase[nodes[i].getParent().getNr()]) {
				// blockCount[i] == 0 but parent colour differs from base colour
				// TODO: confirm there is no contribution ???
			}
		}

		return logP;
	}


//	public static void main(String[] args) {
//		TransmissionTreeLikelihood3 tl = new TransmissionTreeLikelihood3();
//		double Ctr=1.5;
//		double Cs=0.9; // Cs must now be strictly less than 1
//		double p0 = tl.getp0(Cs, Ctr, 0.1);
//		System.err.println("Root found: " + p0);
//		System.err.println("Function value at root: " + tl.f(p0, Cs, Ctr));
//		System.err.println("Number of iterations: " + tl.n);
//	}
}
