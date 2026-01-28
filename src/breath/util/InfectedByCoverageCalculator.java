package breath.util;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;

import javax.imageio.ImageIO;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.core.Input.Validate;
import beast.base.inference.Runnable;
import beast.base.util.HeapSort;
import beastfx.app.beauti.ThemeProvider;
import beastfx.app.tools.Application;
import beastfx.app.tools.LogAnalyser;
import beastfx.app.util.Alert;
import beastfx.app.util.LogFile;
import beastfx.app.util.OutFile;
import javafx.application.Platform;
import javafx.embed.swing.JFXPanel;
import javafx.embed.swing.SwingFXUtils;
import javafx.scene.Node;
import javafx.scene.Scene;
import javafx.scene.SnapshotParameters;
import javafx.scene.chart.BarChart;
import javafx.scene.chart.CategoryAxis;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
import javafx.scene.control.Dialog;
import javafx.scene.control.DialogPane;
import javafx.scene.image.WritableImage;
import javafx.scene.layout.HBox;

@Description("Calculate coverage of who infected who")
public class InfectedByCoverageCalculator extends Runnable {
	final public Input<LogFile> truthInput = new Input<>("truth", "trace file with true infection information", Validate.REQUIRED);
	final public Input<Integer> skipLogLinesInput = new Input<>("skip", "numer of true log file lines to skip", 1);
	final public Input<File> logFilePrefixInput = new Input<>("logFilePrefix", "log file name without the number and '.log' missing. It is assumed there are as many log files as there are entries in the truth file", Validate.REQUIRED);
	final public Input<OutFile> outputInput = new Input<>("out", "output file, or stdout if not specified", new OutFile("[[none]]"));
	final public Input<Integer> burnInPercentageInput = new Input<>("burnin", "percentage of trees to used as burn-in (and will be ignored)", 10);	
	final public Input<Double> coverageInput = new Input<>("coverage", "percentage of coverage to be tested against (between 0 and 100, default 95)", 95.0);	
	final public Input<String> tagInput = new Input<>("tag", "name of the entry in log files containing infector of information. "
			+ "If not specified, assume all entries must be used");	
	final public Input<Boolean> includeUnsampledInput = new Input<>("includeUnsampled", "include unsampled infectors in true-vs-inferred plot", true);	
	final public Input<OutFile> pngFileInput = new Input<>("png", "name of file to write bar-chart plot", new OutFile("[[none]]"));	
	final public Input<Integer> binCountInput = new Input<>("bins", "number of bins=bars to use for the chart", 10);	

	
	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		LogAnalyser trueTrace = new LogAnalyser(truthInput.get().getPath(), 0, true, false);

		PrintStream out = System.out;
		if (outputInput.get() != null && !outputInput.get().getName().equals("[[none]]")) {
			Log.warning("Writing to file " + outputInput.get().getPath());
			out = new PrintStream(outputInput.get());
		}
		String tag = tagInput.get();
	
		int n = trueTrace.getLabels().size();
		if (tag != null) {
			n = 0;
			for (String str : trueTrace.getLabels()) {
				if (str.startsWith(tag + ".")) {
					n++;
				}
			}
		}
		
		out.print("Sample\t");
		for (int i = 0; i < n; i++) {
			out.print("Covered"+i+"\t");
		}
		out.println();
		
				
		int binCount = binCountInput.get();
		int [] truebins = new int [binCount];
		int [] totals = new int [truebins.length];
		double [] probsPerBin = new double [truebins.length];
		boolean includeUnsampled = includeUnsampledInput.get();
		
		int trueOffset = 0;
		if (tag != null) {
			trueOffset = trueTrace.indexof(tag+".1") - 1;
			if (trueOffset == -2) {
				trueOffset = 0;
				List<String> labels = trueTrace.getLabels();
				while (trueOffset < labels.size() && !labels.get(trueOffset).startsWith(tag)) {
					trueOffset++;
				}
			}
		}
		
		for (int i = 0; i < trueTrace.getTrace(0).length - skipLogLinesInput.get(); i++) {
			String filename = logFilePrefixInput.get().getPath() + i + ".log";
			if (new File(filename).exists()) {
				LogAnalyser trace = new LogAnalyser(filename, burnInPercentageInput.get(), true, false);
				if (i % 10 == 0) {
					Log.warning.print("|");
				} else {
					Log.warning.print(".");
				}
				
				out.print(i + "\t");
				int offset = 0;
				if (tag != null) {
					offset = trace.indexof(tag + ".1") - 1;
				}
				for (int j = 1; j <= n; j++) {
					// get true source value
					int trueSource = 0;
					trueSource = (n + 1 + (int)(double)trueTrace.getTrace(trueOffset+j)[i + skipLogLinesInput.get()]) % n;
					
					// collect info from trace
		        	double [] infectedBy = new double[n+1];
		        	Arrays.fill(infectedBy, 0);
		        	Double [] currenttrace = trace.getTrace(offset + j);
		        	for (Double d : currenttrace) {
		        		if (d < -1 || d >= n) {
		        			infectedBy[n]++;
		        		} else {
		        			infectedBy[(n+1+(int)(double)d) % n]++;
		        		}
		        	}
		        	
		        	// determine 95% coverage
		        	int[] index = new int[n+1];
		        	for (int k = 0; k <= n; k++) {
		        		index[k] = k;
		        	}
		        	HeapSort.sort(infectedBy, index);
	
		        	{
//System.out.print(i*10+j + "\t");
//int sum = 0;
//for (int k = 0; k <= n; k++) {
//	sum += infectedBy[index[k]];
//	System.out.print(infectedBy[index[k]] + "\t");
//}
//System.out.println(sum);
				}
	
		        	// is true value in the 95% coverage set?
		        	int threshold = (int)(currenttrace.length * coverageInput.get())/ 100;
		        	int sum = 0;
		        	int k = n;
		        	boolean covered = false;
		        	while (sum < threshold) {
		        		if (index[k] == trueSource) {
		        			covered = true;
		        			break;
		        		}
		        		sum += infectedBy[index[k]];
		        		k--;
		        	}
		        	
		        	// update true vs estimated bins
		        	for (int x = includeUnsampled ? 0 : 1; x < infectedBy.length; x++) {
		        		if (infectedBy[x] > 0) {
				        	int b = (int)(truebins.length*infectedBy[x]/currenttrace.length);
				        	if (b >= truebins.length) {
				        		b = truebins.length-1;
				        	}
				        	probsPerBin[b] += infectedBy[x]/currenttrace.length;
				        	totals[b]++;
		        		}
		        	}
		        	if (includeUnsampled || trueSource > 0) {
		        		if (infectedBy[trueSource] > 0) {
				        	int b = (int)(truebins.length*infectedBy[trueSource]/currenttrace.length);
				        	if (b >= truebins.length) {
				        		b = truebins.length-1;
				        	}
				        	truebins[b]++;
		        		}
		        	}
					out.print((covered?1:0) + "\t");
				}
				out.println();
			}
		}
		
		System.out.println();
		System.out.print("totals: ");
		System.out.println(Arrays.toString(totals));
		System.out.print("true: ");
		System.out.println(Arrays.toString(truebins));
		System.out.print("Observed_probs: ");
    	for (int x = 0; x < truebins.length; x++) {
    		System.out.print((double)truebins[x]/totals[x]);
    		if (x < truebins.length-1) {
    			System.out.print(", ");
    		}
    	}    	
    	for (int x = 0; x < truebins.length; x++) {
    		probsPerBin[x] /= totals[x];
    	}
		System.out.print("\nExpected_probs: ");
    	for (int x = 0; x < truebins.length; x++) {
    		System.out.print(probsPerBin[x]);
    		if (x < truebins.length-1) {
    			System.out.print(", ");
    		}
    	}
		
		if (outputInput.get() != null && !outputInput.get().getName().equals("[[none]]")) {
			out.close();
		}

		if (pngFileInput.get() != null && !pngFileInput.get().getName().equals("[[none]]")) {
			showCoveragePlot(truebins, totals, pngFileInput.get(), includeUnsampledInput.get());
		}
		
		Log.warning("\nDone");
		Platform.exit();
	}
	
	
    static public void showCoveragePlot(int [] truebins, int [] totals, File pngfile, boolean includeUnsampleds) {
		// this initialised the javafx toolkit
		new JFXPanel();
		Platform.runLater(() -> {
	
	    	HBox root = new HBox();
	
	        Scene scene = new Scene(root, 480, 330);
	        CategoryAxis xAxis = new CategoryAxis();
	        xAxis.setLabel("Inferred");
	        
	        NumberAxis yAxis = new NumberAxis();
	        yAxis.setLabel("Actual");
	        yAxis.setAutoRanging(false);
	        yAxis.setLowerBound(0);
	        yAxis.setUpperBound(100);
	
	        BarChart barChart = new BarChart(xAxis, yAxis);
	        barChart.setTitle("Infectors true vs inferred " + 
	        		(includeUnsampleds? "with" : "without") + " unsampled hosts");
	
	        XYChart.Series data = new XYChart.Series<String, Number>();
	
	        for (int i = 0; i < truebins.length; i++) {
	        	data.getData().add(new XYChart.Data<>(
	        			(i * 100)/truebins.length + "-" + ((i+1) * 100)/truebins.length + "\n" + totals[i], 
	        			totals[i] > 0 ? 100.0 * truebins[i]/totals[i] : 0));
	        }
	
	        barChart.getData().add(data);
	        barChart.setLegendVisible(false);
	        
	
	        root.getChildren().add(barChart);
	
			Dialog<Node> alert = new javafx.scene.control.Dialog<>();
			DialogPane pane = new DialogPane();
			pane.setContent(root);
			alert.setDialogPane(pane);
			alert.setHeaderText("coverage");
			alert.getDialogPane().getButtonTypes().addAll(Alert.CLOSED_OPTION);
			pane.setPrefHeight(600);
			pane.setPrefWidth(600);
			alert.setResizable(true);
			ThemeProvider.loadStyleSheet(alert.getDialogPane().getScene());
			

			SnapshotParameters param = new SnapshotParameters();
		    param.setDepthBuffer(true);
		    WritableImage snapshot = root.snapshot(param, null);
		    BufferedImage tempImg = SwingFXUtils.fromFXImage(snapshot, null);
		    try {
			  Graphics g = tempImg.getGraphics();
			  g.setColor(Color.black);
			  g.drawLine(77, 429, 490, 52);
		      ImageIO.write(tempImg, "png", new FileOutputStream(pngfile));
		    } catch (IOException e) {
		      e.printStackTrace();
		    }			
		    System.exit(0);
			// alert.showAndWait();
	    });

    }
	
	
	
	
	public static void main(String[] args) throws Exception {
		new Application(new InfectedByCoverageCalculator(), "InfectedByCoverageCalculator", args);
	}

}
