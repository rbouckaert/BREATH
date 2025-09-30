package breath.util;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Runnable;
import beastfx.app.tools.Application;
import beastfx.app.util.OutFile;

@Description("Create coverage plot from InfectedByCoverageCalculator screen output")
public class CoveragePlot extends Runnable {
	final public Input<OutFile> pngFileInput = new Input<>("png", "name of file to write bar-chart plot", new OutFile("[[none]]"));	
	final public Input<Boolean> includeUnsampledInput = new Input<>("includeUnsampled", "include unsampled infectors in true-vs-inferred plot", true);	
	final public Input<List<String>> truebinsdInput = new Input<>("true:", "java array representation of true bins, e.g. `[102,30,12]`", new ArrayList<>());	
	final public Input<List<String>> totalbinsdInput = new Input<>("totals:", "java array representation of total bins", new ArrayList<>());	


	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
	
		int [] truebins = getArray(truebinsdInput.get()); 
		int [] totals = getArray(totalbinsdInput.get()); 
		
		File pngfile = pngFileInput.get(); 
		boolean includeUnsampleds = includeUnsampledInput.get();
		InfectedByCoverageCalculator.showCoveragePlot(truebins, totals, pngfile, includeUnsampleds);
	}
	

	private int[] getArray(List<String> strs) {
		//String [] strs = str.split(",");
		int [] array = new int[strs.size()];
		for (int i = 0; i < strs.size(); i++) {
			String str = strs.get(i);
			str = str.replaceAll("\\[", "");
			str = str.replaceAll("\\]", "");
			str = str.replaceAll(",", "");
			array[i] = Integer.parseInt(str.trim());
		}
		return array;
	}

	public static void main(String[] args) throws Exception {
		new Application(new CoveragePlot(), "Coverage plot", args);
	}
}
