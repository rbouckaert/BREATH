package breath.test;


import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeParser;
import beast.base.evolution.tree.coalescent.ConstantPopulation;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import breath.distribution.GammaHazardFunction;
import breath.distribution.HazardFunction;
import breath.distribution.TransmissionTreeLikelihood1;

public class SimpleTransmissionLikelihoodTest {

	@Test
	public void testSimpleCase1() {
		TreeParser tree = new TreeParser("(Bob:1.7,Eve:2.1);");
		
        ConstantPopulation cp = new ConstantPopulation();
        cp.initByName("popSize", Double.toString(1.0));

        
        double bs1 = (1.7-0.5)/1.7;
        double bs2 = (2.1-0.65)/2.1;
        RealParameter blockStart = new RealParameter(); blockStart.initByName("dimension", 2, "value",  bs1+" "+bs2);
        RealParameter blockEnd = new RealParameter(); blockEnd.initByName("dimension", 2, "value",  bs1+" "+bs2);
        IntegerParameter blockcount = new IntegerParameter(); blockcount.initByName("dimension", 2, "value", "0 0");
        IntegerParameter colour = new IntegerParameter(); colour.initByName("dimension", 3, "value", "0 1 2");
        
        HazardFunction samplingHazard = new GammaHazardFunction();
        samplingHazard.initByName("C", "0.9", "shape", "2.5", "rate", "10.0");
        
        HazardFunction transmissionHazard = new GammaHazardFunction();
        transmissionHazard.initByName("C", "1.5", "shape", "2.0", "rate", "10.0");
        
        TransmissionTreeLikelihood1 coal = new TransmissionTreeLikelihood1();
        coal.initByName(
        		"tree", tree,
        		"populationModel", cp, 
        		"blockstart", blockStart, 
        		"blockend", blockEnd, 
        		"blockcount", blockcount, 
        		//"colour", colour,
        		"endTime", "0.0",
        		"samplingHazard", samplingHazard,
        		"transmissionHazard", transmissionHazard,
        		"origin", "2.5");
        
        coal.calcColourAtBase();
        double transmissionLikelihood = coal.calcTransmissionLikelihood();
        
        // assertEquals(-24.8279, transmissionLikelihood, 1e-5);
        assertEquals(-31.37552, transmissionLikelihood, 1e-5);
        
	}


	@Test
	public void testSimpleCase2() {
		TreeParser tree = new TreeParser("(Bob:1.7,Eve:2.1);");
		
        ConstantPopulation cp = new ConstantPopulation();
        cp.initByName("popSize", Double.toString(1.0));

        
        double bs1 = (1.7-0.5)/1.7;
        double bs2 = (2.1-0.65)/2.1;
        double bs3 = (2.1-1.8)/2.1;
        RealParameter blockStart = new RealParameter(); blockStart.initByName("dimension", 2, "value",  bs1+" "+bs3);
        blockStart.setID("blockStart");
        RealParameter blockEnd = new RealParameter(); blockEnd.initByName("dimension", 2, "value",  bs1+" "+bs2);
        blockEnd.setID("blockEnd");
        IntegerParameter blockcount = new IntegerParameter(); blockcount.initByName("dimension", 2, "value", "0 4");
        blockcount.setID("blockcount");
        IntegerParameter colour = new IntegerParameter(); colour.initByName("dimension", 3, "value", "0 1 2");
        colour.setID("colour");
        
        HazardFunction samplingHazard = new GammaHazardFunction();
        samplingHazard.initByName("C", "0.9", "shape", "2.5", "rate", "10.0");
        samplingHazard.setID("s");
        
        HazardFunction transmissionHazard = new GammaHazardFunction();
        transmissionHazard.initByName("C", "1.5", "shape", "2.0", "rate", "10.0");
        transmissionHazard.setID("tr");
        
        TransmissionTreeLikelihood1 coal = new TransmissionTreeLikelihood1();
        coal.initByName(
        		"tree", tree,
        		"populationModel", cp, 
        		"blockstart", blockStart, 
        		"blockend", blockEnd, 
        		"blockcount", blockcount, 
        		//"colour", colour,
        		"endTime", "0.0",
        		"samplingHazard", samplingHazard,
        		"transmissionHazard", transmissionHazard,
        		"origin", "3.0");
        
        coal.calcColourAtBase();
        double transmissionLikelihood = coal.calcTransmissionLikelihood();
        
//        assertEquals( -20.70768, transmissionLikelihood, 1e-5);
        assertEquals( -30.39658, transmissionLikelihood, 1e-5);
        
	}
	
	@Test
	public void testSimpleCase3() {
		TreeParser tree = new TreeParser("((t1:0.6587438122,t2:0.22):0.7863448577,(t3:0.3307722867,(t4:0.7084983373,t5:0.6330101104):0.6262222228):0.2);");
		
        ConstantPopulation cp = new ConstantPopulation();
        cp.initByName("popSize", Double.toString(1.0));

        Node [] nodes = tree.getNodesAsArray();
        
        double h = tree.getRoot().getHeight();
        
        double [] end = new double[8];
        double [] start = new double[8];
        int i = 0;
        // blockcount = -1
        start[i] = (h-1.3 - nodes[i].getHeight()) / nodes[i].getLength();
        end[i] = (h-0.94 - nodes[i].getHeight()) / nodes[i].getLength();
        i++;
        // blockcount = 0
        //start[i] = (h-0.92 - nodes[i].getHeight()) / nodes[i].getLength();
        start[i] = (h-0.85 - nodes[i].getHeight()) / nodes[i].getLength();
        end[i] = (h-0.85 - nodes[i].getHeight()) / nodes[i].getLength();
        i++;
        // blockcount = -1
        start[i] = (h-0.47 - nodes[i].getHeight()) / nodes[i].getLength();
        end[i] = (h-0.25 - nodes[i].getHeight()) / nodes[i].getLength();
        i++;
        // blockcount = -1
        //start[i] = (h-1.45 - nodes[i].getHeight()) / nodes[i].getLength();
        start[i] = (h-1.05 - nodes[i].getHeight()) / nodes[i].getLength();
        end[i] = (h-1.05 - nodes[i].getHeight()) / nodes[i].getLength();
        i++;
        // blockcount = 2
        start[i] = (h-1.33 - nodes[i].getHeight()) / nodes[i].getLength();
        end[i] = (h-0.96 - nodes[i].getHeight()) / nodes[i].getLength();
        i++;
        // blockcount = 4
        start[i] = (h-0.5 - nodes[i].getHeight()) / nodes[i].getLength();
        end[i] = (h-0.22 - nodes[i].getHeight()) / nodes[i].getLength();
        i++;
        // blockcount = 3
        start[i] = (h-0.72 - nodes[i].getHeight()) / nodes[i].getLength();
        end[i] = (h-0.35 - nodes[i].getHeight()) / nodes[i].getLength();
        i++;
        // blockcount = 0
        // Note: notes assume infection at `block end` (as defined in notes), not `block start` for test3
        start[i] = (h-0.15 - nodes[i].getHeight()) / nodes[i].getLength();
        //end[i] = (h-0.03 - nodes[i].getHeight()) / nodes[i].getLength();
        end[i] = (h-0.15 - nodes[i].getHeight()) / nodes[i].getLength();
        i++;
 
        
        RealParameter blockStart = new RealParameter(); blockStart.initByName("dimension", 8,       "value", start[0] + " " + start[1] + " " + start[2] + " " + start[3] + " " + start[4] + " " + start[5] + " " + start[6] + " " + start[7]);
        RealParameter blockEnd = new RealParameter(); blockEnd.initByName("dimension", 8,           "value", end[0] + " " + end[1] + " " + end[2] + " " + end[3] + " " + end[4] + " " + end[5] + " " + end[6] + " " + end[7]);
        IntegerParameter blockcount = new IntegerParameter(); blockcount.initByName("dimension", 8, "value", "-1 0 -1 -1 2 4 3 0");
        IntegerParameter colour = new IntegerParameter(); colour.initByName("dimension", 9,         "value", "0 1 2 3 4 0 3 2 8");
        
        HazardFunction samplingHazard = new GammaHazardFunction();
        samplingHazard.initByName("C", "0.9", "shape", "2.5", "rate", "10.0");
        samplingHazard.setID("s");
        
        HazardFunction transmissionHazard = new GammaHazardFunction();
        transmissionHazard.initByName("C", "1.5", "shape", "2.0", "rate", "10.0");
        transmissionHazard.setID("tr");
        
        TransmissionTreeLikelihood1 coal = new TransmissionTreeLikelihood1();
        coal.initByName(
        		"tree", tree,
        		"populationModel", cp, 
        		"blockstart", blockStart, 
        		"blockend", blockEnd, 
        		"blockcount", blockcount, 
        		//"colour", colour,
        		"endTime", (1.5347 - 1.7) + "",
        		"samplingHazard", samplingHazard,
        		"transmissionHazard", transmissionHazard,
        		"origin", "4.0");
        
        coal.calcColourAtBase();
        double transmissionLikelihood = coal.calcTransmissionLikelihood();
        
        // assertEquals( -25.45078, transmissionLikelihood, 1e-5);
        assertEquals( -53.740599, transmissionLikelihood, 1e-5);
	}	

	@Test
	public void testSimpleCase4() {
		TreeParser tree = new TreeParser("(t3:0.3307722867,(t4:0.7084983373,t5:0.6330101104):0.6262222228);");
		
        ConstantPopulation cp = new ConstantPopulation();
        cp.initByName("popSize", Double.toString(1.0));

        Node [] nodes = tree.getNodesAsArray();
        
        double h = tree.getRoot().getHeight();
        double e4 = (h-0.73 - nodes[2].getHeight()) / nodes[2].getLength(); 
        double e3 = (h-0.12 - nodes[0].getHeight()) / nodes[0].getLength(); 

        double s4 = (h-1.05 - nodes[2].getHeight()) / nodes[2].getLength(); 
        double s3 = e3; 
        
        
        RealParameter blockStart = new RealParameter(); blockStart.initByName("dimension", 4,       "value", s3 + " 0.5 " + s4 + " 0.5" );
        RealParameter blockEnd = new RealParameter(); blockEnd.initByName("dimension", 4,           "value", e3 + " 0.5 " + e4 + " 0.5");
        IntegerParameter blockcount = new IntegerParameter(); blockcount.initByName("dimension", 4, "value", "0 -1 2 -1");
        IntegerParameter colour = new IntegerParameter(); colour.initByName("dimension", 5,         "value", "0 1 2 1 1");
        
        HazardFunction samplingHazard = new GammaHazardFunction();
        samplingHazard.initByName("C", "0.9", "shape", "2.5", "rate", "10.0");
        samplingHazard.setID("s");
        
        HazardFunction transmissionHazard = new GammaHazardFunction();
        transmissionHazard.initByName("C", "1.5", "shape", "2.0", "rate", "10.0");
        transmissionHazard.setID("tr");
        
        TransmissionTreeLikelihood1 coal = new TransmissionTreeLikelihood1();
        coal.initByName(
        		"tree", tree,
        		"populationModel", cp, 
        		"blockstart", blockStart, 
        		"blockend", blockEnd, 
        		"blockcount", blockcount, 
        		//"colour", colour,
        		"endTime", (1.3347-1.5)+ "",
        		"samplingHazard", samplingHazard,
        		"transmissionHazard", transmissionHazard,
        		"origin", "4.0");
        
        coal.calcColourAtBase();
        double transmissionLikelihood = coal.calcTransmissionLikelihood();
        
        // assertEquals( -14.726, transmissionLikelihood, 1e-3);
        assertEquals( -84.83451, transmissionLikelihood, 1e-3);
        
	}
}
