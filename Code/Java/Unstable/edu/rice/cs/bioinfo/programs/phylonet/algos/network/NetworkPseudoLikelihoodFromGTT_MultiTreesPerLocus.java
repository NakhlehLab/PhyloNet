/*
 * Copyright (c) 2013 Rice University.
 *
 * This file is part of PhyloNet.
 *
 * PhyloNet is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PhyloNet is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PhyloNet.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.rice.cs.bioinfo.programs.phylonet.algos.network;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 2/11/13
 * Time: 11:40 AM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkPseudoLikelihoodFromGTT_MultiTreesPerLocus extends NetworkPseudoLikelihoodFromGTT {

    //for simulation only
/*
    public void summarizeData(List originalGTs, Map<String,String> allele2species, List allTriplets, List tripletFrequencies){
        computeTripleFrequenciesInGTsForSimulation(originalGTs, allTriplets, tripletFrequencies);
    }
*/

    protected void summarizeData(List originalGTs, Map<String,String> allele2species, List allTriplets, List tripletFrequencies){
        int treeID = 0;
        Map<String, Integer> distinctTree2ID = new HashMap<>();
        List<MutableTuple<Tree,Double>> distinctGTs = new ArrayList<>();
        List<List<MutableTuple<Integer,Double>>> treeCorrespondences = new ArrayList<>();
        for(Object o: originalGTs) {
            List<MutableTuple<Tree,Double>> treesForOneLocus = (List<MutableTuple<Tree,Double>>)o;
            Map<String, Integer> tree2infoIndex = new HashMap<String, Integer>();
            List<MutableTuple<Integer,Double>> infoList = new ArrayList<MutableTuple<Integer, Double>>();
            for (MutableTuple<Tree, Double> gtTuple : treesForOneLocus) {
                for (TNode node : gtTuple.Item1.getNodes()) {
                    node.setParentDistance(TNode.NO_DISTANCE);
                }
                String exp = Trees.getLexicographicNewickString(gtTuple.Item1, allele2species);
                Integer existingTreeID = distinctTree2ID.get(exp);
                if (existingTreeID == null) {
                    existingTreeID = treeID;
                    distinctGTs.add(new MutableTuple<Tree, Double>(gtTuple.Item1,gtTuple.Item2));
                    distinctTree2ID.put(exp, existingTreeID);
                    tree2infoIndex.put(exp, infoList.size());
                    infoList.add(new MutableTuple(treeID, gtTuple.Item2));
                    treeID++;
                } else {
                    distinctGTs.get(existingTreeID).Item2 += gtTuple.Item2;
                    Integer infoID = tree2infoIndex.get(exp);
                    if(infoID == null){
                        tree2infoIndex.put(exp, infoList.size());
                        infoList.add(new MutableTuple(existingTreeID, gtTuple.Item2));
                    }
                    else{
                        infoList.get(infoID).Item2 += gtTuple.Item2;
                    }
                }
            }
            treeCorrespondences.add(infoList);
        }
        computeTripleFrequenciesInGTs_original(distinctGTs, allele2species, treeCorrespondences, allTriplets, tripletFrequencies);
    }


    protected void summarizeData(List originalGTs, Map<String,String> allele2species, List distinctGTs, List allTriplets, List tripletFrequencies){
        int treeID = 0;
        Map<String, Integer> distinctTree2ID = new HashMap<>();
        List<List<MutableTuple<Integer,Double>>> treeCorrespondences = new ArrayList<>();
        int i=0;
        for(Object o: originalGTs) {
            List<MutableTuple<Tree,Double>> treesForOneLocus = (List<MutableTuple<Tree,Double>>)o;
            Map<String, Integer> tree2infoIndex = new HashMap<String, Integer>();
            List<MutableTuple<Integer,Double>> infoList = new ArrayList<MutableTuple<Integer, Double>>();
            for (MutableTuple<Tree, Double> gtTuple : treesForOneLocus) {
                for (TNode node : gtTuple.Item1.getNodes()) {
                    node.setParentDistance(TNode.NO_DISTANCE);
                }
                String exp = Trees.getLexicographicNewickString(gtTuple.Item1, allele2species);
                Integer existingTreeID = distinctTree2ID.get(exp);
                if (existingTreeID == null) {
                    existingTreeID = treeID;
                    distinctGTs.add(new MutableTuple<Tree, Double>(gtTuple.Item1,gtTuple.Item2));
                    distinctTree2ID.put(exp, existingTreeID);
                    tree2infoIndex.put(exp, infoList.size());
                    infoList.add(new MutableTuple(treeID, gtTuple.Item2));
                    treeID++;
                } else {
                    ((MutableTuple<Tree,Double>)(distinctGTs.get(existingTreeID))).Item2 += gtTuple.Item2;
                    Integer infoID = tree2infoIndex.get(exp);
                    if(infoID == null){
                        tree2infoIndex.put(exp, infoList.size());
                        infoList.add(new MutableTuple(existingTreeID, gtTuple.Item2));
                    }
                    else{
                        infoList.get(infoID).Item2 += gtTuple.Item2;
                    }
                }
            }
            treeCorrespondences.add(infoList);
        }
        computeTripleFrequenciesInGTs_original(distinctGTs, allele2species, treeCorrespondences, allTriplets, tripletFrequencies);
    }



/*
    protected void summarizeDataForSimulation(List originalGTs, Map<String,String> allele2species, List distinctGTs, List allTriplets, List tripletFrequencies){
        //Set<String> distinctTree2ID = new HashSet<>();
        for(Object o: originalGTs) {
            List<MutableTuple<Tree,Double>> treesForOneLocus = (List<MutableTuple<Tree,Double>>)o;
            for (MutableTuple<Tree, Double> gtTuple : treesForOneLocus) {
                distinctGTs.add(new MutableTuple<Tree, Double>(gtTuple.Item1,gtTuple.Item2));
            }
        }
        computeTripleFrequenciesInGTsForSimulation(originalGTs, allTriplets, tripletFrequencies);
    }
*/


    private void computeTripleFrequenciesInGTs_original(List<MutableTuple<Tree,Double>> gts, Map<String,String> allele2species, List<List<MutableTuple<Integer,Double>>> treeCorrespondences, List<String> allTriplets, List<List<Tuple<double[][],Double>>> tripletFrequencies){
        List<Map<String, double[]>> triple2countList = new ArrayList<>();
        Set<String> allAlleleSet = new HashSet<>();

        for(MutableTuple<Tree,Double> gt: gts){
            Map<String, double[]> result;

            if(allele2species==null) {
                result = computeTripleFrequenciesFromSingleGT(gt.Item1);

            }else{
                result = computeTripleFrequenciesFromSingleGT(gt.Item1, allele2species);
            }
            allAlleleSet.addAll(result.keySet());
            triple2countList.add(result);
        }

        allTriplets.addAll(allAlleleSet);
        for(List<MutableTuple<Integer,Double>> oneLocus: treeCorrespondences){
            List<Tuple<double[][],Double>> freqList = new ArrayList<>();
            tripletFrequencies.add(freqList);
            for(MutableTuple<Integer,Double> gtInfo: oneLocus){
                double[][] gtTripleFreq = new double[allTriplets.size()][3];
                freqList.add(new Tuple<double[][], Double>(gtTripleFreq, gtInfo.Item2));
                Map<String, double[]> triple2count = triple2countList.get(gtInfo.Item1);
                int index = 0;
                for(String triplet: allTriplets){
                    double[] baseFreq = triple2count.get(triplet);
                    for(int i=0; i<3; i++){
                        gtTripleFreq[index][i] = baseFreq[i];
                    }
                    index++;
                }
            }
        }

    }

/*
    protected void computeTripleFrequenciesInGTsForSimulation(List originalGTs, List allTriplets, List tripletFrequencies){
        //Set<String> distinctTree2ID = new HashSet<>();
        Map<String, double[]> triplet2freq = new HashMap<>();
        //int index = 0;
        for(Object o: originalGTs) {
            //System.out.println(index++);
            List<MutableTuple<Tree,Double>> treesForOneLocus = (List<MutableTuple<Tree,Double>>)o;
            for (MutableTuple<Tree, Double> gtTuple : treesForOneLocus) {
                Map<String, double[]> result = computeTripleFrequenciesFromSingleGT(gtTuple.Item1);
                for(Map.Entry<String,double[]> entry: result.entrySet()) {
                    double[] freq = triplet2freq.get(entry.getKey());
                    if (freq == null) {
                        freq = new double[3];
                        triplet2freq.put(entry.getKey(), freq);
                    }
                    for (int k = 0; k < freq.length; k++) {
                        freq[k] += entry.getValue()[k] * gtTuple.Item2;
                    }

                }
            }

        }
        for(Map.Entry<String, double[]> entry: triplet2freq.entrySet()){
            allTriplets.add(entry.getKey());
            tripletFrequencies.add(entry.getValue());
        }
    }
*/

    private void computeTripleFrequenciesInGTs(List<MutableTuple<Tree,Double>> gts, Map<String,String> allele2species, List<List<MutableTuple<Integer,Double>>> treeCorrespondences, List<String> allTriplets, List<double[]> tripletFrequencies){
        List<Map<String, double[]>> triple2countList = new ArrayList<>();
        Set<String> allAlleleSet = new HashSet<>();

        for(MutableTuple<Tree,Double> gt: gts){
            //System.out.println(k++ + "/" + size);
            Map<String, double[]> result;

            if(allele2species==null) {
                result = computeTripleFrequenciesFromSingleGT(gt.Item1);

            }else{
                result = computeTripleFrequenciesFromSingleGT(gt.Item1, allele2species);
            }
            allAlleleSet.addAll(result.keySet());
            triple2countList.add(result);
        }

        allTriplets.addAll(allAlleleSet);
        for(int i=0; i<allTriplets.size(); i++){
            tripletFrequencies.add(new double[3]);
        }

        for(List<MutableTuple<Integer,Double>> oneLocus: treeCorrespondences){
            double[][] freqList = new double[allAlleleSet.size()][3];
            double totalWeight = 0;
            for(MutableTuple<Integer,Double> gtInfo: oneLocus){
                totalWeight += gtInfo.Item2;
                Map<String, double[]> triple2count = triple2countList.get(gtInfo.Item1);
                int index = 0;
                for(String triplet: allTriplets){
                    double[] baseFreq = triple2count.get(triplet);
                    if(baseFreq == null)continue;
                    for(int i=0; i<3; i++){
                        freqList[index][i] += baseFreq[i];
                    }
                    index++;
                }
            }
            int index = 0;
            for(double[] freq: tripletFrequencies){
                for(int i=0; i<3; i++){
                    freq[i] += freqList[index][i]/totalWeight;
                }
                index++;
            }
            /*
            for(Tuple<double[][],Double> tuple: freqList){
                System.out.println(tuple.Item2);
                for(double[] freq: tuple.Item1){
                    System.out.println(Arrays.toString(freq));
                }
            }
            System.out.println();
            */
        }

/*
        int index = 0;
        for(Object data: allTriplets){
            System.out.println(index++ + ": " + data);
        }
*/
    }

    /*
    protected double calculateFinalLikelihood(double[][] probs, List tripletFrequencies){
        int index = 0;
        double totalProb = 0;
        for(Object o: tripletFrequencies){
            double[] freqs = (double[])o;
            for(int i=0; i<3; i++){
                totalProb += freqs[i] * Math.log(probs[index][i]);
            }
            index++;
        }
        return totalProb;
    }
    */


    protected double calculateFinalLikelihood(double[][] probs, List tripletFrequencies){

        double totalProb = 0;
        for(Object o: tripletFrequencies){
            List<Tuple<double[][],Double>> freqForOneLocus = (List<Tuple<double[][],Double>>)o;
            double probForOneLocus = 0;
            double totalWeight = 0;
            for(Tuple<double[][],Double> freqs: freqForOneLocus){
                double probForOneTree = 1;
                totalWeight += freqs.Item2;
                for(int i=0; i<probs.length; i++){
                    for(int j=0; j<3; j++){
                        if(probs[i][j] == 0)return Double.NEGATIVE_INFINITY;
                        probForOneTree *= Math.pow(probs[i][j],freqs.Item1[i][j]);
                    }
                }
                //System.out.println(probForOneTree);
                probForOneLocus += probForOneTree * freqs.Item2;
            }
            probForOneLocus /= totalWeight;
            totalProb += Math.log(probForOneLocus);

        }
        return totalProb;
    }


/*
    protected void summarizeDataForBird1(String gtFile, List allTriplets, List tripletFrequencies){

        try {
            BufferedReader br = new BufferedReader(new FileReader(gtFile));
            String line = br.readLine();
            while (line!= null) {
                List<String> freqForOneLoci = null;
                //boolean canBeRooted = false;

                for(int i=0; i<200; i++) {
                    if(i==0 || freqForOneLoci!=null) {
                        NewickReader nr = new NewickReader(new StringReader(line));
                        Tree tr = nr.readTree();
                        TNode outgroup = tr.getNode("STRCA");
                        if (outgroup == null) {
                            outgroup = tr.getNode("TINMA");
                        }
                        if (outgroup == null) {
                            if(i!=0) {
                                throw new RuntimeException("Something's wrong");
                            }
                        }
                        else {
                            if(i==0){

                                freqForOneLoci = new ArrayList<>();
                                tripletFrequencies.add(freqForOneLoci);
                            }
                            tr.rerootTreeAtEdge(outgroup.getName());
                            Map<String, double[]> triplet2freq = computeTripleFrequenciesFromSingleGT(tr);
                            if (allTriplets.isEmpty()) {
                                allTriplets.addAll(triplet2freq.keySet());
                                if (allTriplets.size() != 17296) {
                                    throw new RuntimeException("Not all 48 species are included");
                                }
                            }
                            String gtTripleFreq = "";

                            //int index = 0;
                            for (Object triplet : allTriplets) {
                                double[] baseFreq = triplet2freq.get(triplet);
                                if (baseFreq == null){
                                    gtTripleFreq += '-';
                                    continue;
                                }
                                boolean found = false;
                                for (int k = 0; k < 3; k++) {
                                    if (baseFreq[k] == 1) {
                                        gtTripleFreq += k;
                                        found = true;
                                        break;
                                    }
                                }
                                if (!found) {
                                    System.out.println("Wrong");
                                }
                                //index++;
                            }
                            if(gtTripleFreq.length()!=17296){
                                throw new RuntimeException("frequency wrong");
                            }
                            freqForOneLoci.add(gtTripleFreq);
                        }
                    }
                    line = br.readLine();
                }
                if(tripletFrequencies.size()==5){
                    break;
                }
            }
            br.close();


        }catch (Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }

    }





    //summarize triplets frequency inside of bootstrap
    protected void summarizeDataForBird2(String gtFile, List allTriplets, List tripletFrequencies){

        try {
            BufferedReader br = new BufferedReader(new FileReader(gtFile));
            String line = br.readLine();
            Map<String, double[]> triplet2freq = new HashMap<>();
            int index = 0;
            while (line!= null) {
                boolean canBeRooted = false;
                for(int i=0; i<200; i++) {
                    if(i==0 || canBeRooted) {
                        NewickReader nr = new NewickReader(new StringReader(line));
                        Tree tr = nr.readTree();
                        TNode outgroup = tr.getNode("STRCA");
                        if (outgroup == null) {
                            outgroup = tr.getNode("TINMA");
                        }
                        if (outgroup == null) {
                            if(i!=0) {
                                throw new RuntimeException("Something's wrong");
                            }
                        }
                        else {
                            if(i==0){
                                index++;
                                canBeRooted = true;
                            }
                            //System.out.println(index + ": " + i);
                            tr.rerootTreeAtEdge(outgroup.getName());
                            Map<String, double[]> result = computeTripleFrequenciesFromSingleGT(tr);
                            for(Map.Entry<String,double[]> entry: result.entrySet()) {
                                double[] freq = triplet2freq.get(entry.getKey());
                                if (freq == null) {
                                    freq = new double[3];
                                    triplet2freq.put(entry.getKey(), freq);
                                }
                                for (int k = 0; k < freq.length; k++) {
                                    freq[k] += entry.getValue()[k]/200;
                                }

                            }
                        }
                    }
                    line = br.readLine();
                }
            }
            br.close();
            for(Map.Entry<String, double[]> entry: triplet2freq.entrySet()){
                allTriplets.add(entry.getKey());
                tripletFrequencies.add(entry.getValue());
            }
        }catch (Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }

    }


    //reading triplets
    protected void summarizeDataForBird3(String gtFile, List allTriplets, List tripletFrequencies){

        try {
            BufferedReader br = new BufferedReader(new FileReader(gtFile));
            String line;
            while ((line=br.readLine())!= null) {
                String[] field1 = line.split(":");
                allTriplets.add(field1[0]);
                field1[1] = field1[1].substring(1,field1[1].length()-1);
                String[] field2 = field1[1].split(",");
                double[] freq = new double[3];
                for(int k=0; k<3; k++){
                    freq[k] = Double.parseDouble(field2[k].trim());
                }
                tripletFrequencies.add(freq);
            }
            br.close();
        }catch (Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }

    }

    //summarize triplets frequency for best ML gts
    protected void summarizeDataForBird4(String gtFile, List allTriplets, List tripletFrequencies){

        try {
            BufferedReader br = new BufferedReader(new FileReader(gtFile));
            Map<String, double[]> triplet2freq = new HashMap<>();
            String line;
            //int index = 0;
            while ((line = br.readLine() ) != null) {

                NewickReader nr = new NewickReader(new StringReader(line));
                STITree tr = new STITree();
                nr.readTree(tr);
                tr.removeNode("STRCA");
                tr.removeNode("TINMA");
                while(tr.getRoot().getChildCount()==1){
                    TNode node = (TNode)tr.getRoot().getChildren().iterator().next();
                    tr.getRoot().removeChild((TMutableNode)node, true);
                }
                Trees.handleBootStrapInTree(tr, 70);
                Map<String, double[]> result = computeTripleFrequenciesFromSingleGT(tr);
                for (Map.Entry<String, double[]> entry : result.entrySet()) {
                    double[] freq = triplet2freq.get(entry.getKey());
                    if (freq == null) {
                        freq = new double[3];
                        triplet2freq.put(entry.getKey(), freq);
                    }
                    //System.out.println(entry.getKey()+": "+Arrays.toString(entry.getValue()));
                    for (int k = 0; k < freq.length; k++) {
                        freq[k] += entry.getValue()[k];
                    }
                }
                //System.out.println();
            }

            br.close();
            for(Map.Entry<String, double[]> entry: triplet2freq.entrySet()){
                allTriplets.add(entry.getKey());
                tripletFrequencies.add(entry.getValue());
            }
        }catch (Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }

    }


    public double computeProbabilityForBird(Network<Object> speciesNetwork, List allTriplets, List tripleFrequencies) {
        GeneTreeProbabilityPseudo calculator = new GeneTreeProbabilityPseudo();
        calculator.initialize(speciesNetwork);
        double[][] probs = new double[allTriplets.size()][3];
        Thread[] myThreads = new Thread[_numThreads];
        //System.out.println(speciesNetwork);
        if(_numThreads>1) {
            calculator.setParallel(true);
            for (int i = 0; i < _numThreads; i++) {
                myThreads[i] = new MyThread(speciesNetwork, calculator, allTriplets, probs);
                myThreads[i].start();
            }
            for (int i = 0; i < _numThreads; i++) {
                try {
                    myThreads[i].join();
                } catch (InterruptedException ignore) {
                }
            }
        }else{
            try {
                calculator.computePseudoLikelihood(speciesNetwork, allTriplets, probs);
            }catch (Exception e){
                System.out.println(speciesNetwork);
                System.err.println(e.getMessage());
                e.getStackTrace();
                System.exit(-1);

            }
        }
        double totalProb = calculateFinalLikelihoodForBird(probs, tripleFrequencies);
        //System.out.println(speciesNetwork);
        //System.out.println(totalProb);
        return totalProb;
    }


    protected double calculateFinalLikelihoodForBird(double[][] probs, List tripletFrequencies){

        BigDecimal totalProb = new BigDecimal(0);
        for(Object o: tripletFrequencies){
            List<String> freqForOneLocus = (List<String>)o;
            double probForOneLocus = 0;
            int totalWeight = 0;
            for(String freqArrays: freqForOneLocus){
                char[] freqs = freqArrays.toCharArray();
                double probForOneTree = 1;
                totalWeight ++;
                for(int i=0; i<probs.length; i++){
                    int freq = freqs[i] - 48;
                    if(freq<0 || freq>2)continue;
                    probForOneTree *= probs[i][freq];
                    System.out.println(probForOneTree);
                }
                //System.out.println(probForOneTree);
                probForOneLocus += probForOneTree;
            }
            probForOneLocus /= totalWeight;
            totalProb += Math.log(probForOneLocus);

        }
        return totalProb;
    }


    protected double calculateFinalLikelihoodForBird(double[][] probs, List tripletFrequencies){
        int index = 0;
        double totalProb = 0;
        for(Object o: tripletFrequencies){
            double[] freqs = (double[])o;
            for(int i=0; i<3; i++){
                totalProb += freqs[i] * Math.log(probs[index][i]);
            }
            index++;
        }
        return totalProb;
    }

*/


}
