package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.summary;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * Summarize gene tree samples and network samples
 * Created by wendingqiao on 4/12/16.
 */
public class SampleSummary {

    private BufferedWriter _outFile;
    private Utils.SampleType _type;
    private List<String> _samples;

    public SampleSummary(String filename, Utils.SampleType sampleType) {
        try {
            // overwrite file if exists previously
            _outFile = new BufferedWriter(new FileWriter(new File(Utils._OUT_DIRECTORY, filename + ".log"), false));
            _type = sampleType;
            _samples = new ArrayList<>();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void addSample(String sample) {
        try {
            _outFile.write(sample);
            _outFile.write("\n");
            _samples.add(sample);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void summary() {
        String str = (_type == Utils.SampleType.Tree) ? summarizeTreeTopos() :
                (_type == Utils.SampleType.Network) ? summarizeNetTopos() :
                        (_type == Utils.SampleType.DoubleParam) ? summarizeDoubles() : "";
        try {
            _outFile.write("----------------------- Summary -----------------------\n");
            _outFile.write(str);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void close() {
        try {
            _outFile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private String summarizeDoubles() {
        double avg = 0;
        double std = 0;
        double min = Double.MAX_VALUE;
        double max = 0;
        for(String s : _samples) {
            double d = Double.parseDouble(s);
            avg += d;
            std += d * d;
            min = Math.min(min, d);
            max = Math.max(max, d);
        }
        avg /= _samples.size();
        std = Math.sqrt(std / _samples.size() - avg * avg);
        return String.format("avg = %f, std = %f, min = %f, max = %f", avg, std, min, max);
    }

    private String summarizeTreeTopos() {
        StringBuilder sb = new StringBuilder();
        List<MutableTuple<Tree, Integer>> list = new ArrayList<>();
        try {
            for(String s : _samples) {
                NewickReader nr = new NewickReader(new StringReader(s));
                Tree temp = nr.readTree();
                boolean found = false;
                for(MutableTuple<Tree, Integer> tup : list) {
                    if(Trees.haveSameRootedTopology(temp, tup.Item1)) {
                        found = true;
                        tup.Item2 = tup.Item2 + 1;
                    }
                }
                if(!found) {
                    list.add(new MutableTuple<Tree, Integer>(temp, 1));
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        Collections.sort(list, new Comparator<MutableTuple<Tree, Integer>>() {
            @Override
            public int compare(MutableTuple<Tree, Integer> o1, MutableTuple<Tree, Integer> o2) {
                return o2.Item2 - o1.Item2;
            }
        });
        for(MutableTuple<Tree, Integer> tup : list) {
            sb.append(tup.Item2 + "\t" + tup.Item1.toNewick() + "\n");
        }
        return sb.toString();
    }

    public String summarizeNetTopos() {
        StringBuilder sb = new StringBuilder();
        List<MutableTuple<Network, Integer>> list = new ArrayList<>();
        try {
            for(String s : _samples) {
                Network temp = Networks.readNetwork(s);
                boolean found = false;
                for(MutableTuple<Network, Integer> tup : list) {
                    if(Networks.hasTheSameTopology(temp, tup.Item1)) {
                        found = true;
                        tup.Item2 = tup.Item2 + 1;
                    }
                }
                if(!found) {
                    list.add(new MutableTuple<Network, Integer>(temp, 1));
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        Collections.sort(list, new Comparator<MutableTuple<Network, Integer>>() {
            @Override
            public int compare(MutableTuple<Network, Integer> o1, MutableTuple<Network, Integer> o2) {
                return o2.Item2 - o1.Item2;
            }
        });
        for(MutableTuple<Network, Integer> tup : list) {
            sb.append(tup.Item2 + "\t" + tup.Item1.toString() + "\n");
        }
        return sb.toString();
    }
}