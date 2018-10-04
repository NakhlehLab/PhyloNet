package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimGTInNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimSNPInNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.DoubleAdder;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 2/5/18
 * Time: 11:18 AM
 * To change this template use File | Settings | File Templates.
 */
public class ApproximateBayesian {

    List<SimSNPInNetwork> _simulators;
    List<String> _alleles;
    Map<String, List<String>> _species2alleles;
    Map<String, String> _alleles2species;
    ConcurrentHashMap<BitSet, DoubleAdder> _counter;
    Map<BitSet, Double> _counter0;
    boolean _onlyPolymorphic;
    int _simnum = 50000;
    BiAllelicGTR _BAGTRModel;
    boolean _diploid;
    boolean _dominant;

    public ApproximateBayesian(Map<String, String> alleles2species, List<Alignment> alignments, boolean diploid, boolean dominant, boolean onlyPolymorphic, BiAllelicGTR BAGTRModel) {
        if(diploid && !dominant) {
            //throw new RuntimeException("Not support!");
        }

        _diploid = diploid;
        _dominant = dominant;
        _BAGTRModel = BAGTRModel;
        _simulators = new ArrayList<>();
        for(int i = 0 ; i < Utils._NUM_THREADS; i++) {
            _simulators.add(new SimSNPInNetwork(BAGTRModel, Utils._SEED + i * 10000));
            _simulators.get(i)._diploid = diploid;
            _simulators.get(i)._dominant = dominant;
        }

        _onlyPolymorphic = onlyPolymorphic;
        if(_simnum % _simulators.size() != 0) {
            _simnum -= _simnum % _simulators.size();
        }
        _alleles = new ArrayList<>();
        for(String allele : alleles2species.keySet()) {
            _alleles.add(allele);
        }
        Collections.sort(_alleles);

        _counter0 = new HashMap<>();
        double scale = 0.0;
        scale = 1.0 * _simnum / alignments.get(0).getSiteCount();
        for(Alignment aln : alignments) {
            for (int i = 0; i < aln.getSiteCount(); i++) {
                int j = 0;
                BitSet key = new BitSet(_alleles.size());
                for(String allele : _alleles) {
                    char c = aln.getAlignment().get(allele).charAt(i);
                    if(c == '1') {
                        key.set(j);
                    }
                    j++;
                }
                if(!_counter0.containsKey(key)) {
                    _counter0.put(key, 0.0);
                }
                _counter0.put(key, _counter0.get(key) + scale);
            }
        }

        _alleles2species = alleles2species;
        _species2alleles = new HashMap<>();
        for(String allele : _alleles) {
            String species = alleles2species.get(allele);
            if(!_species2alleles.containsKey(species)) {
                _species2alleles.put(species, new ArrayList<>());
            }
            _species2alleles.get(species).add(allele);
        }
    }

    public double computeApproximateBayesianMT(Network network) {
        ExecutorService executor = Executors.newFixedThreadPool(Utils._NUM_THREADS);
        DoubleAdder adder = new DoubleAdder();
        _counter = new ConcurrentHashMap<>();

        _simulators = new ArrayList<>();
        for(int i = 0 ; i < Utils._NUM_THREADS; i++) {
            _simulators.add(new SimSNPInNetwork(_BAGTRModel, Utils._SEED + i * 10000));
            _simulators.get(i)._diploid = _diploid;
            _simulators.get(i)._dominant = _dominant;
        }

        //System.out.println(Networks.getFullString(network));
        //network = Networks.readNetworkWithRootPop("[0.006]((R:0.5128354928838631,C:0.5128354928838631)II1:0.1996645071161368,(L:0.2912372832655849,(A:0.18191724539589205,Q:0.18191724539589205)II3:0.10932003786969285)II2:0.421262716734415)II0;");

        if(Utils._CONST_POP_SIZE) {
            // theta = network.getRoot().getRootPopSize();
            for(Object childObj : network.bfs()) {
                NetNode child = (NetNode) childObj;
                for(Object parentObj : child.getParents()) {
                    NetNode parent = (NetNode) parentObj;
                    child.setParentSupport(parent, network.getRoot().getRootPopSize());
                }
            }
        }

        for(int i = 0 ; i < _simulators.size() ; i++) {
            final int start = _simnum / _simulators.size() * i;
            final int index = i;
            final int perthread = _simnum / _simulators.size();
            final Network net = network.clone();
            executor.execute(new Runnable() {
                public void run() {
                    for(int k = 0 ; k < perthread ; k++) {
                        //System.out.println("index: " + index + " k: " + k);
                        Map<String, String> onesnp;
                        Map<String, List<String>> species2alleles = new HashMap<>();
                        for (String allele : _alleles) {
                            String species = _alleles2species.get(allele);
                            if (!species2alleles.containsKey(species)) {
                                species2alleles.put(species, new ArrayList<>());
                            }
                            species2alleles.get(species).add(allele);
                        }
                        Network net1 = net.clone();
                        onesnp = _simulators.get(index).generateSNPs(net1, species2alleles, 1, true);
                        BitSet key = new BitSet(_alleles.size());
                        int j = 0;
                        for (String allele : _alleles) {
                            char c = onesnp.get(allele).charAt(0);
                            if (c == '1') {
                                key.set(j);
                            }
                            j++;
                        }
                        _counter.putIfAbsent(key, new DoubleAdder());
                        _counter.get(key).add(1.0);
                    }
                }
            });
        }

        try {
            executor.shutdown();
            executor.awaitTermination(1000, TimeUnit.SECONDS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        double modifier = 1.0;
        if(_onlyPolymorphic) {
            BitSet key0 = new BitSet(_alleles.size());
            BitSet key1 = new BitSet(_alleles.size());
            key1.set(0, _alleles.size());
            if(_counter.containsKey(key0)) {
                modifier += _counter.get(key0).sum();
                _counter.remove(key0);
            }
            if(_counter.containsKey(key1)) {
                modifier += _counter.get(key1).sum();
                _counter.remove(key1);
            }
            modifier /= _simnum;
            modifier = 1.0 - modifier;
        }


        double sum = 0.0;
        Set<BitSet> keys = new HashSet<>();
        keys.addAll(_counter0.keySet());
        keys.addAll(_counter.keySet());

        for(BitSet key : keys) {
            double value0 = _counter0.getOrDefault(key, 0.0);
            DoubleAdder value = _counter.getOrDefault(key, new DoubleAdder());
            double value1 = value.sum() / modifier;
            sum += (value0 - value1) * (value0 - value1);
        }

        return -sum;
    }
}
