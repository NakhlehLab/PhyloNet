package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment.MarkerSeq;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimSNPInNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;

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
    ConcurrentHashMap<String, DoubleAdder> _counter;
    Map<String, Double> _counter0;
    boolean _onlyPolymorphic;
    int _simnum = 10000;
    BiAllelicGTR _BAGTRModel;
    boolean _diploid;
    boolean _dominant;
    Integer _polyploid = null;

    public ApproximateBayesian(Map<String, String> alleles2species, List<MarkerSeq> markerSeqs, boolean diploid, boolean dominant, Integer polyploid, boolean onlyPolymorphic, BiAllelicGTR BAGTRModel) {
        if(diploid && !dominant) {
            //throw new RuntimeException("Not support!");
        }

        _diploid = diploid;
        _dominant = dominant;
        _polyploid = polyploid;
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
        scale = 1.0 * _simnum / markerSeqs.get(0).getSiteCount();
        for(MarkerSeq aln : markerSeqs) {
            for (int i = 0; i < aln.getSiteCount(); i++) {
                int j = 0;
//                BitSet key = new BitSet(_alleles.size() * 2);
//                for(String allele : _alleles) {
//                    char c = aln.getAlignment().get(allele).charAt(i);
//                    if(c == '1') {
//                        key.set(j);
//                    }
//                    j++;
//                }
//                if(!_counter0.containsKey(key)) {
//                    _counter0.put(key, 0.0);
//                }
//                _counter0.put(key, _counter0.get(key) + scale);

                for(int j1 = 0 ; j1 < _alleles.size() ; j1++) {
                    for(int j2 = j1 + 1 ; j2 < _alleles.size() ; j2++) {
                        for(int j3 = j2 + 1 ; j3 < _alleles.size() ; j3++) {
                            char c1 = aln.getAlignment().get(_alleles.get(j1)).charAt(i);
                            char c2 = aln.getAlignment().get(_alleles.get(j2)).charAt(i);
                            char c3 = aln.getAlignment().get(_alleles.get(j3)).charAt(i);
                            StringBuilder keybuilder = new StringBuilder();
                            for(int ii = 0 ; ii < _alleles.size() * 2 ; ii++) keybuilder.append("0");

                            //BitSet key = new BitSet(_alleles.size() * 2);
                            keybuilder.setCharAt(j1, '1');
                            keybuilder.setCharAt(j2, '1');
                            keybuilder.setCharAt(j3, '1');

//                            key.flip(j1);
//                            key.flip(j2);
//                            key.flip(j3);
//
//                            if(c1 == '1') key.flip(j1 + _alleles.size());
//                            if(c2 == '1') key.flip(j2 + _alleles.size());
//                            if(c3 == '1') key.flip(j3 + _alleles.size());

                            keybuilder.setCharAt(j1 + _alleles.size(), c1);
                            keybuilder.setCharAt(j2 + _alleles.size(), c2);
                            keybuilder.setCharAt(j3 + _alleles.size(), c3);

                            String key = keybuilder.toString();

                            if(!_counter0.containsKey(key)) {
                                _counter0.put(key, 0.0);
                            }
                            _counter0.put(key, _counter0.get(key) + scale);
                        }
                    }
                }



//                for(int k = 0 ; k < _alleles.size() ; k++) {
//                    BitSet newkey = (BitSet) key.clone();
//                    if(k == -1) {
//                        newkey.flip(0, _alleles.size() - 1);
//                    } else {
//                        newkey.flip(k);
//                    }
//                    if(!_counter0.containsKey(newkey)) {
//                        _counter0.put(newkey, 0.0);
//                    }
//                    _counter0.put(newkey, _counter0.get(newkey) + scale);
//                }
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

        //Utils._SEED = new Random().nextLong();//4385882401632407363L;
        //System.out.println(Utils._SEED);
        _simulators = new ArrayList<>();
        for(int i = 0 ; i < Utils._NUM_THREADS; i++) {
            _simulators.add(new SimSNPInNetwork(_BAGTRModel, Utils._SEED + i * 10000));
            _simulators.get(i)._diploid = _diploid;
            _simulators.get(i)._dominant = _dominant;
            _simulators.get(i)._polyploid = _polyploid;
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
//                        BitSet key = new BitSet(_alleles.size());
//                        int j = 0;
//                        for (String allele : _alleles) {
//                            char c = onesnp.get(allele).charAt(0);
//                            if (c == '1') {
//                                key.set(j);
//                            }
//                            j++;
//                        }
//                        _counter.putIfAbsent(key, new DoubleAdder());
//                        _counter.get(key).add(1.0);

                        for(int j1 = 0 ; j1 < _alleles.size() ; j1++) {
                            for(int j2 = j1 + 1 ; j2 < _alleles.size() ; j2++) {
                                for(int j3 = j2 + 1 ; j3 < _alleles.size() ; j3++) {
                                    char c1 = onesnp.get(_alleles.get(j1)).charAt(0);
                                    char c2 = onesnp.get(_alleles.get(j2)).charAt(0);
                                    char c3 = onesnp.get(_alleles.get(j3)).charAt(0);

                                    StringBuilder keybuilder = new StringBuilder();
                                    for(int ii = 0 ; ii < _alleles.size() * 2 ; ii++) keybuilder.append("0");

                                    keybuilder.setCharAt(j1, '1');
                                    keybuilder.setCharAt(j2, '1');
                                    keybuilder.setCharAt(j3, '1');

                                    keybuilder.setCharAt(j1 + _alleles.size(), c1);
                                    keybuilder.setCharAt(j2 + _alleles.size(), c2);
                                    keybuilder.setCharAt(j3 + _alleles.size(), c3);

                                    String key = keybuilder.toString();
//                                    BitSet key = new BitSet(_alleles.size() * 2);

//                                    key.flip(j1);
//                                    key.flip(j2);
//                                    key.flip(j3);
//
//                                    if(c1 == '1') key.flip(j1 + _alleles.size());
//                                    if(c2 == '1') key.flip(j2 + _alleles.size());
//                                    if(c3 == '1') key.flip(j3 + _alleles.size());

                                    _counter.putIfAbsent(key, new DoubleAdder());
                                    _counter.get(key).add(1.0);
                                }
                            }
                        }

//                        for(int kk = 0 ; kk < _alleles.size() ; kk++) {
//                            BitSet newkey = (BitSet) key.clone();
//                            if(kk == -1) {
//                                newkey.flip(0, _alleles.size() - 1);
//                            } else {
//                                newkey.flip(kk);
//                            }
//
//                            _counter.putIfAbsent(newkey, new DoubleAdder());
//                            _counter.get(newkey).add(1.0);
//                        }
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
        Set<String> keys = new HashSet<>();
        keys.addAll(_counter.keySet());
        keys.addAll(_counter0.keySet());

        for(String key : keys) {
            double value0 = _counter0.getOrDefault(key, 0.0);
            DoubleAdder value = _counter.getOrDefault(key, new DoubleAdder());
            double value1 = value.sum() / modifier;
            if(value1 == 0.0) value1 = 1.0;
//            if(value1 / _simnum < 0.01) {
//                continue;
//            }

            sum += (value0 - value1) * (value0 - value1) / value1;
        }

        return -sum;
    }
}
