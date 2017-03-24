package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.summary;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.MutableTuple3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.summary.SummaryBranch;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.*;

/**
 * Created by wendingqiao on 7/20/15.
 */
public class Convergence {

    private List<String> _files;

    private double threshold = 0.95;
    private List<List<Summ>> summaryLists;

    private Map<Network, Double> networkCounts = new HashMap<>();
    private PriorityQueue<MutableTuple<Double,Network>> pqueue;
    private List<MutableTuple<Double,Network>> pqList = new ArrayList<>();


    public Convergence(List<String> files) {
        _files = files;
        summaryLists = new ArrayList<>();
        for(String f : _files) {
            processFile(f);
        }
        pqueue = new PriorityQueue<>(networkCounts.size(), new TupleComparator());
        for(Network key : networkCounts.keySet()) {
            double percent = networkCounts.get(key) / summaryLists.size() / summaryLists.get(0).size();
            pqueue.add(new MutableTuple<Double, Network>(percent, key));
        }
        double tmp = 0.0;
        while(!pqueue.isEmpty()) {
            pqList.add(pqueue.poll());
            tmp += pqList.get(pqList.size()-1).Item1;
            System.out.printf("%.4f - %s\n", pqList.get(pqList.size()-1).Item1, pqList.get(pqList.size()-1).Item2);
            if(tmp > threshold) break;
        }
    }

    public void summarizeTopo(boolean withBL) {
        for(MutableTuple<Double, Network> tup : pqList) {
            if(withBL) {
                SummaryBranch sum = new SummaryBranch(tup.Item2.toString());
                for(String f : _files) {
                    sum.addFile(f);
                }
                sum.report(1, 1);
            } else {
                SummaryBL sum = new SummaryBL(tup.Item2.toString());
                for(String f: _files) {
                    sum.addFile(f);
                }
                sum.report();
            }
        }
    }


    public void computePSRF() {
        List<List<Double>> vars = new ArrayList<>(); // each file, each param  E(X^2)
        List<List<Double>> aves = new ArrayList<>(); // each file, each param  E(X)
        double nsize = summaryLists.get(0).size(); // number of samples
        double msize = summaryLists.size(); // number of files

        for(List<Summ> list : summaryLists) { // each file
            List<Double> subvar = new ArrayList<>();
            List<Double> subave = new ArrayList<>();
            for( int i = 0; i < list.get(0).params.size(); i++ ) { // each param
                double var = 0, ave = 0;
                for(Summ su : list) { // each sample
                    var += su.params.get(i) / nsize * su.params.get(i);
                    ave += su.params.get(i) / nsize;
                }
                subvar.add(var);
                subave.add(ave);
            }
            vars.add(subvar);
            aves.add(subave);
        }

        List<Double> wvars = new ArrayList<>(); // each param  W
        List<Double> waves = new ArrayList<>(); // each param  ave-ave
        List<Double> bvars = new ArrayList<>(); // each param  B
        for(int i = 0; i < vars.get(0).size(); i++ ) { // for each param
            double wv = 0, wa = 0;
            for(int j = 0; j < vars.size(); j++) { // for each file
                wv += vars.get(j).get(i) - aves.get(j).get(i) * aves.get(j).get(i); // sj^2
                wa += aves.get(j).get(i);
            }
            wvars.add(wv * nsize / (nsize-1) / msize); // w
            waves.add(wa / msize); // ave-ave
        }
        for(int i = 0; i < vars.get(0).size(); i++) { // for each param
            double bv = 0;
            for(int j = 0; j < vars.size(); j++) {
                bv += (aves.get(j).get(i) - waves.get(i)) * (aves.get(j).get(i) - waves.get(i));
            }
            bvars.add(bv / (msize-1));
        }
        System.out.println("---------------- PSRF ---------------- ");
        for(int i = 0; i < wvars.size(); i++) { // each param
            double res = Math.sqrt( 1 - 1/nsize + bvars.get(i)/wvars.get(i) );
            System.out.println("PSRF - " + res);
        }
        System.out.println();
    }


    public void generateTracePlot() {
        System.out.println("---------------- Trace Plot ---------------- ");
        for(List<Summ> list : summaryLists) {
            generateTracePlot(list);
        }
    }


    public void generateSRQ() {
        List<MutableTuple<Network, MutableTuple3<List<Double>,List<Double>,List<Double> >>> coutsList = new ArrayList<>();
        for(MutableTuple<Double,Network> mt : pqList) {
            MutableTuple3<List<Double>, List<Double>, List<Double>> mt3 = new MutableTuple3<List<Double>, List<Double>, List<Double>>(
                    new ArrayList<Double>(), new ArrayList<Double>(), new ArrayList<Double>());
            coutsList.add(new MutableTuple<Network, MutableTuple3<List<Double>, List<Double>, List<Double>>>(
                    mt.Item2, mt3 ) );
        }
        for(List<Summ> list : summaryLists) {
            generateSRQ(list, coutsList);
        }
        System.out.println("---------------- Summary ---------------- ");
        for(int fidx = 0; fidx < summaryLists.size(); fidx++) {
            System.out.printf("Topology  Frequency Posterior #sojourns max(sjs)  ave(sjs)\n");
            int tidx = 0;
            for(MutableTuple<Network, MutableTuple3<List<Double>,List<Double>,List<Double> >> mt : coutsList) {
                StringBuilder sb = new StringBuilder("");
                sb.append(String.format("%6s   &", ++tidx));
                sb.append(String.format("%7.0f &", mt.Item2.Item1.get(fidx)));
                sb.append(String.format("  %1.4f &", mt.Item2.Item1.get(fidx) / summaryLists.get(fidx).size() ));
                sb.append(String.format("%7.0f &", mt.Item2.Item2.get(fidx) ));
                sb.append(String.format("%7.0f &", mt.Item2.Item3.get(fidx)));
                sb.append(String.format("    %1.2f \\\\", mt.Item2.Item1.get(fidx) / mt.Item2.Item2.get(fidx) ));
                System.out.println(sb.toString());
            }
        }
        System.out.println("---------------- SRQ Plot ---------------- ");
        for(List<Summ> list : summaryLists) {
            generateSRQPlot( list, pqList.get(0).Item2 );
        }
    }


    private void generateTracePlot( List<Summ> list ) {
        for(int i = 0; i < list.get(0).params.size(); i++) {
            StringBuilder sby = new StringBuilder("y" + i + " = [");
            for(Summ s : list) {
                sby.append(s.params.get(i) + ",");
            }
            sby.append("];");
            System.out.println(sby.toString());
        }
    }


    private void generateSRQPlot( List<Summ> list, Network net) {
        StringBuilder sbx = new StringBuilder("x = [");
        StringBuilder sby = new StringBuilder("y = [");
        int x = 0, y = 0;
        for(Summ s : list) {
            x++;
            if (Networks.hasTheSameTopology(net, Networks.readNetwork(s.net))) {
                y++;
                sbx.append(x + ",");
                sby.append(y + ",");
            }
        }
        sbx.deleteCharAt(sbx.length()-1);
        sby.deleteCharAt(sby.length()-1);
        sbx.append("];");
        sby.append("];");
        System.out.println(sbx);
        System.out.println(sby);
    }


    private void generateSRQ(List<Summ> list,
                             List<MutableTuple<Network, MutableTuple3<List<Double>,List<Double>,List<Double>>>> coutsList)
    {
        Map<Network, Double> netcouts = new HashMap<>();
        Map<Network, Double> sojournsMap = new HashMap<>();
        Map<Network, Double> sojournsCount = new HashMap<>();
        for(MutableTuple<Double, Network> mt : pqList) {
            netcouts.put(mt.Item2, 0.0);
            sojournsMap.put(mt.Item2, 0.0);
            sojournsCount.put(mt.Item2, 0.0);
        }
        int maxSojourn = 1;

        Network prev = null;
        for(Summ s : list) {
            Network test = Networks.readNetwork(s.net);

            boolean sj = false;
            if(prev != null && Networks.hasTheSameTopology(test, prev)) {
                maxSojourn++;
            } else {
                maxSojourn = 1;
                sj = true;
            }
            prev = test;

            for(Network net : netcouts.keySet()) {
                if(Networks.hasTheSameTopology(net, test)) {
                    netcouts.put(net, netcouts.get(net) + 1.0);
                    sojournsMap.put(net, Math.max(sojournsMap.get(net), maxSojourn));
                    if(sj) sojournsCount.put(net, sojournsCount.get(net)+1.0);
                    break;
                }
            }

        }
        for(MutableTuple<Network, MutableTuple3<List<Double>,List<Double>,List<Double>>> mt : coutsList) {
            if(!netcouts.containsKey(mt.Item1)){
                mt.Item2.Item1.add(0.0);
                mt.Item2.Item2.add(0.0);
                mt.Item2.Item3.add(0.0);
            } else {
                mt.Item2.Item1.add( netcouts.get(mt.Item1) );
                mt.Item2.Item2.add( sojournsCount.get(mt.Item1) );  // count
                mt.Item2.Item3.add( sojournsMap.get(mt.Item1)); // map
            }
        }
    }


    private void processFile(String file) {
        List<Summ> list = new ArrayList<>();
        try{
            BufferedReader br = new BufferedReader((new FileReader(file)));
            br.readLine(); br.readLine();
            String s;
            String[] ss;
            boolean burnin = true;
            while( (s = br.readLine()) != null) {
                if(s.contains("Summarization")) break;
                ss = s.split(";");
                if( ss.length != 7 || ss[0].contains("Iteration") || !Character.isDigit(ss[0].charAt(0))) {
                    continue;
                }
                if(Double.parseDouble(ss[2]) != 0 && Double.parseDouble(ss[5]) != 0) {
                    burnin = false;
                }
                if(!burnin) {
                    br.readLine();
                    continue;
                }
                List<Double> params = new ArrayList<>();
                params.add(Double.parseDouble(ss[1]));
                params.add(Double.parseDouble(ss[3]));
                params.add(Double.parseDouble(ss[4]));
                String net = br.readLine();
                addNet(net);
                Summ summ = new Summ( params, net );
                list.add(summ);
            }
            br.close();

        } catch (Exception ex) {
            ex.printStackTrace();
        }
        System.out.println(list.size());
        summaryLists.add(list);
    }

    private void addNet(String n) {
        Network net = Networks.readNetwork(n);
        boolean found = false;
        for(Network key : networkCounts.keySet()) {
            if(Networks.hasTheSameTopology(key, net)) {
                networkCounts.put( key, networkCounts.get(key)+1.0 );
                found = true;
                break;
            }
        }
        if(!found) networkCounts.put(net, 1.0);
    }


    class TupleComparator implements Comparator<MutableTuple<Double, Network>> {
        @Override
        public int compare(MutableTuple<Double, Network> x, MutableTuple<Double, Network> y) {
            return (x.Item1 < y.Item1) ? 1 : -1;
        }
    }


    class Summ{

        List<Double> params;
        String net;

        public Summ(List<Double> list, String s) {
            params = list;
            net = s;
        }
    }
}
