package edu.rice.cs.bioinfo.programs.phylonet.algos.network;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;

/**
 * Created by zhiyan on 7/5/20.
 */
public class MDCOnAllopolyploidNetwork {
    boolean _printDetail = false;
    boolean _parallel = false;
    int _currentTreeID = -1;

    /**
     * This function is for setting the option of printing details
     */
    public void setParallel(boolean parallel) { _parallel = parallel; }


    /**
     * Returns the next gene tree to compute, which is used for parallel computing
     */
    public synchronized int getNextTreeID() { return ++_currentTreeID; }


    /**
     * This function is for setting the option of printing details
     */
    public void setPrintDetails(boolean p){
        _printDetail = p;
    }


    /**
     * This is the main function for counting the minimal number of extra lineages required for reconciling a collection of gene trees within the branches of a given species network
     *
     * @param	net 	the given species network
     * @param 	gts		the collection of gene trees
     * @param	alleles2species		the mapping from alleles to the species they are sampled from
     *
     * @return	the resulting number of extra lineages which have one-to-one correspondence with gene trees in gts
     */
    public void countExtraCoal(Network net, List<MutableTuple<Tree,Double>> gts, Map<String, String> alleles2species, int[] xls){
        List<Integer> xlList = new ArrayList<Integer>();
        Map<String,Integer> nname2tamount = new HashMap<String,Integer>();
        Tree superst = networkToTree(net, nname2tamount);
        Map<String,String> tname2nname = new HashMap<String,String>();

        for(Map.Entry<String, Integer> entry: nname2tamount.entrySet()){
            if(entry.getValue()>1){
                for(int i=1; i<=entry.getValue(); i++){
                    tname2nname.put(entry.getKey()+"_"+i, entry.getKey());
                }
            }
        }

        int treeID = 0;
        if(_parallel) {
            treeID = getNextTreeID();
        }

        while(treeID < gts.size()){
            MutableTuple<Tree, Double> tuple = gts.get(treeID);
            Tree gt = tuple.Item1;
            List<String> gtTaxa = Arrays.asList(gt.getLeaves());
            if(alleles2species == null){
                alleles2species = new HashMap<String,String>();
                for(String taxon: gtTaxa){
                    alleles2species.put(taxon, taxon);
                }
            }

            List<String> netTaxa = new ArrayList<String>();
            List<List<String>> allelesList = new ArrayList<List<String>>();
            List<Integer> upper = new ArrayList<Integer>();
            for(String gtleaf: gtTaxa){
                String nleaf = alleles2species.get(gtleaf);
                int index = netTaxa.indexOf(nleaf);
                if(index==-1){
                    netTaxa.add(nleaf);
                    List<String> alleles = new ArrayList<String>();
                    alleles.add(gtleaf);
                    allelesList.add(alleles);
                    upper.add(nname2tamount.get(nleaf));
                }
                else{
                    allelesList.get(index).add(gtleaf);
                }
            }

            List<int[]> mergeNumber = new ArrayList<int[]>();
            for(List<String> alleles: allelesList){
                int[] first = new int[alleles.size()];
                Arrays.fill(first, 1);
                mergeNumber.add(first);
            }
            int minCoal = Integer.MAX_VALUE;
            do{
                Map<String, String> mapping = new HashMap<String,String>();
                for(int i=0; i<netTaxa.size(); i++){
                    String baseName = netTaxa.get(i);
                    List<String> alleles = allelesList.get(i);
                    int[] subscribes = mergeNumber.get(i);
                    for(int j=0; j<alleles.size(); j++){
                        mapping.put(alleles.get(j), baseName+"_"+subscribes[j]);
                    }
                }
                int coal = countExtraCoal(gt, superst, mapping);
                if(_printDetail){
                    System.out.println(mapping + ":  " + coal);
                }
                minCoal = Math.min(coal, minCoal);

            }while(mergeNumberAddOne(mergeNumber,upper));

            xls[treeID] = minCoal;
            if (_parallel) {
                treeID = getNextTreeID();
            }
            else {
                treeID++;
            }
        }
        return;
    }



    /**
     * This function is to help enumerating allele mappings
     */
    private boolean mergeNumberAddOne(List<int[]> mergeNumber, List<Integer> upper){
        for(int i=0; i<mergeNumber.size(); i++){
            int[] partNumber = mergeNumber.get(i);
            int max = upper.get(i);
            for(int j=0; j<partNumber.length; j++){
                if(partNumber[j]==max){
                    partNumber[j] = 1;
                }
                else{
                    partNumber[j] = partNumber[j]+1;
                    return true;
                }
            }
            Arrays.fill(partNumber, 1);
        }
        return false;
    }



    /**
     * The function is to convert a network to a multilabel tree
     *
     * @param net 	            the species network to convert
     * @param nname2tamount     mappings from taxa in the species network to the number of copies of that leaf in the resulting mul-trees
     */
    private Tree networkToTree(Network<Double> net, Map<String,Integer> nname2tamount){
        removeBinaryNodes(net);
        Tree st = new STITree<Double>();
        Queue<NetNode<Double>> source = new LinkedList<NetNode<Double>>();
        Queue<TMutableNode> dest = new LinkedList<TMutableNode>();
        source.offer(net.getRoot());
        dest.offer((TMutableNode)st.getRoot());
        long nameid = System.currentTimeMillis();
        while(!source.isEmpty()){
            NetNode<Double> parent = source.poll();
            TMutableNode peer = dest.poll();
            for (NetNode<Double> child : parent.getChildren()) {
                TMutableNode copy;
                if (child.getName().equals(NetNode.NO_NAME)) {
                    child.setName("hnode" + (nameid++));
                }
                String name = child.getName();
                if(child.isNetworkNode()){
                    name = child.getName()+"TO"+parent.getName();
                }
                Integer amount = nname2tamount.get(name);
                if(amount==null){
                    amount = 0;
                }
                nname2tamount.put(name, ++amount);
                String newname = name + "_" + amount;
                copy = peer.createChild(newname);

                double distance = child.getParentDistance(parent);
                if (distance == NetNode.NO_DISTANCE) {
                    copy.setParentDistance(0);
                }
                else {
                    copy.setParentDistance(distance);
                }

                source.offer(child);
                dest.offer(copy);
            }
        }
        return st;
    }




    /**
     * This function is to count the minimal number of extra lineages required for reconciling a gene tree within the branches of a given species network
     *
     * @param gt   	        the given gene tree
     * @param st		    the multree
     * @param taxonMap	    the mapping from alleles to the species they are sampled from
     *
     * @return	the number of extra lineages
     */
    private int countExtraCoal(Tree gt, Tree st, Map<String, String> taxonMap){
        int sum = 0;
        String[] stTaxa = st.getLeaves();

        Map<TNode, BitSet> map = new HashMap<TNode, BitSet>();
//        Map<String, Integer> nname2xl = new HashMap<String, Integer>();

        for (TNode node : st.postTraverse()) {
            BitSet bs = new BitSet();
            if (node.isLeaf()) {
                for (int i = 0; i < stTaxa.length; i++) {
                    if (node.getName().equals(stTaxa[i])) {
                        bs.set(i);
                        break;
                    }
                }
                map.put(node, bs);
            }
            else {
                for (TNode child : node.getChildren()) {
                    BitSet childCluster = map.get(child);
                    bs.or(childCluster);
                }
                map.put(node, bs);
            }

            if(node.getChildCount()==1 && node.getParentDistance()==0){
                continue;
            }

            STITreeCluster c = new STITreeCluster(stTaxa);
            c.setCluster(bs);

            if(c.getClusterSize()<stTaxa.length){
                int el = getClusterCoalNum(gt, c, taxonMap);
                sum += Math.max(0, el-1);
//                String tname = node.getName();

//                if(tname!=null && tname2nname.containsKey(tname)){
//                    String nname = tname2nname.get(tname);
//                    Integer xl = nname2xl.get(nname);
//                    if(xl==null){
//                        xl = 0;
//                    }
//                    xl += el;
//                    nname2xl.put(nname, xl);
//                }
//                else{
//                    sum += Math.max(0, el-1);
//                }
            }
        }

//        for(Map.Entry<String, Integer> entry: nname2xl.entrySet()){
//            sum += Math.max(0, entry.getValue()-1);
//        }

        return sum;
    }



    /**
     * This function is to count the minimal number of extra lineages given a gene tree and a cluster in the species tree
     *
     * @param tr   	        the given gene tree
     * @param cluster		the cluster in the species tree
     * @param taxonMap	    the mapping from alleles to the species they are sampled from
     *
     * @return	the number of extra lineages
     */
    private int getClusterCoalNum(Tree tr, STITreeCluster cluster, Map<String, String> taxonMap) {
        Map<TNode, BitSet> map = new HashMap<TNode, BitSet>();
        List<String> taxa = new LinkedList<String>();	// List of species taxa.

        Collections.addAll(taxa, cluster.getTaxa());

        int count = 0;
        for (TNode node : tr.postTraverse()) {
            if (node.isLeaf()) {
                String stTaxon = taxonMap.get(node.getName());	// Get the corresponding species name.
                int index = taxa.indexOf(stTaxon);
                BitSet bs = new BitSet(taxa.size());
                bs.set(index);
                if (cluster.containsCluster(bs)) {
                    count++;
                }
                map.put(node, bs);
            }
            else {
                BitSet bs = new BitSet(taxa.size());
                int intersect = 0;
                int childCount = node.getChildCount();
                for (TNode child : node.getChildren()) {
                    BitSet v = map.get(child);
                    bs.or(v);
                    if(childCount>2){
                        if(cluster.containsCluster(v)){
                            intersect ++;
                        }
                    }
                }

                if (cluster.containsCluster(bs)) {
                    count -= node.getChildCount();
                    count++;
                }
                else if(intersect>1){
                    count -= intersect;
                    count ++;
                }

                map.put(node, bs);
            }
        }
        return count;
    }




    /**
     * The function is to remove binary nodes of a species network
     */
    private void removeBinaryNodes(Network<Double> net)
    {
        // Find all binary nodes.
        List<NetNode<Double>> binaryNodes = new LinkedList<NetNode<Double>>();
        for (NetNode<Double> node : net.bfs()) {
            if (node.getIndeg() == 1 && node.getOutdeg() == 1) {
                binaryNodes.add(node);
            }
        }

        // Remove them.
        for (NetNode<Double> node : binaryNodes) {
            NetNode<Double> child = node.getChildren().iterator().next();	// Node's only child.
            if(child.getIndeg() != 1){
                continue;
            }
            NetNode<Double> parent = node.getParents().iterator().next();	// Node's only parent.
            double distance = node.getParentDistance(parent) + child.getParentDistance(node);
            double gamma = node.getParentProbability(parent) * child.getParentProbability(node);
            parent.removeChild(node);
            node.removeChild(child);
            parent.adoptChild(child, distance);
            child.setParentProbability(parent, gamma);
        }
    }

    private static void testCountExtraCoal1() {
        Tree geneTree = Trees.readTree("((a,b),(c,d));");
        List<MutableTuple<Tree,Double>> geneTrees = new ArrayList<>();
        geneTrees.add(new MutableTuple<>(geneTree, 0.0));

        Network speciesNetwork = Networks.readNetwork("(((B,C)#H1,A),(#H1,D));");
        Map<String, String> allele2species = new HashMap<>();
        allele2species.put("a", "A");
        allele2species.put("b", "B");
        allele2species.put("c", "C");
        allele2species.put("d", "D");

        MDCOnAllopolyploidNetwork mdc = new MDCOnAllopolyploidNetwork();
        int xls[] = new int[geneTrees.size()];
        mdc.countExtraCoal(speciesNetwork, geneTrees, allele2species, xls);
        System.out.println(Arrays.toString(xls));
    }

    private static void countExtraCoalBioData() {
        boolean useVanillaMdc = false;
        List<String> spNewicks = Lists.newArrayList(
                "(Castrilanthemum_debeauxii,((((Leucanthemopsis_pulverulenta,(Leucanthemopsis_pallida_subsp_spathulifolia,Leucanthemopsis_pallida_var_bilbilitanum)),((Leucanthemopsis_alpina_subsp_cuneata,Leucanthemopsis_pallida_var_virescens),Leucanthemopsis_pectinata)),(((Leucanthemopsis_alpina4x,Leucanthemopsis_alpina),Leucanthemopsis_alpina_subsp_tomentosa),Leucanthemopsis_longipectinata)),(Hymenostemma_pseudanthemis,Prolongoa_hispanica)));",
                "(Castrilanthemum_debeauxii,(((((Leucanthemopsis_alpina4x,Leucanthemopsis_alpina),Leucanthemopsis_alpina_subsp_tomentosa),((Leucanthemopsis_alpina_subsp_cuneata)#H1,Leucanthemopsis_longipectinata)),(((Leucanthemopsis_pectinata,Leucanthemopsis_pallida_subsp_spathulifolia),(Leucanthemopsis_pallida_var_virescens,(Leucanthemopsis_pallida_var_bilbilitanum,#H1))),Leucanthemopsis_pulverulenta)),(Prolongoa_hispanica,Hymenostemma_pseudanthemis)));",
                "((((((Leucanthemopsis_pallida_subsp_spathulifolia,Leucanthemopsis_pectinata),(Leucanthemopsis_pallida_var_virescens,(Leucanthemopsis_pallida_var_bilbilitanum,(((Leucanthemopsis_longipectinata,Leucanthemopsis_alpina_subsp_cuneata))#H2)#H1))),Leucanthemopsis_pulverulenta),(#H1,(Leucanthemopsis_alpina_subsp_tomentosa,(Leucanthemopsis_alpina4x,Leucanthemopsis_alpina)))),(Hymenostemma_pseudanthemis,Prolongoa_hispanica)),(Castrilanthemum_debeauxii,#H2));",
                "(Castrilanthemum_debeauxii,((Hymenostemma_pseudanthemis,Prolongoa_hispanica)I3,(((Leucanthemopsis_alpina_subsp_cuneata)I8#H1,Leucanthemopsis_longipectinata)I5,((((Leucanthemopsis_pallida_subsp_spathulifolia,(I8#H1,(Leucanthemopsis_pallida_var_bilbilitanum,Leucanthemopsis_pulverulenta)I14)I13)I12)I11,(Leucanthemopsis_pallida_var_virescens,Leucanthemopsis_pectinata)I10)I7,(Leucanthemopsis_alpina4x,(Leucanthemopsis_alpina,Leucanthemopsis_alpina_subsp_tomentosa)I9)I6)I4)I2)I1)I0;"
        );

        Tree gt1 = Trees.readTree("((Castrilanthemum_debeauxii_IA2170-04A1:0.003325333847,Castrilanthemum_debeauxii_IA2170-04A2:0.003325333847):0.142353004,(((((((Leucanthemopsis_pulverulenta_LPS134-01A2:0.003626327751,Leucanthemopsis_pallida_var_virescens_LPS185:0.003626327751):0.005732667214,((Leucanthemopsis_longipectinata_LPS189A1:0.005427109399,(Leucanthemopsis_alpina_subsp_tomentosa_LPS181-03A1:0.004195710171,Leucanthemopsis_alpina_subsp_tomentosa_LPS182-01A1:0.004195710171):0.001231399228):0.0004250762879,(((Leucanthemopsis_alpina_subsp_tomentosa_LPS181-03A2:0.0006163708318,Leucanthemopsis_longipectinata_LPS189A21:0.0006163708318):0.00282209522,((Leucanthemopsis_alpina_LPS064-01A1:0.0001322033992,Leucanthemopsis_alpina_LPS074-01:0.0001322033992):0.0004577158846,(Leucanthemopsis_alpina4x_LPS119-07:0.0001834718317,Leucanthemopsis_alpina_subsp_tomentosa_LPS182-01A2:0.0001834718317):0.0004064474521):0.002848546768):0.001983275192,(Leucanthemopsis_alpina_LPS064-01A2:0.002555193731,Leucanthemopsis_pulverulenta_LPS137-01A2:0.002555193731):0.002866547513):0.0004304444425):0.003506809277):0.004976661308,(((((Leucanthemopsis_pulverulenta_LPS142-02A2:0.0001967565991,Leucanthemopsis_pallida_var_bilbilitanum_LPS138-01A2:0.0001967565991):0.0007749569925,(Leucanthemopsis_pallida_var_bilbilitanum_LPS138-01A1:0.0002883081808,Leucanthemopsis_pulverulenta_LPS142-02A1:0.0002883081808):0.0006834054108):0.008987272541,Leucanthemopsis_pallida_subsp_spathulifolia_LPS150-01A5:0.009958986133):0.001303451602,((Leucanthemopsis_alpina_subsp_cuneata_LPS168-09A5:0.004182603449,(Leucanthemopsis_pallida_subsp_spathulifolia_LPS150-01A1:0.001608517103,Leucanthemopsis_pulverulenta_LPS137-01A1:0.001608517103):0.002574086347):0.004276193868,Leucanthemopsis_pulverulenta_LPS134-01A1:0.008458797317):0.002803640418):0.002212801875,(Leucanthemopsis_pallida_var_virescens_LPS163-02A2:0.004848982407,((Leucanthemopsis_pectinata_LPS166-11A1:0.0005958001694,Leucanthemopsis_pectinata_LPS167-01A2:0.0005958001694):0.003985571633,Leucanthemopsis_pallida_var_virescens_LPS163-02A1:0.004581371802):0.0002676106044):0.008626257203):0.0008604166624):0.003902336802,(Leucanthemopsis_pulverulenta_LPS159-01:0.007168022582,Leucanthemopsis_alpina_subsp_cuneata_LPS168-09A6:0.007168022582):0.01106997049):0.01303679816,(Leucanthemopsis_pectinata_LPS166-11A2:0.002353634317,Leucanthemopsis_pectinata_LPS167-01A1:0.002353634317):0.02892115692):0.01939852134,(Hymenostemma_pseudanthemis_LPS131-09A2:0.01132131926,Hymenostemma_pseudanthemis_LPS131-09A1:0.01132131926):0.03935199331):0.009563715348,(Prolongoa_hispanica_LPS135-10A1:0.005158182569,Prolongoa_hispanica_LPS135-10A2:0.005158182569):0.05507884535):0.08544130993);");
        Tree gt2 = Trees.readTree("((Castrilanthemum_debeauxii_IA2170-04A11:0.04866272479,Leucanthemopsis_longipectinata_LPS189A21:0.04866272479):0.00425007537,(((((Leucanthemopsis_alpina_LPS064-01A1:0.002606593801,Leucanthemopsis_longipectinata_LPS189A11:0.002606593801):0.0006122252751,(Leucanthemopsis_pulverulenta_LPS159-01A1:0.001817342664,Leucanthemopsis_alpina_subsp_cuneata_LPS168-09A4:0.001817342664):0.001401476412):0.009330555655,(Prolongoa_hispanica_LPS135-10A2:0.01253677216,Hymenostemma_pseudanthemis_LPS131-09:0.01253677216):1.260256857e-05):0.009829912582,(Prolongoa_hispanica_LPS135-10A1:0.01343877391,(((Leucanthemopsis_pectinata_LPS167-01A2:0.00146393232,Leucanthemopsis_pallida_subsp_spathulifolia_LPS150-01A1:0.00146393232):0.005635175451,((Leucanthemopsis_pectinata_LPS166-11A1:0.003132441282,(((Leucanthemopsis_pulverulenta_LPS134-01A:0.0004834146955,(Leucanthemopsis_pallida_var_virescens_LPS163-02:0.0002166129901,Leucanthemopsis_alpina_subsp_cuneata_LPS168-09A1:0.0002166129901):0.0002668017054):0.0008088488071,(Leucanthemopsis_pallida_var_virescens_LPS185:0.0002106351844,(Leucanthemopsis_pectinata_LPS167-01A1:0.0001657485932,Leucanthemopsis_pallida_var_bilbilitanum_LPS138-01:0.0001657485932):4.488659113e-05):0.001081628318):0.0007017906379,Leucanthemopsis_pectinata_LPS166-11A2:0.001994054141):0.001138387141):0.003085603239,((((Leucanthemopsis_pulverulenta_LPS142-02A1:0.0008830504479,Leucanthemopsis_pulverulenta_LPS159-01A2:0.0008830504479):3.439938081e-05,Leucanthemopsis_pallida_subsp_spathulifolia_LPS150-01A2:0.0009174498287):0.004223963518,Leucanthemopsis_alpina4x_LPS119-07A4:0.005141413347):0.0008928160864,((Leucanthemopsis_pulverulenta_LPS137-01A1:0.0004971846145,Leucanthemopsis_pulverulenta_LPS142-02A2:0.0004971846145):0.0003906973394,Leucanthemopsis_pulverulenta_LPS137-01A2:0.0008878819539):0.00514634748):0.0001838150875):0.0008810632498):0.0009475983809,((Leucanthemopsis_alpina_subsp_tomentosa_LPS181-03:0.004085814047,Leucanthemopsis_alpina_subsp_tomentosa_LPS182-01A1:0.004085814047):0.002415237713,((Leucanthemopsis_alpina_LPS064-01A2:0.0003130142879,Leucanthemopsis_alpina_LPS074-01:0.0003130142879):0.002414345952,(Leucanthemopsis_alpina4x_LPS119-07A1:0.000248527692,Leucanthemopsis_alpina_subsp_tomentosa_LPS182-01A2:0.000248527692):0.002478832548):0.00377369152):0.001545654392):0.005392067758):0.008940513403):0.01507757003,Castrilanthemum_debeauxii_IA2170-04A2:0.03745685734):0.01545594282);");
        Tree gt3 = Trees.readTree("(((Leucanthemopsis_longipectinata_LPS189A12:0.02942384723,(((Hymenostemma_pseudanthemis_LPS131-09A1:0.005121227087,(Leucanthemopsis_alpina4x_LPS119-07A1:0.004249140785,Leucanthemopsis_pectinata_LPS167-01A2:0.004249140785):0.0008720863018):0.005910869764,(Leucanthemopsis_pectinata_LPS167-01A1:0.002711380859,Leucanthemopsis_alpina_subsp_tomentosa_LPS181-03A1:0.002711380859):0.008320715991):0.006587357294,(Leucanthemopsis_pallida_var_bilbilitanum_LPS138-01A4:0.01188722887,Leucanthemopsis_pulverulenta_LPS159-01A2:0.01188722887):0.005732225276):0.01180439308):0.007193756791,((Castrilanthemum_debeauxii_IA2170-04:0.02757442249,(((Leucanthemopsis_pulverulenta_LPS137-01A2:0.009781590877,(Leucanthemopsis_pallida_subsp_spathulifolia_LPS150-01A1:0.009460392143,Leucanthemopsis_pulverulenta_LPS142-02A1:0.009460392143):0.0003211987335):0.005435965569,(Leucanthemopsis_pallida_var_virescens_LPS163-02A3:0.004419019272,Leucanthemopsis_pulverulenta_LPS134-01A1:0.004419019272):0.01079853717):0.009990145798,(Leucanthemopsis_pulverulenta_LPS137-01A1:0.01332397888,((Leucanthemopsis_pectinata_LPS166-11A2:0.0008822368353,Leucanthemopsis_pectinata_LPS166-11A1:0.0008822368353):0.00532231554,Leucanthemopsis_pallida_subsp_spathulifolia_LPS150-01A6:0.006204552375):0.007119426504):0.01188372337):0.002366720247):0.003790501511,(((Leucanthemopsis_longipectinata_LPS189A42:0.002051865123,(Leucanthemopsis_pulverulenta_LPS142-02A2:0.001831061895,(Leucanthemopsis_alpina_LPS064-01:0.001070863298,(Leucanthemopsis_pulverulenta_LPS134-01A2:0.0001943217482,Leucanthemopsis_alpina_subsp_tomentosa_LPS182-01A4:0.0001943217482):0.00087654155):0.0007601985972):0.0002208032278):0.001383699673,Leucanthemopsis_pulverulenta_LPS159-01A1:0.003435564796):0.01482616929,(Leucanthemopsis_alpina_subsp_cuneata_LPS168-09A5:0.008509598898,(((Leucanthemopsis_alpina_LPS074-01:0.0004763460806,Leucanthemopsis_alpina_subsp_tomentosa_LPS181-03A2:0.0004763460806):0.002877132792,Leucanthemopsis_alpina4x_LPS119-07A2:0.003353478873):0.002460211498,Leucanthemopsis_alpina_subsp_tomentosa_LPS182-01A1:0.005813690371):0.002695908527):0.009752135184):0.01310318992):0.005252680018):0.0007390005173,((Leucanthemopsis_pallida_var_bilbilitanum_LPS138-01A11:0.01578965214,((Hymenostemma_pseudanthemis_LPS131-09A2:0.01395213652,(Leucanthemopsis_pallida_var_virescens_LPS185:0.0002264593978,Leucanthemopsis_pallida_var_virescens_LPS163-02A1:0.0002264593978):0.01372567713):0.001353511553,Leucanthemopsis_alpina_subsp_cuneata_LPS168-09A11:0.01530564808):0.0004840040661):0.01984254592,Prolongoa_hispanica_LPS135-10:0.03563219807):0.001724406471);");
        Tree gt4 = Trees.readTree("((((((Leucanthemopsis_pulverulenta_LPS137-01:0.0009759867208,(Leucanthemopsis_pulverulenta_LPS159-01:0.0007202646041,Leucanthemopsis_pulverulenta_LPS134:0.0007202646041):0.0002557221167):0.003585236511,(Leucanthemopsis_pulverulenta_LPS142-02A4:0.002620176499,Leucanthemopsis_pulverulenta_LPS142-02A6:0.002620176499):0.001941046733):0.009015354046,(Leucanthemopsis_alpina_subsp_cuneata_LPS168-09A4:0.01314511895,((((Leucanthemopsis_alpina_subsp_tomentosa_LPS182-01:0.001066157523,Leucanthemopsis_alpina_subsp_tomentosa_LPS181:0.001066157523):0.003253584059,Leucanthemopsis_pallida_subsp_spathulifolia_LPS150-01:0.004319741582):0.002885653407,(Leucanthemopsis_pallida_var_virescens_LPS185_C4:0.005819664237,((Leucanthemopsis_pallida_var_virescens_LPS185_C6:0.0005506671305,Leucanthemopsis_pallida_var_bilbilitanum_LPS138:0.0005506671305):0.0009250292883,Leucanthemopsis_pallida_var_virescens_LPS163:0.001475696419):0.004343967818):0.001385730752):0.0003204652022,(Leucanthemopsis_pectinata_LPS166:0.001378467433,Leucanthemopsis_pectinata_LPS167-01:0.001378467433):0.006147392759):0.005619258761):0.0004314583247):0.02154447274,(Prolongoa_hispanica_LPS135:0.02816100267,(((Leucanthemopsis_alpina_LPS074:0.0004035051437,Leucanthemopsis_alpina4x_LPS119-07:0.0004035051437):0.001813891909,Leucanthemopsis_alpina_LPS064-01:0.002217397052):0.0006130343519,(Leucanthemopsis_alpina_subsp_cuneata_LPS168-09A3:0.0008099839403,Leucanthemopsis_longipectinata_LPS189:0.0008099839403):0.002020447464):0.02533057126):0.006960047348):0.005043748014,Hymenostemma_pseudanthemis_LPS131:0.04016479803):0.05685574807,(Castrilanthemum_debeauxii_IA2170_C8:0.004548153706,Castrilanthemum_debeauxii_IA2170_C2:0.004548153706):0.0924723924);");
        Tree gt5 = Trees.readTree("((Castrilanthemum_debeauxii_IA2170-trnC-petN:0.03821754715,((((Leucanthemopsis_alpina_LPS074-trnC-petN:0.002725686961,Leucanthemopsis_alpina_subsp_tomentosa_LPS181-trnC-petN:0.002725686961):0.003368944125,((Leucanthemopsis_alpina_LPS182-01-trnC-petN:0.000307793498,Leucanthemopsis_alpina4x_LPS119-07-trnC-petN:0.000307793498):0.0008922428686,Leucanthemopsis_alpina_LPS064-01-trnC-petN:0.001200036367):0.004894594719):0.005981311397,((((Leucanthemopsis_pulverulenta_LPS137-01-trnC-petN:0.001200523315,(Leucanthemopsis_pulverulenta_LPS142-02-trnC-petN:0.0004760222264,(Leucanthemopsis_pulverulenta_LPS134-trnC-petN:0.0003675459711,Leucanthemopsis_pulverulenta_LPS159-01-trnC-petN:0.0003675459711):0.0001084762553):0.0007245010885):0.001738370843,Leucanthemopsis_pallida_var_virescens_LPS163-trnC-petN:0.002938894158):0.000528057165,(Leucanthemopsis_pallida_var_bilbilitanum_LPS138-trnC-petN:0.002909545293,(Leucanthemopsis_pallida_subsp_spathulifolia_LPS150-01-trnC-petN:0.001644829845,((Leucanthemopsis_pallida_var_virescens_LPS185-trnC-petN:0.00110483341,Leucanthemopsis_alpina_subsp_cuneata_LPS168-09-trnC-petN:0.00110483341):0.0001619728811,(Leucanthemopsis_pectinata_LPS166-trnC-petN:0.0005518747866,Leucanthemopsis_pectinata_LPS167-01-trnC-petN:0.0005518747866):0.0007149315049):0.0003780235533):0.001264715449):0.00055740603):0.004452858173,Leucanthemopsis_longipectinata_LPS189-trnC-petN:0.007919809496):0.004156132986):0.008785136533,(Prolongoa_hispanica_LPS135-trnC-petN:0.01346883896,Hymenostemma_pseudanthemis_LPS131-trnC-petN:0.01346883896):0.007392240059):0.01735646813):0.5200943247,(((((Leucanthemopsis_alpina_LPS064-01-psbA-trnH:0.005385418899,Leucanthemopsis_alpina4x_LPS119-07-psbA-trnH:0.005385418899):0.001373554361,(Leucanthemopsis_alpina_LPS074-psbA-trnH:0.003310263725,(Leucanthemopsis_alpina_subsp_tomentosa_LPS181-psbA-trnH:0.002438705017,Leucanthemopsis_alpina_LPS182-01-psbA-trnH:0.002438705017):0.0008715587074):0.003448709535):0.007911082656,((((Leucanthemopsis_pectinata_LPS166-psbA-trnH:0.001916082047,Leucanthemopsis_pallida_var_virescens_LPS185-psbA-trnH:0.001916082047):0.00141607777,(Leucanthemopsis_alpina_subsp_cuneata_LPS168-09-psbA-trnH:0.0003979169398,Leucanthemopsis_pectinata_LPS167-01-psbA-trnH:0.0003979169398):0.002934242877):0.0001217030979,(((Leucanthemopsis_pulverulenta_LPS142-02-psbA-trnH:0.0001377467017,Leucanthemopsis_pulverulenta_LPS134-psbA-trnH:0.0001377467017):0.0006848453096,Leucanthemopsis_pallida_var_virescens_LPS163-psbA-trnH:0.0008225920114):0.0008838782513,((Leucanthemopsis_pulverulenta_LPS159-01-psbA-trnH:5.967921798e-05,Leucanthemopsis_pulverulenta_LPS137-01-psbA-trnH:5.967921798e-05):0.00100138753,Leucanthemopsis_pallida_var_bilbilitanum_LPS138-psbA-trnH:0.001061066748):0.0006454035149):0.001747392652):0.005920620617,(Leucanthemopsis_pallida_subsp_spathulifolia_LPS150-01-psbA-trnH:0.006800185404,Leucanthemopsis_longipectinata_LPS189-psbA-trnH:0.006800185404):0.002574298127):0.005295572385):0.02718771694,(Prolongoa_hispanica_LPS135-psbA-trnH:0.02351248166,Hymenostemma_pseudanthemis_LPS131-psbA-trnH:0.02351248166):0.0183452912):0.01199426788,Castrilanthemum_debeauxii_IA2170-psbA-trnH:0.05385204074):0.5044598311);");
        List<MutableTuple<Tree,Double>> geneTrees = new ArrayList<>();
        geneTrees.add(new MutableTuple<>(gt1, 0.0));
        geneTrees.add(new MutableTuple<>(gt2, 0.0));
        geneTrees.add(new MutableTuple<>(gt3, 0.0));
        geneTrees.add(new MutableTuple<>(gt4, 0.0));
        geneTrees.add(new MutableTuple<>(gt5, 0.0));
        List<String> speciesNames = Lists.newArrayList("Leucanthemopsis_alpina_subsp_tomentosa",
                "Leucanthemopsis_alpina_subsp_cuneata",
                "Leucanthemopsis_alpina4x",
                "Leucanthemopsis_alpina",
                "Leucanthemopsis_pectinata",
                "Leucanthemopsis_pallida_var_virescens",
                "Leucanthemopsis_pulverulenta",
                "Leucanthemopsis_pallida_var_bilbilitanum",
                "Leucanthemopsis_pallida_subsp_spathulifolia",
                "Leucanthemopsis_longipectinata",
                "Prolongoa_hispanica",
                "Hymenostemma_pseudanthemis",
                "Castrilanthemum_debeauxii");

        Map<String, String> allele2species = new HashMap<>();
        for (MutableTuple<Tree, Double> tuple : geneTrees) {
            for (String allele : tuple.Item1.getLeaves()) {
                for (String species : speciesNames) {
                    if (allele.startsWith(species)) {
                        allele2species.put(allele, species);
                        break;
                    }
                }
            }
        }
        Multimap<String, String> species2alleles = HashMultimap.create();
        for (String allele : allele2species.keySet()) {
            species2alleles.put(allele2species.get(allele), allele);
        }
        StringBuilder sb = new StringBuilder();
        sb.append("<");
        for (String s : species2alleles.keySet()) {
            sb.append(s + ":");
            boolean first = true;
            for (String a : species2alleles.get(s)) {
                if (first) {
                    first = false;
                    sb.append(a);
                } else {
                    sb.append(",").append(a);
                }

            }
            sb.append(";");
        }
        sb.append(">");
        System.out.println(sb.toString());

        for (String s : spNewicks) {
            Network speciesNetwork = Networks.readNetwork(s);
            for (Object obj : speciesNetwork.dfs()) {
                NetNode node = (NetNode) obj;
                int indeg = node.getIndeg();
                double prob = 1.0 / indeg;
                for (Object par : node.getParents()) {
                    if (node.getParentProbability((NetNode)par) == NetNode.NO_PROBABILITY)
                        node.setParentProbability((NetNode) par, prob);
                }

            }
            int xls[] = new int[geneTrees.size()];
            if (useVanillaMdc) {
                MDCOnNetwork mdc = new MDCOnNetwork();
                int index = 0;
                for (int xl : mdc.countExtraCoal(speciesNetwork, geneTrees, allele2species)) {
                    xls[index++] = xl;
                }
            } else {
                MDCOnAllopolyploidNetwork mdc = new MDCOnAllopolyploidNetwork();
                mdc.countExtraCoal(speciesNetwork, geneTrees, allele2species, xls);
            }



            System.out.print(Arrays.toString(xls) + ": " + Arrays.stream(xls).sum());
            System.out.println();

        }

    }

    private static void testCountExtraCoal3() {
        Tree geneTree = Trees.readTree("(((a,b),(x1,y1)),(c,((x2,z2),d)));");
        List<MutableTuple<Tree,Double>> geneTrees = new ArrayList<>();
        geneTrees.add(new MutableTuple<>(geneTree, 0.0));

        Network speciesNetwork = Networks.readNetwork("(((((X,Y),Z)#H1,B),A),((#H1,C),D));");
        Map<String, String> allele2species = new HashMap<>();
        allele2species.put("x1", "X");
        allele2species.put("x2", "X");
        allele2species.put("y1", "Y");
        allele2species.put("z2", "Z");
        allele2species.put("a", "A");
        allele2species.put("b", "B");
        allele2species.put("c", "C");
        allele2species.put("d", "D");
        MDCOnAllopolyploidNetwork mdc = new MDCOnAllopolyploidNetwork();
        int xls[] = new int[geneTrees.size()];
        mdc.countExtraCoal(speciesNetwork, geneTrees, allele2species, xls);
        System.out.println(Arrays.toString(xls));
    }

    private static void testCountExtraCoal2() {
        Tree geneTree = Trees.readTree("(((((x1,y1),z1),b),a),((((x2,y2),z2),c),d));");
        List<MutableTuple<Tree,Double>> geneTrees = new ArrayList<>();
        geneTrees.add(new MutableTuple<>(geneTree, 0.0));

        Network speciesNetwork = Networks.readNetwork("(((((X,Y),Z)#H1,B),A),((#H1,C),D));");
        Map<String, String> allele2species = new HashMap<>();
        allele2species.put("x1", "X");
        allele2species.put("x2", "X");
        allele2species.put("y1", "Y");
        allele2species.put("y2", "Y");
        allele2species.put("z1", "Z");
        allele2species.put("z2", "Z");
        allele2species.put("a", "A");
        allele2species.put("b", "B");
        allele2species.put("c", "C");
        allele2species.put("d", "D");
        MDCOnAllopolyploidNetwork mdc = new MDCOnAllopolyploidNetwork();
        int xls[] = new int[geneTrees.size()];
        mdc.countExtraCoal(speciesNetwork, geneTrees, allele2species, xls);
        System.out.println(Arrays.toString(xls));
    }

    private static void countExtraCoalSimData() {
        List<String> spNewicks = Lists.newArrayList(
                "(O,(C,((B,(T)#H1),(A,#H1))));", //scenario 1
                "(O,((B,(A,(T)#H1)),(C,#H1)));",         //scenario 2
                "(O,(((B,A),(T)#H1),(#H1,C)));"          //scenario 3
        );

    }



    public static void main(String[] args) {
    }
}
