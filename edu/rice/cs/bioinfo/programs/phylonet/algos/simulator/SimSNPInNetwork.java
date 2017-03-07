package edu.rice.cs.bioinfo.programs.phylonet.algos.simulator;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.distribution.SNAPPLikelihood;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.RPattern;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.SNAPPAlgorithm;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.OneNucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.*;
import java.util.concurrent.TimeUnit;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 11/27/16
 * Time: 11:50 AM
 * To change this template use File | Settings | File Templates.
 */
public class SimSNPInNetwork {
    public boolean _diploid = true;
    private Long _seed = null;
    private BiAllelicGTR _model;
    private SimGTInNetworkWithTheta _gtsim;
    private SimSNPInGT _snpsim;

    public SimSNPInNetwork(BiAllelicGTR model, Long seed) {
        _model = model;
        _seed = seed;
        _gtsim = new SimGTInNetworkWithTheta();
        _gtsim.setSeed(_seed);
        _snpsim = new SimSNPInGT(model);
        _snpsim.setSeed(_seed);
    }

    public Map<String, String> generateSNPs(Network network, Map<String, List<String>> species2alleles, int numGTs, boolean allowConstant){

        BiAllelicGTR model = _model;
        double pi0 = model.getEquilibriumVector().get(0, 0);
        double pi1 = model.getEquilibriumVector().get(1, 0);
        double u = 1.0 / (2.0 * pi0);
        double v = 1.0 / (2.0 * pi1);
        double mu = 2.0 * u * v / (u + v);
        Map<String, StringBuilder> res = new HashMap<>();

        if(species2alleles==null){
            species2alleles = new HashMap<String, List<String>>();
            for(Object leaf: network.getLeaves()){
                String species = ((NetNode)leaf).getName();
                List<String> alleles = new ArrayList<String>();
                alleles.add(species);
                species2alleles.put(species, alleles);
            }
        }

        if(_diploid) {
            for(String species : species2alleles.keySet()) {
                List<String> newalleles = new ArrayList<>();
                for(String allele : species2alleles.get(species)) {
                    newalleles.add(allele + "_0");
                    newalleles.add(allele + "_1");
                }
                species2alleles.put(species, newalleles);
            }
        }

        for(int i = 0 ; i < numGTs ; i++) {
            List<Tree> gts = _gtsim.generateGTs(network, species2alleles, mu, 1);
            Tree gt = gts.get(0);
            Map<String, String> snp = _snpsim.generateSingeSite(gt);

            boolean all1 = true;
            boolean all0 = true;
            for (String name : snp.keySet()) {
                if(snp.get(name).equals("0")) all1 = false;
                if(snp.get(name).equals("1")) all0 = false;
            }

            if(!allowConstant && (all1 || all0)) {
                i--;
                continue;
            }

            Map<String, Character> actualSite = new HashMap<>();
            for(String name : snp.keySet()) {
                if(_diploid) {
                    String actualName = name.substring(0, name.length() - 2);
                    if(!actualSite.containsKey(actualName))
                        actualSite.put(actualName, '0');

                    if(snp.get(name).equals("1")) {
                        actualSite.put(actualName, (char)(actualSite.get(actualName).charValue() + 1));
                    }
                } else {
                    actualSite.put(name, snp.get(name).charAt(0));
                }
            }

            for (String name : actualSite.keySet()) {
                if (!res.containsKey(name))
                    res.put(name, new StringBuilder());
                res.get(name).append(actualSite.get(name));
                //res.put(name, res.get(name) + snp.get(name));
            }

        }

        Map<String, String> ret = new HashMap<>();
        for(String s : res.keySet()) {
            ret.put(s, res.get(s).toString());
        }
        return ret;

    }

    public static void testDiploid(){
        double pi0 = 0.9;
        double pi1 = 1- pi0;

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {pi0, pi1}, new double[] {1.0/ (2.0 * pi0)});
        boolean useOnlyPolymorphic = false;
        int numSites = 1000000;

        double gamma = 1 - 0.3;

        double y = 1.5;
        double x = 0.5;

        Network<Object> trueNetwork;
        trueNetwork = Networks.readNetwork("((((b:1.0,c:1.0)I4:" + y + ")I3#H1:0.1::" + (1 - gamma) + ",a:" + (1.1 + y) + ")I1:" + x + ",(I3#H1:0.1::" + gamma + ",d:" + (1.1 + y) + ")I2:"+ x +")I0;");
        //trueNetwork = Networks.readNetwork("(((((A:0.7)I6#H1:1.3::0.8,Q:2.0)I4:1.0,L:3.0)I3:1.0,R:4.0)I2:1.0,(G:2.0,(I6#H1:0.7::0.2,C:1.4)I5:0.6)I1:3.0)I0;");
        //trueNetwork = Networks.readNetwork("(((((Q:2.0,A:2.0)I4:1.0,L:3.0)I3:0.5)I8#H1:0.5::0.7,R:4.0)I2:1.0,(I8#H1:0.5::0.3,(G:2.0,C:2.0)I1:2.0)I7:1.0)I0;");
        //trueNetwork = Networks.readNetwork("(((((Q:0.5)I8#H1:0.5::0.7,A:1.0)I4:1.0,L:2.0)I3:2.0,(I8#H1:1.0::0.3,R:1.5)I7:2.5)I2:1.0,(G:1.0,C:1.0)I1:4.0)I0;");
        trueNetwork = Networks.readNetwork("(((((Q:0.5)I8#H1:0.5::0.7,(A:0.5)I6#H2:0.5::0.8)I4:1.0,L:2.0)I3:2.0,(I8#H1:1.0::0.3,R:1.5)I7:2.5)I2:1.0,((I6#H2:0.5::0.2,C:1.0)I5:1.0,G:2.0)I1:3.0)I0;");

        double constTheta = 0.036;
        double mu = 1.0;
        for(Object nodeObject : Networks.postTraversal(trueNetwork)) {
            NetNode node = (NetNode) nodeObject;
            for(Object parentObject :  node.getParents()) {
                NetNode parent = (NetNode) parentObject;
                //node.setParentSupport(parent, constTheta);
                node.setParentDistance(parent, node.getParentDistance(parent) * constTheta / 2.0);
            }
        }
        trueNetwork.getRoot().setRootPopSize(constTheta);
//        trueNetwork = Networks.readNetwork("((C:0.09604283044864592:0.036,G:0.09604283044864592:0.036)II1:0.08821096143809048:0.036,(R:0.16796441222747166:0.036,(L:0.13918714276096708:0.036,(A:0.10079837440250675:0.036,Q:0.10079837440250675:0.036)II4:0.038388768358460335:0.036)II3:0.02877726946650458:0.036)II2:0.016289379659264747:0.036)II0;");
//        trueNetwork.getRoot().setRootPopSize(0.036);
//        trueNetwork = Networks.readNetwork("(((G:0.03271098551352612,(C:0.014581262376039935,(A:0.007848436317467452)#H2:0.006732826058572483::0.20092242303662247):0.01812972313748619):6.589767080962955E-4)#H1:0.056974310762192545::0.973670280875714,(((Q:0.01575520654534508)#H3:0.008437888036381732::0.29205403332422386,R:0.02419309458172681):0.04760724352481251,(#H1:0.012622246801820981::0.026329719124286055,((#H3:1.8671708033869347E-4::0.7079459666757761,#H2:0.00809348730821632::0.7990775769633776):0.01963878547064636,L:0.03558070909633013):0.010411499927113266):0.025808129083095925):0.01854393487727564);");
//        trueNetwork.getRoot().setRootPopSize(0.036);

        int nameCount = 0;
        for(Object nodeObject : Networks.postTraversal(trueNetwork)) {
            NetNode node = (NetNode) nodeObject;
            for(Object parentObject :  node.getParents()) {
                NetNode parent = (NetNode) parentObject;
                node.setParentSupport(parent, trueNetwork.getRoot().getRootPopSize());
            }
        }
        for(Object node : trueNetwork.dfs()) {
            NetNode mynode = (NetNode) node;
            if(mynode.getName().equals("")) {
                mynode.setName("I" + nameCount);
                nameCount++;
            }
        }

        SimSNPInNetwork simulator = new SimSNPInNetwork(BAGTRModel, null);
        Map<String, String> onesnp = simulator.generateSNPs(trueNetwork, null, numSites, !useOnlyPolymorphic);

        List<Alignment> alns = new ArrayList<>();
        Alignment aln = new Alignment(onesnp);
        alns.add(aln);
        aln._RPatterns = SNAPPLikelihood.diploidSequenceToPatterns(null, alns);

        int count = 0;

        for(RPattern pattern : aln._RPatterns.keySet()) {
            //if(count > 10) break;
            count++;
            Network cloneNetwork = Networks.readNetwork(trueNetwork.toString());
            cloneNetwork.getRoot().setRootPopSize(trueNetwork.getRoot().getRootPopSize());
            SNAPPAlgorithm run = new SNAPPAlgorithm(cloneNetwork, BAGTRModel, null);
            double likelihood = 0;
            try {
                long start = System.currentTimeMillis();
                likelihood = run.getProbability(pattern);
                System.out.println(pattern + " " + likelihood * Math.exp(aln._RPatterns.get(pattern)[1]) * numSites  + " " + aln._RPatterns.get(pattern)[0] );
                System.out.println("Time: " + (System.currentTimeMillis()-start)/1000.0);
            } catch(Exception e) {
                e.printStackTrace();
                System.out.println("Exceptional network");
            }
        }

        return;
    }

    public static void basic() {

        double pi0 = 0.9;
        double pi1 = 1- pi0;

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {pi0, pi1}, new double[] {1.0/ (2.0 * pi0)});
        boolean useOnlyPolymorphic = false;
        int numSites = 100000;

        double gamma = 1 - 0.3;

        double y = 1.5;
        double x = 0.5;

        Network<Object> trueNetwork;
        trueNetwork = Networks.readNetwork("((((b:1.0,c:1.0)I4:" + y + ")I3#H1:0.1::" + (1 - gamma) + ",a:" + (1.1 + y) + ")I1:" + x + ",(I3#H1:0.1::" + gamma + ",d:" + (1.1 + y) + ")I2:"+ x +")I0;");
        //trueNetwork = Networks.readNetwork("(((((A:0.7)I6#H1:1.3::0.8,Q:2.0)I4:1.0,L:3.0)I3:1.0,R:4.0)I2:1.0,(G:2.0,(I6#H1:0.7::0.2,C:1.4)I5:0.6)I1:3.0)I0;");
        //trueNetwork = Networks.readNetwork("(((((Q:2.0,A:2.0)I4:1.0,L:3.0)I3:0.5)I8#H1:0.5::0.7,R:4.0)I2:1.0,(I8#H1:0.5::0.3,(G:2.0,C:2.0)I1:2.0)I7:1.0)I0;");
        //trueNetwork = Networks.readNetwork("(((((Q:0.5)I8#H1:0.5::0.7,A:1.0)I4:1.0,L:2.0)I3:2.0,(I8#H1:1.0::0.3,R:1.5)I7:2.5)I2:1.0,(G:1.0,C:1.0)I1:4.0)I0;");
        //trueNetwork = Networks.readNetwork("(((((Q:0.5)I8#H1:0.5::0.7,(A:0.5)I6#H2:0.5::0.8)I4:1.0,L:2.0)I3:2.0,(I8#H1:1.0::0.3,R:1.5)I7:2.5)I2:1.0,((I6#H2:0.5::0.2,C:1.0)I5:1.0,G:2.0)I1:3.0)I0;");
        trueNetwork = Networks.readNetwork("((R:1.0)I10#H1:4.0::0.6,(((Q:2.0,A:2.0)I4:1.0,(L:2.0,I10#H1:1.0::0.4)I9:1.0)I3:1.0,(G:2.0,C:2.0)I1:2.0)I7:1.0)I0;");

        System.out.println(trueNetwork.toString());

        double constTheta = 0.036;
        double mu = 1.0;
        for(Object nodeObject : Networks.postTraversal(trueNetwork)) {
            NetNode node = (NetNode) nodeObject;
            for(Object parentObject :  node.getParents()) {
                NetNode parent = (NetNode) parentObject;
                //node.setParentSupport(parent, constTheta);
                node.setParentDistance(parent, node.getParentDistance(parent) * constTheta / 2.0);
            }
        }
        trueNetwork.getRoot().setRootPopSize(constTheta);


        //trueNetwork = Networks.readNetwork("((G:0.018380841366896713,C:0.018380841366896713):0.07259787745284114,(((Q:0.019223930643805762)#H1:0.0079382721121616::0.3027236780708331,R:0.02716220275596736):0.04612672867584792,((L:0.02760957365672652)#H2:0.01933595272319038::0.07225913233066261,((#H1:3.112852137000452E-4::0.6972763219291669,A:0.019535215857505808):0.017171505371586866,#H2:0.009097147572366154::0.9277408676693374):0.010238805150824225):0.02634340505189838):0.01768978738792258);");
        //trueNetwork.getRoot().setRootPopSize(0.03434332011291962);
        //trueNetwork = Networks.readNetwork("(((((R:0.03692297297059491)#H1:0.026335083580234356::0.32659382420676164,(A:0.048296112704057884,(Q:0.0016061601679680012)#H2:0.04668995253608988::0.8824921981923378):0.01496194384677138):4.711142085103842E-4,((L:0.03025867321849643)#H3:0.03261745957608034::0.675454288991292,((G:0.0514837139295224,C:0.0514837139295224):6.708696808103431E-4)#H4:0.010721549184244027::0.16076035943533917):8.530379647628816E-4):0.04434852129641825,#H4:0.05592310844542516::0.8392396405646608):9.399206670118038E-4,((#H3:0.003242697377805359::0.324545711008708,#H2:0.031895210428333785::0.11750780180766218):0.058420991394148414,#H1:0.054999389019855294::0.6734061757932384):0.017095250732319503);");
        //trueNetwork.getRoot().setRootPopSize(6.14445682191189E-4);
        //trueNetwork = Networks.readNetwork("(((G:0.040360738566830895,C:0.040360738566830895):0.002333190642093877)#H1:5.486184649373872::0.04647906647548805,(((R:0.08110800117484995,((L:0.015242059440912224)#H4:0.06507842315816675::0.4169585812227872,((Q:0.01964782667669728)#H2:0.01030261251091685::0.49982865039323054,(#H2:0.010276565369475205::0.5001713496067695,(A:8.140857204844255E-5)#H3:0.029842983474124046::0.44717984052179527):2.6047141441644384E-5):0.05037004341146485):7.875185757709735E-4):1.3144339721918308E-4,#H4:0.06599738513115691::0.5830414187772128):0.020102587039841108,(#H1:0.037216724597387355::0.953520933524512,#H3:0.07982924523426368::0.5528201594782047):0.021431377805598117):5.427536546970886);");
        //trueNetwork = Networks.readNetwork("(((G:0.1,C:0.1):0.002)#H1:5.5::0.5,(((R:0.08,(L:0.09,(Q:0.03,A:0.03):0.06):8E-4):1E-4):0.02,#H1:0.06::0.5):5);");
        //trueNetwork = Networks.readNetwork("(((C:0.5,G:0.5)I3:0.5)I2#H1:6.0::0.5,((((A:0.2,Q:0.2)I4:0.4,L:0.6)I5:4E-4,R:1.0)I6:5E-4,I2#H1:0.5::0.5)I1:5.5)I0;");
        //trueNetwork.getRoot().setRootPopSize(0.00817528575792145);
        //trueNetwork = Networks.readNetwork("((((C:0.006,(B:0.004)I7#H1:0.002::0.2)I4:0.004,(I7#H1:0.004::0.8)I6#H2:0.002::0.2)I3:0.004,(I6#H2:0.004::0.8)I5#H3:0.002::0.2)I2:0.006,(I5#H3:0.004::0.8,A:0.016)I1:0.004)I0;");
        //trueNetwork = Networks.readNetwork("((((C:0.07,(B:0.06)I7#H1:0.01::0.2)I4:0.02,(I7#H1:0.02::0.8)I6#H2:0.01::0.2)I3:0.02,(I6#H2:0.02::0.8)I5#H3:0.01::0.2)I2:0.07,(I5#H3:0.06::0.8,A:0.16)I1:0.02)I0;");
        //trueNetwork.getRoot().setRootPopSize(0.04);
        //trueNetwork = Networks.readNetwork("(((G:0.05211670360238044,C:0.05211670360238044):5.458297681763952E-4)#H1:0.05908373214397336::0.7705581304452586,(((#H1:1.5831098796992293E-5::0.22944186955474144,(L:0.007831161802081543)#H2:0.04484720266727228::0.7415692509081577):0.01812915575416195,(A:0.051817335095614515,(Q:0.005097739344889214)#H3:0.046719595750725304::0.8907488990731456):0.01899018512790126):0.02668594384866703,(R:0.08512519218544201,(#H3:0.007550537153160525::0.10925110092685442,#H2:0.004817114695968195::0.2584307490918423):0.07247691568739227):0.012368271886740792):0.01425280144234739);");
        //trueNetwork.getRoot().setRootPopSize(3.3418436760954384E-4);
        //trueNetwork = Networks.readNetwork("((((L:0.009304825035069735)#H1:14.370818397517256::0.020382388550110098)#H2:25.83351345488807::0.28805683793703096,((#H1:0.0013174924892173801::0.9796176114498899,(Q:0.009166552336431677,(A:0.005296210274212645,(G:0.0014429534014067207,C:0.0014429534014067207):0.003853256872805924):0.0038703420622190326):0.0014557651878554373):0.0047109417058845745,R:0.015333259230171689):40.19830341821022):75.28824647371039,#H2:101.12175992859846::0.711943162062969);");
        //trueNetwork.getRoot().setRootPopSize(0.011);
        //trueNetwork = Networks.readNetwork("((R:0.036413848502975396)#H1:0.05317572836552176::0.4832649910028073,((#H1:0.018089072459181774::0.5167350089971927,((Q:0.032272732571886384,A:0.032272732571886384):0.022057943031171354,L:0.05433067560305774):1.7224535909943217E-4):0.023327475186427166,(C:0.030256525309546957,G:0.030256525309546957):0.04757387083903738):0.011759180719912818);");
        //trueNetwork.getRoot().setRootPopSize(0.047915145174597676);
        //trueNetwork = Networks.readNetwork("((((((((C:6.071885919301598E-4,B:6.071885919301598E-4):0.004084977900907154)#H2:1.9180763597067933E-4::0.9984510547505954,A:0.004883974128807993):0.0461067287811829,#H2:0.04629853641715358::0.0015489452494046319):12.119610091552316)#H1:1.178512708785492::0.12210260828789488)#H3:0.0696553828071309::0.9838799395493458,#H1:1.248168091592623::0.8778973917121051):92.1852811291236,#H3:92.25493651193074::0.016120060450654172);");
        //trueNetwork.getRoot().setRootPopSize(9.057126525354496E-12);
        //trueNetwork = Networks.readNetwork("(C:920.8473410747589,(B:626.3380424743829,A:626.3380424743829):294.50929860037604);");
        //trueNetwork.getRoot().setRootPopSize(0.04);
        //trueNetwork = Networks.readNetwork("(((C:0.01659501739858635,((A:0.005800146445648392)#H2:5.156005619451585E-4::0.2744671496877781)#H1:0.0102792703909928::0.7233366864278274):0.017378455723520436,G:0.03397347312210679):0.052230584880487166,(((Q:0.006777148338408591)#H3:0.01823506491202251::0.2984279870355976,R:0.025012213250431098):0.04370832343492151,(L:0.03386930402729607,((#H3:0.007964099302371526::0.7015720129644024,#H2:0.008941101195131725::0.7255328503122219):0.018542674482596277,#H1:0.026968175115782842::0.27666331357217255):5.853819039196781E-4):0.034851232658056536):0.017483521317241346);");
        //trueNetwork.getRoot().setRootPopSize(0.03776047660800021);
        //trueNetwork = Networks.readNetwork("((C:0.09604283044864592:0.036,G:0.09604283044864592:0.036)II1:0.08821096143809048:0.036,(R:0.16796441222747166:0.036,(L:0.13918714276096708:0.036,(A:0.10079837440250675:0.036,Q:0.10079837440250675:0.036)II4:0.038388768358460335:0.036)II3:0.02877726946650458:0.036)II2:0.016289379659264747:0.036)II0;");
        //trueNetwork.getRoot().setRootPopSize(0.036);
        trueNetwork = Networks.readNetwork("((R:0.01)#H1:0.7::0.3,(L:0.13,((A:0.05,((C:0.01,G:0.01):0.01,#H1:0.01::0.7):0.03):0.01,Q:0.06):0.07):0.71);");
        trueNetwork.getRoot().setRootPopSize(0.01);
        //trueNetwork = Networks.readNetwork("((R:0.01)#H1:1.7::0.3,(L:0.13,((A:0.05,((C:0.01,G:0.01):0.01,#H1:0.01::0.7):0.03):0.01,Q:0.06):0.07):1.71);");
        //trueNetwork.getRoot().setRootPopSize(0.01);


        int nameCount = 0;
        for(Object nodeObject : Networks.postTraversal(trueNetwork)) {
            NetNode node = (NetNode) nodeObject;
            for(Object parentObject :  node.getParents()) {
                NetNode parent = (NetNode) parentObject;
                node.setParentSupport(parent, trueNetwork.getRoot().getRootPopSize());
            }
        }
        for(Object node : trueNetwork.dfs()) {
            NetNode mynode = (NetNode) node;
            if(mynode.getName().equals("")) {
                mynode.setName("I" + nameCount);
                nameCount++;
            }
        }


        SimSNPInNetwork simulator = new SimSNPInNetwork(BAGTRModel, null);
        simulator._diploid = false;
        Map<String, String> onesnp = simulator.generateSNPs(trueNetwork, null, numSites, !useOnlyPolymorphic);

        List<Alignment> alns = new ArrayList<>();
        Alignment aln = new Alignment(onesnp);
        alns.add(aln);
        aln._RPatterns = SNAPPLikelihood.haploidSequenceToPatterns(null, alns);

        int count = 0;

        for(RPattern pattern : aln._RPatterns.keySet()) {
            //if(count > 10) break;
            count++;
            Network cloneNetwork = Networks.readNetwork(trueNetwork.toString());
            cloneNetwork.getRoot().setRootPopSize(trueNetwork.getRoot().getRootPopSize());
            SNAPPAlgorithm run = new SNAPPAlgorithm(cloneNetwork, BAGTRModel, null);
            double likelihood = 0;
            try {
                long start = System.currentTimeMillis();
                likelihood = run.getProbability(pattern);
                System.out.println(pattern + " " + likelihood * Math.exp(aln._RPatterns.get(pattern)[1]) * numSites  + " " + aln._RPatterns.get(pattern)[0] );
                System.out.println("Time: " + (System.currentTimeMillis()-start)/1000.0);
            } catch(Exception e) {
                e.printStackTrace();
                System.out.println("Exceptional network");
            }
        }

//        List<Map<String, String>> snpdata = new ArrayList<>();
//        snpdata.add(onesnp);
//        List<Alignment> alns = new ArrayList<>();
//        for(Map<String, String> input : snpdata) {
//            for(String allele : input.keySet()) {
//                System.out.println(allele + " " + input.get(allele));
//            }
//            System.out.println();
//            Alignment aln = new Alignment(input);
//
//            Map<Integer, Integer> cache = new HashMap<>();
//
//            for(int i = 0 ; i < aln.getSiteCount() ; i++) {
//                Map<String, Character> colorMap = new TreeMap<String, Character>();
//                for(String taxon : aln.getAlignment().keySet()) {
//                    colorMap.put(taxon, aln.getAlignment().get(taxon).charAt(i));
//                }
//                Integer represent = 0;
//                for(String s : colorMap.keySet()) {
//                    represent = (represent << 1) + (colorMap.get(s) == '0' ? 0 : 1);
//                }
//                if(!cache.containsKey(represent)) {
//                    cache.put(represent, 0);
//                }
//                cache.put(represent, cache.get(represent) + 1);
//            }
//
//            if(useOnlyPolymorphic) {
//                cache.put(0, 0);
//                cache.put((1 << input.keySet().size()) - 1, 0);
//            }
//
//            aln.setCache(cache);
//            alns.add(aln);
//
//
//            for(Alignment alg : alns) {
//                double P0 = 0;
//                double P1 = 0;
//                Map<Integer, Double> likelihoods = new TreeMap<>();
//                List<String> names = alg.getTaxaNames();
//                double sum = 0;
//
//                for(Integer represent : alg.getCache().keySet()) {
//                    int cur = represent;
//                    Integer count = alg.getCache().get(cur);
//                    Map<String, Character> colorMap = new HashMap<String, Character>();
//                    for(int i = names.size() - 1 ; i >= 0 ; i--) {
//                        Character site = cur % 2 == 1 ? '1' : '0';
//                        colorMap.put(names.get(i), site);
//                        cur /= 2;
//                    }
//                    OneNucleotideObservation converter = new OneNucleotideObservation(colorMap);
//                    Network cloneNetwork = Networks.readNetwork(trueNetwork.toString());
//                    cloneNetwork.getRoot().setRootPopSize(trueNetwork.getRoot().getRootPopSize());
//                    SNAPPAlgorithm run = new SNAPPAlgorithm(cloneNetwork, BAGTRModel, null);
//                    double likelihood = 0;
//                    try {
//                        likelihood = run.getProbability(converter);
//                        if(represent == 0) P0 = likelihood;
//                        if(represent == ((1 << names.size()) - 1)) P1 = likelihood;
//                        sum += likelihood;
//                        likelihoods.put(represent, likelihood);
//                        //System.out.println(represent + " " + likelihood  + " " + count );
//                    } catch(Exception e) {
//                        e.printStackTrace();
//                        System.out.println("Exceptional network");
//                    }
//
//
//                }
//
//                for(Integer represent : likelihoods.keySet()) {
//                    Integer count = alg.getCache().get(represent);
//                    double likelihood = likelihoods.get(represent);
//                    if(useOnlyPolymorphic)  likelihood /= (1.0 - P0 - P1);
//                    System.out.println(represent + " " + likelihood  + " " + count );
//                }
//                System.out.println("sum " + sum);
//
//
//            }
//        }
    }

    public static void main(String []args) {
        //testDiploid();
        basic();
    }
}
