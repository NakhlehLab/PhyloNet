package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import org.jblas.util.Random;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 9/1/18
 * Time: 5:28 PM
 * To change this template use File | Settings | File Templates.
 */
public class Test {
    static void testSamplingEmbedding() {
        double alpha1 = 0.5;
        double alpha2 = 0.4;
        double alpha3 = 0.6;
        double alpha4 = 0.8;
        double a = 0;
        double b = 0;
        double c = 0;
        double theta = 0.1;
        double gamma = 0.5;

        String netstring = String.format("[%f](((B:%f)I3#H1:%f::%f,A:%f)I1:%f,(C:%f,I3#H1:%f::%f)I2:%f)I0;", theta, alpha2, alpha1-alpha2, gamma, alpha1, alpha4-alpha1, alpha3, alpha3-alpha2, 1.0-gamma, alpha4-alpha3);
        //String netstring = "[0.02](((((C:0.01)#H1:0.02::0.6,B:0.03):0.01)#H2:0.02::0.7,A:0.06):0.03,((D:0.02,#H1:0.01::0.4):0.03,#H2:0.01::0.3):0.04);";
        //String netstring = "[0.03478971320729695](((((D:0.018773841591172784,C:0.018773841591172784):0.01665803647774957,B:0.035431878068922354):5.193091455444887E-4)#H1:0.042646185645307325::0.490471497495645,((#H1:0.01768017414602177::0.509528502504355,A:0.053631361360488614):2.981418077369069E-5)#H2:0.024936197318511863::0.28337456507014525):2.9602647955248584E-4,#H2:0.02523222379806435::0.7166254349298548);";
        //String netstring = "[0.07138571458425186]((((B:0.021545725769406308,(D:0.02072866724785186,C:0.02072866724785186):8.170585215544479E-4):0.012235537519323995)#H1:0.02015659627482603::0.17371727563554662,A:0.05393785956355633):0.005245366442141658,#H1:0.025401962716967687::0.8262827243644534);";
        Network<Object> trueNetwork;
        trueNetwork = Networks.readNetworkWithRootPop(netstring);
        System.out.println(Networks.getFullString(trueNetwork));

        String gtstring = "((A:1.0,B:1.0):0.3,C:1.3);";
        //String gtstring = "(((A:0.2,B:0.2):0.3,C:0.5):0.2,D:0.7);";
        //String gtstring = "((B:0.07,A:0.07)I1:0.03,(D:0.08,C:0.08)I2:0.02)I0;";
        String gtstring1 = "((B:0.12,A:0.12)I1:0.03,(D:0.13,C:0.13)I2:0.02)I0;";
        Tree gt0 = Trees.readTree(gtstring);
        List<UltrametricTree> gts = new ArrayList<>();
        for(int i = 0 ; i < 2 ; i++)
            gts.add(new UltrametricTree(Trees.readTree(gtstring)));
        //gts.add(new UltrametricTree(gt0));
        //gts.add(new UltrametricTree(gt0));

        //Utils._NUM_THREADS = 1;
        Utils._SEED = Random.nextInt(99999999);
        Utils._ESTIMATE_POP_SIZE = false;
        Utils._CONST_POP_SIZE = true;
        Utils.SAMPLE_EMBEDDINGS = true;
        UltrametricNetwork network = new UltrametricNetwork(trueNetwork, gts);
        network.setDirty(true);
        double prevLikelihood = network.logDensity();
        network._logGeneTreeNetwork = network.logDensity1();
        int acc = 0;

        if(Utils.SAMPLE_EMBEDDINGS) {
            Map<Tuple<TreeEmbedding, TreeEmbedding>, Double> likelihood2 = new HashMap<>();
            Map<Tuple<TreeEmbedding, TreeEmbedding>, Integer> embeddingCount2 = new HashMap<>();

            Map<TreeEmbedding, Double> likelihood = new HashMap<>();
            Map<TreeEmbedding, Integer> embeddingCount = new HashMap<>();
            Map<TreeEmbedding, TreeEmbedding> embeddings = new HashMap<>();
            for (int i = 0; i < 10000; i++) {
                network.setDirty(true);
                double embeddingLogHr = network.rebuildEmbeddings();
                //System.out.println(embeddingLogHr);
                double nextLikelihood = network.logDensity();
                //System.out.println(embeddingLogHr);
                double logAlpha = nextLikelihood - prevLikelihood + embeddingLogHr;
                System.out.println(embeddingLogHr + " " + nextLikelihood + " " + logAlpha);

                if (logAlpha >= Math.log(Randomizer.getRandomDouble())) {
                    network.accept();
                    acc++;
                    prevLikelihood = nextLikelihood;
                } else {
                    network.undo();
                    network.reject();
                }

                TreeEmbedding embedding = network._embeddings.get(0).clone();
                if (!embeddingCount.containsKey(embedding)) {
                    embeddingCount.put(embedding, 0);
                }
                embeddingCount.put(embedding, embeddingCount.get(embedding) + 1);
                likelihood.put(embedding, network._logGeneTreeNetwork[0]);

                if(gts.size() >= 2) {
                    Tuple tuple = new Tuple<>(network._embeddings.get(0).clone(), network._embeddings.get(1).clone());
                    likelihood2.put(tuple, network._logGeneTreeNetwork[0] + network._logGeneTreeNetwork[1]);
                    if (!embeddingCount2.containsKey(tuple)) {
                        embeddingCount2.put(tuple, 0);
                    }
                    embeddingCount2.put(tuple, embeddingCount2.get(tuple) + 1);
                }

                if(embeddings.containsKey(embedding) && embedding.probsum != embeddings.get(embedding).probsum) {
                    throw new RuntimeException("!!!");
                }
                if(embeddings.containsKey(embedding) && embedding.prob != embeddings.get(embedding).prob) {
                    throw new RuntimeException("!!!");
                }
                embeddings.put(embedding, embedding);

                // System.out.println("!!!!! " + network.logDensity());
            }

            double sum = 0;
            for (TreeEmbedding embedding : embeddingCount.keySet()) {
                System.out.println(embeddingCount.get(embedding) + " " + likelihood.get(embedding) + " " + embeddings.get(embedding).prob + " " + embeddings.get(embedding).probsum);
                sum += Math.exp(likelihood.get(embedding));
            }

            System.out.println("Embedding2");
            for (Tuple tuple : embeddingCount2.keySet()) {
                System.out.println(embeddingCount2.get(tuple) + " " + likelihood2.get(tuple));
            }
            System.out.println("Embedding2");

            System.out.println(Math.log(sum));
            System.out.println("prev likelihood " + prevLikelihood);
            System.out.println("Acc " + acc);
        } else {
            System.out.println(prevLikelihood);
        }
    }

    public static void main(String[] args) {
        testSamplingEmbedding();

    }
}
