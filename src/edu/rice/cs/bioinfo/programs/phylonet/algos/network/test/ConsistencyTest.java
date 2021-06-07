package edu.rice.cs.bioinfo.programs.phylonet.algos.network.test;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by hunter on 8/16/18.
 */
public class ConsistencyTest {
    public static void main(String[] args) {
        ConsistencyTest cTest = new ConsistencyTest();
        cTest.runTest();
    }

    public static void runTest() {
        TestIntegratedProbability test = new TestIntegratedProbability();
        Network netForTest = test.getScenario(2, 2);
        List<Network> netsForTest = new ArrayList<>();
        netsForTest.add(netForTest);
        for (Network n: netsForTest) {
            System.out.println("testing " + n);
            List<Tree> trees = getAllTrees(n);
        }
    }

    private static List<Tree> getAllTrees(Network n) {
        return null;
    }
}
