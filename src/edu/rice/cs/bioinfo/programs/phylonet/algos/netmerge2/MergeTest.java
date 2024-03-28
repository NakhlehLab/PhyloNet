package edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge2;


import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

public class MergeTest {

    private final ByteArrayOutputStream outContent = new ByteArrayOutputStream();
    private final ByteArrayOutputStream errContent = new ByteArrayOutputStream();
    private final PrintStream originalOut = System.out;
    private final PrintStream originalErr = System.err;

    @Test
    public void mergeNetsViaNJ() {
    }

    @Before
//    public void mapSubNets(){
//        Network net1 = Networks.readNetwork("((F:3.0,((D,E))#H1:1.0):2.0,((C,(B,A)),#H1:2.0):1.0);");
//        Network net2 = Networks.readNetwork("((F,(E,D)),(((A,B))#H1:2.0,(C:2.0,#H1:1.0):1.0):1.0);");
//
//        Map<NetNode, NetNode> nodemap = Merge.mapSubNets(net1, net2);
//        for(NetNode n1: nodemap.keySet()){
//            try{
//                System.setOut(new PrintStream(n1+"--"+nodemap.get(n1)));
//            }catch (Exception e){
//
//            }
//
//        }
//
//    }
    @After
    public void restoreStreams() {
        System.setOut(originalOut);
        System.setErr(originalErr);
    }
}
