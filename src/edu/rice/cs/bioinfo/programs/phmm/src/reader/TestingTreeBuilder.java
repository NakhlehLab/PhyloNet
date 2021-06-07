/**
 * This file is part of PhyloNet-HMM.
 *
 * Copyright © 2013-2014 Kevin Liu, Jingxuan Dai, Kathy Truong, 
 * Ying Song, Michael H. Kohn, and Luay Nakhleh. <http://bioinfo.cs.rice.edu/>
 * 
 * PhyloNet-HMM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * PhyloNet-HMM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.rice.cs.bioinfo.programs.phmm.src.reader;

import java.util.ArrayList;

import edu.rice.cs.bioinfo.programs.phmm.src.phylogeny.EvoTree;
import edu.rice.cs.bioinfo.programs.phmm.src.phylogeny.Node;

public class TestingTreeBuilder {

    /**
     * @param args
     */
    public static void main(String[] args) {
//		HashMap<Integer, String> seqTypes = new HashMap<Integer, String>();
//		seqTypes.put(0, "human");
//		seqTypes.put(1, "chimp");
//		seqTypes.put(2, "gorilla");
//		seqTypes.put(3, "mouse");
//
//		AllTrees genTrees = new AllTrees(seqTypes);
//		ArrayList<EvoTree> alltrees = genTrees.getTrees();
//		for (int i = 0; i < alltrees.size(); i++) {
//			System.out.println(alltrees.get(i));
//		}
//
//		ArrayList<LeafNode> leaves = alltrees.get(0).getLeaves();
//		for (int i = 0; i < seqTypes.size(); i++) {
//			System.out.println(leaves.get(i));
//		}

        //make a tree
        EvoTree mytree = new EvoTree();
        Node root = new Node();


        Node in1 = new Node();
        Node leaf1 = new Node("human", in1, 10);
        Node leaf2 = new Node("chimp", in1, 5);
        Node leaf3 = new Node("gorilla", root, 3);
        // leaf1.setObs("A");
        // leaf2.setObs("T");
        // leaf3.setObs("G");
        in1.setParent(root);

        in1.addChild(leaf2);

        in1.addChild(leaf3);

        in1.setTbranch(20);

        root.addChildren(leaf1, in1);

        mytree.setRoot(root);
        mytree.setID(1);


        System.out.println("hey this is my tree : " + mytree);


        ArrayList<Node> leaves = mytree.getLeaves();
        for (int i = 0; i < leaves.size(); i++) {
            System.out.println(leaves.get(i));
        }

        //System.out.println("my likelihood is = " + mytree.getLikelihood());


    }


}
