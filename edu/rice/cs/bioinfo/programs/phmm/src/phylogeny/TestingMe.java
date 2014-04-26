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

package phylogeny;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;

public class TestingMe {


    public static void main (String[] args) throws FileNotFoundException {

        BufferedReader r = new BufferedReader(new FileReader("treeoflife.tree"));
        TreeParser tp = new TreeParser(r);
        ArrayList<EvoTree> treeoflife = tp.nexusFileTreeNames("treeoflife.tree");
        for (int i = 0; i < treeoflife.size(); i ++) {
            System.out.println(((EvoTree)treeoflife.get(i)));
            //recursive_print(0, 0, (EvoTree) treeoflife.get(i));
            System.out.println("\n\n");
        }



        //treeoflife = tp.tokenize(1, "treeoflife", null);
        //int tree_height = treeoflife.getHeight();
        //System.out.println("largest tree height is: " + tree_height);

    }

//    static void recursive_print (int currkey, int currdepth, Tree atree) {
//        Node currNode = atree.getNodeByKey(currkey);
//        //int numChildren = currNode.numberChildren();
//        for (int i = 0; i < numChildren; i++) {
//            int childkey = currNode.getChild(i).key;
//            Node childnode = atree.getNodeByKey(childkey);
//            System.out.println("child name is: " + childnode.getName()
//                                 + " depth is: " + currdepth);
//            recursive_print(childkey, currdepth+1, atree);
//        }
//    }

}
