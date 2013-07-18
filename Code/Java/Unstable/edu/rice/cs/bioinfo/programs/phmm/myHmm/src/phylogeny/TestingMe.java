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
