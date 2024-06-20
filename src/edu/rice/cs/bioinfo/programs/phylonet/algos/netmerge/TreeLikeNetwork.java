//package edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge;
///*
// * @ClassName:   TreeLikeNetwork
// * @Description:
// * @Author:      Zhen Cao
// * @Date:        3/24/24 9:13 PM
// */
//
//import javax.swing.*;
//import java.util.ArrayList;
//import java.util.LinkedList;
//
//public class TreeLikeNetwork {
//    /* Constructor */
//    public TreeLikeNetwork() {
//    }
//
//    public void run(){
//        long startTime = System.currentTimeMillis();
//        boolean mergeSubTrees=true;
//        try{
//            if(root==null){ return; } // Unrooted tree cannot be height sorted
//            if(root.getHeight()==-1){ heightLabel(); }
//            LinkedList[] heightLinkedList = new LinkedList[root.getHeight()+1]; // Not going to use height 0
//            root.computeIsomorphismCode(); // Will label all children recursively
//            root.computeClusterCode(); // Will label all children recursively
//
//            // Perform breadth first enumeration
//            heightLinkedList[root.getHeight()] = new LinkedList();
//            heightLinkedList[root.getHeight()].add(root);
//
//            for(int h=root.getHeight();h>=0;h--){
//                //		System.out.println("Analysing height " + h + " Num of subtrees " + heightLinkedList[h].size());
//                while(heightLinkedList[h].isEmpty()==false){
//                    Vertex currentVertex = (Vertex)heightLinkedList[h].removeFirst();
//                    Vertex otherVertex = null;
//
//                    // If it's isomorphic to any of the others then merge it
//                    ArrayList isomorphicTrees=null; // This represents trees T2...Tn
//
//                    for(int i=0;i<heightLinkedList[h].size();i++){
//                        otherVertex = (Vertex)heightLinkedList[h].get(i);
//                        if(currentVertex.equals(otherVertex)){
//                            // Check that the user wishes to merge these
//                            //			System.out.println("Isomorphic subtree found!");
//
//                            boolean merge=true;
//                            if(askBeforeMerging){
//                                currentVertex.setFlagWholeSubTree(true); otherVertex.setFlagWholeSubTree(true);
//                                if(thisFrame!=null){ thisFrame.repaint(); }
//
//                                final Object[] options = {"Yes, merge them","No, leave separate"};
//                                merge=(JOptionPane.showOptionDialog(thisFrame,"Would you like to merge the flagged pair of inextendible subtrees?","Isomorphic subtrees identified",JOptionPane.YES_NO_OPTION,JOptionPane.QUESTION_MESSAGE,null,options,options[0])==0);
//
//                                currentVertex.setFlagWholeSubTree(false); otherVertex.setFlagWholeSubTree(false);
//                            }
//
//                            if(merge){
//                                if(isomorphicTrees==null){
//                                    isomorphicTrees = new ArrayList();
//                                    isomorphicTrees.add(currentVertex);
//                                }
//                                isomorphicTrees.add(otherVertex);
//
//                                heightLinkedList[h].remove(i--);
//                            }
//                        }
//                    }
//
//                    if(mergeSubTrees){
//                        if(isomorphicTrees!=null){
//                            Vertex[] identifyVertices = new Vertex[isomorphicTrees.size()];
//                            Vertex youngestVertex=null; // The youngest parent ie max(distance_to_root)
//
//                            for(int i=0;i<isomorphicTrees.size();i++){
//                                Vertex currentIsomorphicVertex = (Vertex)isomorphicTrees.get(i);
//
//                                // Keep the vertex furthest ahead in time to ensure all parents are ancestors to their children
//                                if(youngestVertex==null){ youngestVertex=currentIsomorphicVertex; }
//                                else{
//                                    //	System.out.println("DistToRoot" + youngestVertex.getParent().getX());
//                                    if(youngestVertex.getParent().getX()<currentIsomorphicVertex.getParent().getX()){
//                                        youngestVertex=currentIsomorphicVertex;
//                                    }
//                                }
//
//                                // Subdivide
//                                identifyVertices[i]=currentIsomorphicVertex.getParentEdge().subdivide(thisNetwork);
//                            }
//
//                            // Add youngest parent to newVertex as we're keeping it
//                            currentVertex = youngestVertex;
//
//                            // for each identifyVertex which isn't youngest parent indentify with youngestParent
//                            for(int i=0;i<identifyVertices.length;i++){
//                                if(identifyVertices[i]!=youngestVertex.getParent()){
//                                    youngestVertex.getParent().identifyWith(identifyVertices[i]);
//                                    identifyVertices[i].destroy(thisNetwork); // Prune
//                                }
//                            }
//                            youngestVertex.getParent().ensureExistsAfterParents(thisNetwork);
//
//                            if(stepOver){ heightLabel(); depthLabel(); parent.repaint(); suspendThread(); }
//                            while(threadSuspended)
//                                sleep(100);
//                        }
//                    }
//
//
//                    // This subtree may still have sub-trees which are isomorphic to other sub-trees
//                    for(int i=0;i<currentVertex.getOutDegree();i++){
//                        Vertex vertex = (Vertex)currentVertex.getChildVertexAt(i);
//                        int height = vertex.getHeight();
//
//                        if(heightLinkedList[height]==null){ heightLinkedList[height] = new LinkedList(); }
//                        heightLinkedList[height].add(vertex);
//                    }
//
//                    convertToArray();
//                }
//            }
//
//            convertToArray();
//
//            // Because we subdivided edges and introduced new vertices the height labels will be wrong so need to height label again
//            heightLabel();
//
//            endThread();
//            parent.repaint();
//            // Give the user a little prompt to say no more found
//            if(thisFrame!=null){ JOptionPane.showMessageDialog(null, "No more inextendible subtrees found!", "Algorithm finished", JOptionPane.PLAIN_MESSAGE); }
//            //	System.out.println("Running time: " + (System.currentTimeMillis() - startTime));
//        }
//        catch(Exception e){ e.printStackTrace(); }
//    }
//
//    public static void main(String[] args) {
//
//    }
//
//
//}
