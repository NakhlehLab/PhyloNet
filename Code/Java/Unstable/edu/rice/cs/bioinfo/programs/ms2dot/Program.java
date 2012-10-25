//package edu.rice.cs.bioinfo.programs.ms2dot;
//
//import java.util.HashMap;
//import java.util.HashSet;
//import java.util.LinkedList;
//import java.util.Map;
//
///**
// * Created by IntelliJ IDEA.
// * User: Matt
// * Date: 9/11/12
// * Time: 4:17 PM
// * To change this template use File | Settings | File Templates.
// */
//public class Program
//{
//    static class Edge
//    {
//        Object Tail;
//        Object Tip;
//        String popNumber;
//    }
//
//    public static void main(String[] args)
//    {
//        String msString = args[0];
//
//        String[] commandParts = msString.split("\\s+");
//
//        HashSet<String> knownPopulations = new HashSet<String>();
//
//        StringBuffer dotResult = new StringBuffer();
//        int nextNodeNum = 0;
//        dotResult.append("digraph graph {");
//
//
//        Map<String,Object> populationToLeadNode = new HashMap<String, Object>();
//        LinkedList<Edge> edges = new LinkedList<Edge>();
//
//        for(int i = 1; i<commandParts; i++)
//        {
//            if(commandParts[i] == "ej")
//            {
//                i++;
//                String joinPop1 = commandParts[i];
//                i++;
//                String joinPop2 = commandParts[i];
//
//                if(!knownPopulations.contains(joinPop1) && !knownPopulations.contains(joinPop2))
//                {
//                    Edge edge1 = new Edge();
//                    edge1.Tail = new Object();
//                    edge1.Tip = new Object();
//                    edge1.popNumber = joinPop1;
//                    edges.add(edge1);
//                    knownPopulations.add(joinPop1);
//
//                    Edge edge2 = new Edge();
//                    edge2.Tail = new Object();
//                    edge2.Tip = edge1.Tip;
//                    edge2.popNumber = joinPop2;
//                    edges.add(edge2);
//                    knownPopulations.add(joinPop2);
//
//                    populationToLeadNode.put(joinPop2, edge2.Tip);
//                }
//                else if(knownPopulations.contains(joinPop1) && !knownPopulations.contains(joinPop2))
//                {
//                    Object popLeadNode = populationToLeadNode.remove(joinPop1);
//
//                    Edge edge = new Edge();
//                    edge.Tail = new Object();
//                    edge.Tip = popLeadNode;
//                    edge.popNumber = joinPop2;
//                    edges.add(edge);
//                    knownPopulations.add(joinPop2);
//
//                    populationToLeadNode.put(joinPop2, popLeadNode);
//
//                }
//                else if(knownPopulations.contains(joinPop1) && knownPopulations.contains(joinPop2))
//                {
//                    Edge edge1 = new Edge();
//                    edge1.Tail = populationToLeadNode.remove(joinPop1);
//                    edge1.Tip = new Object();
//                    edge1.popNumber = joinPop1;
//                    knownPopulations.add(joinPop1);
//
//                    Edge edge2 = new Edge();
//                    edge2.Tail = populationToLeadNode.remove(joinPop2);
//                    edge2.Tip = edge1.Tip;
//                    edge2.popNumber = joinPop2;
//                    knownPopulations.add(joinPop2);
//                }
//            }
//        }
//
//        dotResult.append("}");
//    }
//}
